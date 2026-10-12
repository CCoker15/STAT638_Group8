###############################################################################
# STAT 638 – Group 6 – Bayesian Piecewise Exponential Model (KIRC overall survival)
# Covariates: Stage (I–IV), Age, Sex  +  Histologic grade (tumor aggressiveness)
# Output: S(t) for the four stages in one plot – one plot for a 60-yr-old male,
#         one for a 60-yr-old female (grade held at G1–G2, the reference level)
#
# Posterior sampled with a hand-coded Gibbs sampler (Ch. 6 style):
#   - lambda_k  : exact Gibbs draws (Gamma prior is conjugate)
#   - beta_j    : Metropolis step inside the Gibbs loop (no closed-form FCD)
# No Stan / brms / compiler needed.
###############################################################################

#install.packages(c("readxl","dplyr","survival","coda","ggplot2","patchwork"))
library(readxl)
library(dplyr)
library(survival)
library(coda)        # MCMC diagnostics (effectiveSize, gelman.diag, traceplot), as in Ch. 6
library(ggplot2)
library(patchwork)

# ---------------------------------------------------------------------------
# STEP 1. Load and clean the data (time in YEARS)
# ---------------------------------------------------------------------------
raw <- read_excel("KIRC - Filtered Data.xlsx")   # the 529-patient filtered file

kirc <- raw %>%
  # 8 patients have grade GX (5) or missing (3) -> drop them (529 -> 521)
  filter(histologic_grade %in% c("G1", "G2", "G3", "G4")) %>%
  transmute(
    time   = overall_survival_time / 365.25,
    status = overall_survival,                      # 1 = died, 0 = censored
    stage  = factor(sub("Stage ", "", pathologic_stage),
                    levels = c("I", "II", "III", "IV")),
    age10  = (age_at_diagnosis - 60) / 10,          # same centering as Weibull model
    male   = as.integer(gender == "MALE"),
    # G1 has only 13 patients and ZERO deaths -> its coefficient would not be
    # identifiable, so G1 and G2 are combined into one "low grade" reference group
    grade  = factor(ifelse(histologic_grade %in% c("G1", "G2"), "G1&2", histologic_grade),
                    levels = c("G1&2", "G3", "G4"))
  )

nrow(kirc)          # 521
sum(kirc$status)    # 173 deaths
table(kirc$grade, kirc$status)

# ---------------------------------------------------------------------------
# STEP 2. Choose cutpoints (enough deaths in every interval)
# ---------------------------------------------------------------------------
cuts    <- c(1, 2, 3, 5)                 # intervals: 0-1, 1-2, 2-3, 3-5, 5+ years
cut_all <- c(0, cuts, Inf)

table(cut(kirc$time[kirc$status == 1], cut_all, right = FALSE))
# 48 / 32 / 29 / 39 / 25 deaths -> every interval is well populated.
# (Only ~2 deaths after year 8, which is why everything past 5 years is one interval.)

# ---------------------------------------------------------------------------
# STEP 3. Wide -> long with survSplit (one row per patient per interval)
# ---------------------------------------------------------------------------
long <- survSplit(Surv(time, status) ~ ., data = kirc, cut = cuts,
                  start = "tstart", end = "tstop", episode = "interval") %>%
  mutate(exposure = tstop - tstart,              # time at risk inside the interval
         interval = factor(interval))            # levels "1".."5"

head(long[, c("tstart", "tstop", "exposure", "status", "interval")], 10)
nrow(long)              # 1746 rows
sum(long$status)        # must still be 173 deaths

# ---------------------------------------------------------------------------
# STEP 4. Model, priors and full conditionals
# ---------------------------------------------------------------------------
# Row r of the long data = patient i in interval k(r):
#   d_r | lambda, beta ~ Poisson( lambda_k(r) * e_r * exp(x_i' beta) )     ("Poisson trick")
#   d_r = died in this interval (0/1),  e_r = exposure (years at risk in the interval)
#
# Priors:
#   lambda_k ~ Gamma(a0, b0) = Gamma(0.01, 0.01)     k = 1..5  (baseline hazard per year)
#   beta_j   ~ Normal(0, 5^2)                         j = 1..7  (log hazard ratios)
#
# Likelihood (kernel):
#   prod_r [lambda_k e^{x'beta}]^{d_r} exp(-lambda_k e_r e^{x'beta})
#
# FCD of lambda_k (collect the lambda_k terms; same algebra as Poisson-Gamma, Ch. 3):
#   lambda_k | beta, data ~ Gamma( a0 + D_k ,  b0 + sum_{r in k} e_r exp(x_r' beta) )
#   D_k = deaths in interval k; the rate is the "covariate-weighted exposure".
#   -> exactly like Quiz 2 Q1: shape = a + total count, rate = b + total exposure.
#
# FCD of beta_j (not a known distribution):
#   log pi(beta_j | rest) = beta_j * sum_r d_r x_rj
#                           - sum_r lambda_k e_r exp(x_r' beta)
#                           - beta_j^2 / (2 * 5^2)  + const
#   -> sample it with a random-walk Metropolis step inside the Gibbs loop.

a0  <- 0.01; b0 <- 0.01     # Gamma prior on each lambda_k
sd0 <- 5                    # Normal(0, 5^2) prior on each beta_j

# Pieces the sampler needs (computed once)
X    <- model.matrix(~ stage + age10 + male + grade, data = long)[, -1]  # 7 covariate columns
Z    <- model.matrix(~ 0 + interval, data = long)                        # 1746 x 5 interval indicators
kid  <- as.integer(long$interval)   # interval number of each row
e    <- long$exposure
d    <- long$status
K    <- ncol(Z); p <- ncol(X)
D_k  <- colSums(Z * d)              # deaths per interval: 48 32 29 39 25
Xd   <- colSums(X * d)              # sum_r d_r x_rj  (sufficient statistic for each beta_j)
colnames(X)                         # "stageII" "stageIII" "stageIV" "age10" "male" "gradeG3" "gradeG4"

# Metropolis proposal SDs: 2 x the frequentist standard errors of the same model.
# (Only sets the step size; any reasonable value still targets the right posterior.)
freq_pwe <- glm(status ~ 0 + interval + stage + age10 + male + grade + offset(log(exposure)),
                family = poisson, data = long)
step <- 2 * sqrt(diag(vcov(freq_pwe)))[colnames(X)]

# ---------------------------------------------------------------------------
# STEP 5. Metropolis-within-Gibbs sampler
# ---------------------------------------------------------------------------
gibbs_pwe <- function(S, beta_init, seed) {
  set.seed(seed)
  beta  <- beta_init
  eta   <- drop(X %*% beta)                      # x_r' beta for every row
  acc   <- rep(0, p)
  store <- matrix(NA, S, K + p,
                  dimnames = list(NULL, c(paste0("lambda", 1:K), colnames(X))))

  for (s in 1:S) {
    # (1) Gibbs step: all 5 lambda_k from their Gamma FCDs, given current beta
    rate_k <- b0 + drop(crossprod(Z, e * exp(eta)))
    lambda <- rgamma(K, shape = a0 + D_k, rate = rate_k)

    # (2) Metropolis step for each beta_j, given lambda and the other betas
    w <- lambda[kid] * e                          # lambda_k(r) * e_r
    for (j in 1:p) {
      prop  <- beta[j] + rnorm(1, 0, step[j])     # random-walk proposal
      eta_p <- eta + X[, j] * (prop - beta[j])
      log_r <- Xd[j] * (prop - beta[j]) -                     # log-likelihood ratio ...
               sum(w * (exp(eta_p) - exp(eta))) -
               (prop^2 - beta[j]^2) / (2 * sd0^2)             # ... + log-prior ratio
      if (log(runif(1)) < log_r) {                # accept with prob min(1, r)
        beta[j] <- prop; eta <- eta_p; acc[j] <- acc[j] + 1
      }
    }
    store[s, ] <- c(lambda, beta)
  }
  list(draws = store, acc = acc / S)
}

# 3 chains from different random starting values (to check convergence)
S <- 30000
fits <- lapply(1:3, function(ch) gibbs_pwe(S, beta_init = rnorm(p, 0, 1), seed = 638 + ch))
# roughly 20-60 seconds per chain, depending on the computer

# ---------------------------------------------------------------------------
# STEP 6. Diagnostics (Ch. 6: trace plot, ACF, ESS, burn-in, thinning)
# ---------------------------------------------------------------------------
round(sapply(fits, `[[`, "acc"), 2)   # acceptance rates: ~0.3-0.5 is good for these steps

# Look at the start of the chains to choose the burn-in
matplot(sapply(fits, function(f) f$draws[1:3000, "stageIV"]), type = "l", lty = 1,
        ylab = "beta stageIV", xlab = "iteration", main = "First 3000 iterations, 3 chains")

burn  <- 2000
chains <- mcmc.list(lapply(fits, function(f) mcmc(f$draws[-(1:burn), ])))

plot(chains[, c("lambda1", "stageIV", "age10", "gradeG4")])  # trace + density, 3 chains overlapping
effectiveSize(chains)          # ESS summed over chains; want well over 1000 per parameter
gelman.diag(chains)            # R-hat (potential scale reduction); want < 1.01-1.05
acf(chains[[1]][, "stageIV"])  # autocorrelation dies out after ~10-20 lags
autocorr.plot(chains[[1]][, c("stageIV", "gradeG4")])

# Pool the 3 chains after burn-in
post <- as.matrix(chains)      # ~84,000 x 12 posterior draws
# Thin by 5 for the survival curves (draws 5 apart are much less correlated + saves memory)
draws <- post[seq(1, nrow(post), by = 5), ]

# ---------------------------------------------------------------------------
# STEP 7. Hazard ratios with 95% credible intervals
# ---------------------------------------------------------------------------
hr_vars <- colnames(X)
hr_tab <- t(sapply(hr_vars, function(v)
  c(quantile(exp(post[, v]), c(0.5, 0.025, 0.975)),
    "P(HR>1)" = mean(post[, v] > 0))))           # posterior prob the variable raises the hazard
colnames(hr_tab)[1:3] <- c("HR (median)", "2.5%", "97.5%")
round(hr_tab, 3)
# Expected (approx.):
#   stageII 1.17 (0.60, 2.12)   stageIII 2.03 (1.32, 3.09)   stageIV 4.75 (3.07, 7.38)
#   age10   1.37 (1.18, 1.58)   male     0.94 (0.69, 1.28)
#   gradeG3 1.43 (0.97, 2.13)   gradeG4  2.24 (1.43, 3.58)

# Baseline hazards (deaths per person-year for the reference patient)
round(t(apply(post[, paste0("lambda", 1:K)], 2, quantile, c(0.5, 0.025, 0.975))), 4)

# ---------------------------------------------------------------------------
# STEP 8. Posterior survival curves
#   S(t | x) = exp( - exp(x'beta) * sum_k lambda_k * (time spent in interval k up to t) )
#   computed for EVERY posterior draw, then summarised (median + 95% band)  (Ch. 4 idea)
# ---------------------------------------------------------------------------
tgrid <- seq(0, 10, by = 0.05)                       # years

# E[j, k] = time spent in interval k by time tgrid[j]
E <- sapply(1:5, function(k) pmax(0, pmin(tgrid, cut_all[k + 1]) - cut_all[k]))

lambda <- draws[, paste0("lambda", 1:5)]             # draws x 5
H0     <- E %*% t(lambda)                            # cumulative baseline hazard: time x draws

surv_curve <- function(stage, male, age = 60, grade = "G1&2") {
  lp <- (if (stage != "I")    draws[, paste0("stage", stage)] else 0) +
        draws[, "age10"] * (age - 60) / 10 +
        draws[, "male"]  * male +
        (if (grade != "G1&2") draws[, paste0("grade", grade)] else 0)
  S <- exp(-sweep(H0, 2, exp(lp), "*"))              # time x draws
  data.frame(time  = tgrid,
             stage = paste("Stage", stage),
             est   = apply(S, 1, median),
             lo    = apply(S, 1, quantile, 0.025),
             hi    = apply(S, 1, quantile, 0.975))
}

plot_sex <- function(male, title) {
  df <- bind_rows(lapply(c("I", "II", "III", "IV"), surv_curve, male = male))
  ggplot(df, aes(time, est, colour = stage, fill = stage)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = cuts, linetype = "dotted", colour = "grey60") +
    scale_x_continuous(breaks = 0:10) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = title, x = "Years since diagnosis", y = "S(t)",
         colour = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

p_male   <- plot_sex(1, "Male, age 60 (grade G1–G2)")
p_female <- plot_sex(0, "Female, age 60 (grade G1–G2)")

p_male + p_female
ggsave("pwe_survival_male60.png",   p_male,   width = 6, height = 5, dpi = 300)
ggsave("pwe_survival_female60.png", p_female, width = 6, height = 5, dpi = 300)

# 5-year survival table for the meeting
i5 <- which.min(abs(tgrid - 5))
s5 <- bind_rows(
  lapply(c("I","II","III","IV"), function(s) cbind(sex = "Male",   surv_curve(s, 1)[i5, ])),
  lapply(c("I","II","III","IV"), function(s) cbind(sex = "Female", surv_curve(s, 0)[i5, ]))
)
s5
# Expected (approx.), median (95% CrI):
#   Male   I 0.85 (0.79, 0.89)  II 0.82 (0.71, 0.91)  III 0.71 (0.59, 0.81)  IV 0.45 (0.30, 0.61)
#   Female I 0.84 (0.77, 0.89)  II 0.81 (0.68, 0.90)  III 0.70 (0.57, 0.80)  IV 0.43 (0.26, 0.60)

# ---------------------------------------------------------------------------
# STEP 9. Frequentist check (same likelihood, no priors) + Cox comparison
# ---------------------------------------------------------------------------
round(exp(coef(freq_pwe)), 3)        # should be close to the Bayesian medians in Step 7

cox <- coxph(Surv(time, status) ~ stage + age10 + male + grade, data = kirc)
summary(cox)                         # HRs should be close to the PWE HRs
cox.zph(cox)                         # proportional-hazards check (PWE makes the same assumption)
