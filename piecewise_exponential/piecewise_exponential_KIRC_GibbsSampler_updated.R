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
#
# Priors: shared with the group's Weibull model (same literature hazard ratios and
# the same calibrated baseline survival for the reference patient), translated to
# the piecewise-exponential (log-hazard) scale; derived in PWE_Prior_Derivation.R.
# Set prior_type <- "vague" in Step 4 to rerun with the original vague priors
# as a sensitivity analysis.
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
    time   = overall_survival_time / 365,
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
# Priors (shared with the Weibull model; values derived in PWE_Prior_Derivation.R):
#   lambda_k ~ Gamma(a0_k, b0_k)          k = 1..5  baseline hazard per year, calibrated so the
#                                                    reference patient's S(1), S(5), S(10) are plausible
#   beta_j   ~ Normal(mu0_j, sd0_j^2)     j = 1..7  literature-based log hazard ratios
#
# Likelihood (kernel):
#   prod_r [lambda_k e^{x'beta}]^{d_r} exp(-lambda_k e_r e^{x'beta})
#
# FCD of lambda_k (collect the lambda_k terms; same algebra as Poisson-Gamma, Ch. 3):
#   lambda_k | beta, data ~ Gamma( a0_k + D_k ,  b0_k + sum_{r in k} e_r exp(x_r' beta) )
#   D_k = deaths in interval k; the rate is the "covariate-weighted exposure".
#   -> exactly like Quiz 2 Q1: shape = a + total count, rate = b + total exposure.
#
# FCD of beta_j (not a known distribution):
#   log pi(beta_j | rest) = beta_j * sum_r d_r x_rj
#                           - sum_r lambda_k e_r exp(x_r' beta)
#                           - (beta_j - mu0_j)^2 / (2 * sd0_j^2)  + const
#   -> sample it with a random-walk Metropolis step inside the Gibbs loop.

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

# ---------------------------------------------------------------------------
# PRIOR VALUES (copied from the output of PWE_Prior_Derivation.R)
#   "informative" = shared literature-based priors (main analysis, comparable to Weibull)
#   "vague"       = original weak priors (prior sensitivity analysis)
# ---------------------------------------------------------------------------
prior_type <- "informative"

if (prior_type == "informative") {
  # Baseline hazards lambda_k ~ Gamma(shape a0_k, rate b0_k): 0-1, 1-2, 2-3, 3-5, 5+ yrs
  a0 <- c(2.74, 5.99, 8.52, 11.63, 14.94)
  b0 <- c(112.8, 208.3, 272.6, 342.8, 401.4)
  # Coefficients beta_j = log(HR) ~ Normal(mu0_j, sd0_j^2)
  mu0 <- c(stageII  = 0.405,   # HR 1.5  (0.8-2.8)
           stageIII = 0.916,   # HR 2.5  (1.3-4.8)
           stageIV  = 2.303,   # HR 10   (4-25)
           age10    = 0.392,   # HR 1.04 per year (1.02-1.06), x10
           male     = 0.000,   # HR 1    (0.67-1.5)
           gradeG3  = 0.336,   # HR 1.4  (0.9-2.2)
           gradeG4  = 0.916)   # HR 2.5  (1.4-4.5)
  sd0 <- c(0.320, 0.333, 0.467, 0.098, 0.206, 0.228, 0.298)
} else {
  a0  <- rep(0.01, 5); b0 <- rep(0.01, 5)
  mu0 <- setNames(rep(0, 7), colnames(X)); sd0 <- rep(5, 7)
}
stopifnot(identical(names(mu0), colnames(X)))   # priors line up with the X columns

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
               ((prop - mu0[j])^2 - (beta[j] - mu0[j])^2) / (2 * sd0[j]^2)   # ... + log-prior ratio
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

saveRDS(list(fits = fits, burn = burn, cut_all = cut_all, time_divisor = 365,
             prior_type = prior_type, a0 = a0, b0 = b0, mu0 = mu0, sd0 = sd0),
        paste0("PWE_Fit_", prior_type, ".rds"))

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
# Expected (approx.), informative priors:
#   stageII 1.37 (0.89, 2.08)   stageIII 2.26 (1.64, 3.11)   stageIV 5.45 (3.86, 7.67)
#   age10   1.41 (1.26, 1.58)   male     0.95 (0.75, 1.20)
#   gradeG3 1.39 (1.06, 1.81)   gradeG4  2.19 (1.58, 3.02)
# Vague priors (sensitivity run): II 1.16, III 2.01, IV 4.69, age 1.37, male 0.93, G3 1.44, G4 2.26

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

(p_male + p_female) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
#ggsave("pwe_survival_both.png", width = 11, height = 5, dpi = 300)
#ggsave("pwe_survival_male60.png",   p_male,   width = 6, height = 5, dpi = 300)
#ggsave("pwe_survival_female60.png", p_female, width = 6, height = 5, dpi = 300)

# 5-year survival table for the meeting
i5 <- which.min(abs(tgrid - 5))
s5 <- bind_rows(
  lapply(c("I","II","III","IV"), function(s) cbind(sex = "Male",   surv_curve(s, 1)[i5, ])),
  lapply(c("I","II","III","IV"), function(s) cbind(sex = "Female", surv_curve(s, 0)[i5, ]))
)
s5
# Expected medians, informative priors: Male I-IV ~0.86/0.81/0.71/0.43, Female ~0.85/0.80/0.69/0.41
# Vague priors (sensitivity run):        Male I-IV  0.85/0.82/0.71/0.46, Female  0.84/0.81/0.70/0.43

# ---------------------------------------------------------------------------
# STEP 8b. Posterior survival at 1, 5 and 10 years (age 60, grade G1-G2)
# ---------------------------------------------------------------------------
times_interest <- c(1, 5, 10)
lambda_post    <- post[, paste0("lambda", 1:5)]          # all 84,000 draws x 5

# S(t) for every posterior draw at one time point
surv_draws <- function(stage, male, t, age = 60, grade = "G1&2") {
  lp <- (if (stage != "I")    post[, paste0("stage", stage)] else 0) +
    post[, "age10"] * (age - 60) / 10 +
    post[, "male"]  * male +
    (if (grade != "G1&2") post[, paste0("grade", grade)] else 0)
  Et <- pmax(0, pmin(t, cut_all[-1]) - cut_all[-length(cut_all)])   # years in each interval by t
  exp(-drop(lambda_post %*% Et) * exp(lp))
}

surv_table <- do.call(rbind, lapply(c(0, 1), function(m)
  do.call(rbind, lapply(c("I", "II", "III", "IV"), function(s)
    do.call(rbind, lapply(times_interest, function(t) {
      S <- surv_draws(s, m, t)
      data.frame(Sex              = ifelse(m == 1, "Male", "Female"),
                 Stage            = paste("Stage", s),
                 Time_years       = t,
                 Posterior_mean   = mean(S),
                 Posterior_median = median(S),
                 Lower_95_CrI     = quantile(S, 0.025, names = FALSE),
                 Upper_95_CrI     = quantile(S, 0.975, names = FALSE))
    }))))))

surv_table <- surv_table %>%
  mutate(Median_percent   = 100 * Posterior_median,
         Lower_95_percent = 100 * Lower_95_CrI,
         Upper_95_percent = 100 * Upper_95_CrI)

print(surv_table, digits = 3, row.names = FALSE)
#write.csv(surv_table,
#          paste0("PWE_Survival_1_5_10_age60_", prior_type, ".csv"),
#          row.names = FALSE)


# ---------------------------------------------------------------------------
# STEP 9. Frequentist check (same likelihood, no priors) + Cox comparison
# ---------------------------------------------------------------------------
round(exp(coef(freq_pwe)), 3)        # should be close to the Bayesian medians in Step 7

cox <- coxph(Surv(time, status) ~ stage + age10 + male + grade, data = kirc)
summary(cox)                         # HRs should be close to the PWE HRs
cox.zph(cox)                         # proportional-hazards check (PWE makes the same assumption)
