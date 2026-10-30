###############################################################################
# STAT 638 – Group 6 – Piecewise Exponential Model
# LOG-LIKELIHOOD METRICS + LOG MARGINAL LIKELIHOOD (for Bayes factors)
#
# Run AFTER piecewise_exponential_KIRC_GibbsSampler_updated.R (needs PWE_Fit.rds).
# Uses the SURVIVAL likelihood (same scale as the Weibull/lognormal log_lik):
#   death:    log f(t_i) = log h(t_i) - H(t_i)
#   censored: log S(t_i) = -H(t_i)
#
# Outputs
#   PiecewiseExp_LogLik_Metrics.csv   one-row table (same idea as the metrics CSV)
#   PiecewiseExp_Bridge.rds           bridge-sampling object (log marginal likelihood)
#
# Bayes factor between two models A and B (same data, proper priors):
#   log BF_AB = logML_A - logML_B          (> 0 favours A)
###############################################################################

#install.packages(c("readxl","dplyr","loo","bridgesampling"))
library(readxl)
library(dplyr)
library(loo)
library(bridgesampling)

set.seed(638)
fit <- readRDS("PWE_Fit.rds")          # saved chains from the Gibbs sampler (informative priors)
stopifnot(fit$time_divisor == 365)

# ---------------------------------------------------------------------------
# Data: identical to KIRC_PiecewiseExp_Section10.R (same patients, order, coding)
# ---------------------------------------------------------------------------
dat <- read_excel("KIRC - Filtered Data.xlsx") %>%
  filter(toupper(trimws(histologic_grade)) %in% c("G1", "G2", "G3", "G4")) %>%
  arrange(bcr_patient_barcode) %>%
  mutate(OS_years = as.numeric(overall_survival_time) / 365,
         event    = as.integer(overall_survival),
         Age10    = (as.numeric(age_at_diagnosis) - 60) / 10,
         Male     = as.integer(toupper(trimws(gender)) == "MALE"),
         stage    = sub("Stage ", "", trimws(pathologic_stage)))
N <- nrow(dat); stopifnot(N == 521, sum(dat$event) == 173)

X_pat <- cbind(stageII  = as.numeric(dat$stage == "II"),
               stageIII = as.numeric(dat$stage == "III"),
               stageIV  = as.numeric(dat$stage == "IV"),
               age10    = dat$Age10,
               male     = dat$Male,
               gradeG3  = as.numeric(dat$histologic_grade == "G3"),
               gradeG4  = as.numeric(dat$histologic_grade == "G4"))
beta_names <- colnames(X_pat)
lam_names  <- paste0("lambda", 1:5)

cut_all <- fit$cut_all
lo <- head(cut_all, -1); hi <- tail(cut_all, -1)
E_pat <- t(sapply(dat$OS_years, function(t) pmax(0, pmin(t, hi) - lo)))      # N x 5
k_pat <- findInterval(dat$OS_years, cut_all[2:5], left.open = TRUE) + 1

# Priors used in the fit (saved in the fit file)
a0 <- fit$a0; b0 <- fit$b0; mu0 <- fit$mu0; sd0 <- fit$sd0

# Total log-likelihood for one parameter vector
loglik_total <- function(lam, beta) {
  lp <- drop(X_pat %*% beta)
  sum(dat$event * (log(lam[k_pat]) + lp) - drop(E_pat %*% lam) * exp(lp))
}

# ---------------------------------------------------------------------------
# Posterior draws (all chains, burn-in removed)
# ---------------------------------------------------------------------------
chains <- lapply(fit$fits, function(f) f$draws[-(1:fit$burn), c(lam_names, beta_names)])
post   <- do.call(rbind, chains)
thin   <- post[seq(1, nrow(post), by = 5), ]                  # 16,800 draws

# ---------------------------------------------------------------------------
# 1. Log-likelihood summaries
# ---------------------------------------------------------------------------
ll_draws <- apply(thin, 1, function(p) loglik_total(p[lam_names], p[beta_names]))

ll_mean    <- mean(ll_draws)                                   # posterior mean of log L
ll_at_mean <- loglik_total(colMeans(post[, lam_names]), colMeans(post[, beta_names]))
ll_max     <- max(ll_draws)                                    # highest log L among draws

p_DIC <- 2 * (ll_at_mean - ll_mean)
DIC   <- -2 * ll_at_mean + 2 * p_DIC

# pointwise log-likelihood matrix (draws x patients) for lppd and WAIC
lp_mat  <- thin[, beta_names] %*% t(X_pat)
log_lik <- sweep(log(thin[, lam_names][, k_pat]) + lp_mat, 2, dat$event, "*") -
           (thin[, lam_names] %*% t(E_pat)) * exp(lp_mat)
lppd    <- sum(log(colMeans(exp(log_lik))))
waic_o  <- waic(log_lik)
rm(lp_mat); invisible(gc())

# ---------------------------------------------------------------------------
# 2. Log marginal likelihood by bridge sampling  ->  Bayes factors
#    log p(y) = log of the integral of likelihood x prior over all parameters
# ---------------------------------------------------------------------------
log_posterior <- function(pars, data) {
  lam  <- pars[lam_names]; beta <- pars[beta_names]
  loglik_total(lam, beta) +
    sum(dgamma(lam, shape = a0, rate = b0, log = TRUE)) +
    sum(dnorm(beta, mean = mu0, sd = sd0, log = TRUE))
}

bs_draws <- post[round(seq(1, nrow(post), length.out = 20000)), ]   # 20,000 draws
lb <- setNames(c(rep(0, 5), rep(-Inf, 7)), c(lam_names, beta_names))
ub <- setNames(rep(Inf, 12), c(lam_names, beta_names))

bridge <- bridge_sampler(samples = bs_draws, log_posterior = log_posterior,
                         data = NULL, lb = lb, ub = ub,
                         method = "normal", repetitions = 1, silent = TRUE)
print(bridge)
logML     <- bridge$logml
logML_err <- error_measures(bridge)$percentage
saveRDS(bridge, "PiecewiseExp_Bridge.rds")

# ---------------------------------------------------------------------------
# 3. Save one-row table
# ---------------------------------------------------------------------------
ll_metrics <- data.frame(
  Model = "Piecewise exponential", N = N, N_params = 12,
  LogLik_posterior_mean = ll_mean,        # E[log L | data]
  LogLik_at_post_mean   = ll_at_mean,     # log L(theta-bar)
  LogLik_max_draw       = ll_max,
  lppd = lppd,
  WAIC = waic_o$estimates["waic", "Estimate"], p_WAIC = waic_o$estimates["p_waic", "Estimate"],
  DIC = DIC, p_DIC = p_DIC,
  Log_Marginal_Likelihood = logML, LogML_error = logML_err
)
print(ll_metrics, digits = 6)
write.csv(ll_metrics, "PiecewiseExp_LogLik_Metrics.csv", row.names = FALSE)
# Expected (approx.): mean -529.4, at posterior mean -525.1, DIC ~1067, logML ~ -535.8

# ---------------------------------------------------------------------------
# 4. Bayes factors (once the other models' log marginal likelihoods are available)
# ---------------------------------------------------------------------------
# logML_weibull   <- ...   # from the Weibull model (bridge sampling on its fit)
# logML_lognormal <- ...
# logBF <- c(PWE_vs_Weibull   = logML - logML_weibull,
#            PWE_vs_Lognormal = logML - logML_lognormal)
# exp(logBF)   # BF > 1 favours the piecewise exponential model
