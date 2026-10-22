# ============================================================
# WEIBULL MODEL COMPARISON METRICS - KIRC
#
# Run AFTER:
# KIRC_Weibull_Posterior_Stage_Age_Grade.R
#
# Reports:
#   ELPD-LOO and SE
#   LOOIC
#   Pareto-k diagnostics
#   Brier score at 1, 5, 10 years
#   Integrated Brier Score (0-10 years)
#   max R-hat, min ESS, divergences
#   log marginal likelihood for Bayes-factor comparison
#   bridge-sampling stability diagnostics
#
# PRIMARY MODEL-SELECTION METRIC:
#   larger ELPD-LOO = better
#
# SUPPORTING METRICS:
#   smaller Brier / IBS = better
# ============================================================

library(readxl)
library(dplyr)
library(survival)
library(cmdstanr)
library(posterior)
library(loo)

if (!requireNamespace("bridgesampling", quietly = TRUE)) {
  stop(
    "Package 'bridgesampling' is required. Install it with:\n",
    "install.packages('bridgesampling')"
  )
}

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_dir <- "D:/OneDrive - Texas A&M University/PhD/Courses/STAT 638/Project"
output_dir  <- file.path(project_dir, "Stage_3_Grade")

data_file <- file.path(project_dir, "KIRC - Filtered Data.xlsx")
fit_file  <- file.path(
  output_dir,
  "KIRC_Weibull_Stage_Age_Sex_Grade_Fit.rds"
)

if (!file.exists(data_file)) stop("Data file not found.")
if (!file.exists(fit_file))  stop("Saved Weibull fit not found.")

# ------------------------------------------------------------
# 2. Recreate EXACT model dataset
# ------------------------------------------------------------

dat0 <- read_excel(data_file)

grade_raw <- toupper(trimws(as.character(dat0$histologic_grade)))

dat0$GradeGroup <- case_when(
  grade_raw %in% c("G1", "G2") ~ "G1-G2",
  grade_raw == "G3" ~ "G3",
  grade_raw == "G4" ~ "G4",
  TRUE ~ NA_character_
)

# Same exclusion used in the fitted model:
# remove GX and missing/unusable grade
dat <- dat0 %>%
  filter(!is.na(GradeGroup))

dat$OS_years <- as.numeric(dat$overall_survival_time) / 365
dat$event    <- as.integer(dat$overall_survival)

dat$Age_years <- as.numeric(dat$age_at_diagnosis)
dat$Age10     <- (dat$Age_years - 60) / 10

sex_text <- toupper(trimws(as.character(dat$gender)))
dat$Male <- as.integer(sex_text == "MALE")

stage_text <- trimws(as.character(dat$pathologic_stage))
dat$StageII  <- as.integer(stage_text == "Stage II")
dat$StageIII <- as.integer(stage_text == "Stage III")
dat$StageIV  <- as.integer(stage_text == "Stage IV")

dat$GradeG3 <- as.integer(dat$GradeGroup == "G3")
dat$GradeG4 <- as.integer(dat$GradeGroup == "G4")

keep <- complete.cases(
  dat$OS_years, dat$event,
  dat$StageII, dat$StageIII, dat$StageIV,
  dat$Age10, dat$Male,
  dat$GradeG3, dat$GradeG4
)

dat <- dat[keep, , drop = FALSE]
N <- nrow(dat)

cat("\nFinal comparison N =", N, "\n")

# ------------------------------------------------------------
# 3. Load fitted Weibull model
# ------------------------------------------------------------

fit <- readRDS(fit_file)

par_names <- c(
  "beta0", "beta2", "beta3", "beta4",
  "betaA", "betaS",
  "gamma3", "gamma4",
  "log_k", "k"
)

# ------------------------------------------------------------
# 4. MCMC diagnostics
# ------------------------------------------------------------

sum_par <- fit$summary(variables = par_names)

max_rhat    <- max(sum_par$rhat, na.rm = TRUE)
min_ess_bulk <- min(sum_par$ess_bulk, na.rm = TRUE)
min_ess_tail <- min(sum_par$ess_tail, na.rm = TRUE)

diag_df <- posterior::as_draws_df(fit$sampler_diagnostics())

n_divergent <- sum(diag_df$divergent__)
n_treedepth <- sum(diag_df$treedepth__ >= 12)

cat("\nMCMC diagnostics\n")
cat("Max R-hat       =", max_rhat, "\n")
cat("Min bulk ESS    =", min_ess_bulk, "\n")
cat("Min tail ESS    =", min_ess_tail, "\n")
cat("Divergences     =", n_divergent, "\n")
cat("Treedepth hits  =", n_treedepth, "\n")

# ------------------------------------------------------------
# 5. PSIS-LOO
#
# This is the PRIMARY model comparison.
# Your Stan model already stored pointwise log_lik[i].
# ------------------------------------------------------------

log_lik_array <- fit$draws(
  variables = "log_lik",
  format = "draws_array"
)

if (dim(log_lik_array)[3] != N) {
  stop(
    "log_lik contains ", dim(log_lik_array)[3],
    " patients but recreated data contains ", N,
    ". All models must use exactly the same patients/order."
  )
}

r_eff <- loo::relative_eff(exp(log_lik_array))

loo_weibull <- loo::loo(
  log_lik_array,
  r_eff = r_eff,
  cores = min(4, parallel::detectCores())
)

print(loo_weibull)

saveRDS(
  loo_weibull,
  file.path(output_dir, "LOO_Weibull.rds")
)

elpd_loo <- loo_weibull$estimates["elpd_loo", "Estimate"]
se_elpd  <- loo_weibull$estimates["elpd_loo", "SE"]
p_loo    <- loo_weibull$estimates["p_loo", "Estimate"]
looic    <- loo_weibull$estimates["looic", "Estimate"]
se_looic <- loo_weibull$estimates["looic", "SE"]

# Pareto-k diagnostic
pk <- loo::pareto_k_values(loo_weibull)

max_pk   <- max(pk, na.rm = TRUE)
n_pk_07  <- sum(pk > 0.7, na.rm = TRUE)
n_pk_10  <- sum(pk > 1.0, na.rm = TRUE)

cat("\nPSIS-LOO\n")
cat("ELPD-LOO        =", elpd_loo, "\n")
cat("SE(ELPD-LOO)    =", se_elpd, "\n")
cat("LOOIC           =", looic, "\n")
cat("Max Pareto-k    =", max_pk, "\n")
cat("Pareto-k > 0.7  =", n_pk_07, "\n")
cat("Pareto-k > 1.0  =", n_pk_10, "\n")

# ------------------------------------------------------------
# 6. Posterior draws for survival predictions
# ------------------------------------------------------------

post <- as.data.frame(
  posterior::as_draws_df(
    fit$draws(variables = par_names)
  )
)

M <- nrow(post)

# ------------------------------------------------------------
# 6A. Log marginal likelihood for Bayes-factor comparison
#
# IMPORTANT:
# This uses the MATRIX interface to bridgesampling rather than
# bridge_sampler(fit). Therefore it does NOT require the old
# CmdStan JSON file or CmdStan model-method recompilation.
#
# Quantity estimated:
#
#   log p(y | M_Weibull)
#
# Each teammate should estimate the same quantity for their own
# model using the same patients and their fitted priors.
#
# The priors below are the SAME priors used for this final
# Weibull Stage + Age + Sex + Grade model.
# ------------------------------------------------------------

bridge_par_names <- c(
  "beta0",
  "beta2",
  "beta3",
  "beta4",
  "betaA",
  "betaS",
  "gamma3",
  "gamma4",
  "log_k"
)

# Extract ONLY the sampled parameters.
# Do not include k because k = exp(log_k) is a transformed
# parameter and is not an independent posterior parameter.
bridge_samples <- fit$draws(
  variables = bridge_par_names,
  format = "matrix"
)

bridge_samples <- as.matrix(bridge_samples)

if (!all(bridge_par_names %in% colnames(bridge_samples))) {
  stop(
    "Bridge-sampling posterior matrix is missing parameter(s): ",
    paste(
      setdiff(bridge_par_names, colnames(bridge_samples)),
      collapse = ", "
    )
  )
}

# Preserve the exact parameter order.
bridge_samples <- bridge_samples[, bridge_par_names, drop = FALSE]

cat(
  "\nBridge-sampling posterior draws =",
  nrow(bridge_samples),
  "\n"
)

# Exact dataset used by the log-posterior function.
bridge_data <- list(
  time = as.numeric(dat$OS_years),
  event = as.integer(dat$event),
  StageII = as.numeric(dat$StageII),
  StageIII = as.numeric(dat$StageIII),
  StageIV = as.numeric(dat$StageIV),
  Age10 = as.numeric(dat$Age10),
  Male = as.numeric(dat$Male),
  GradeG3 = as.numeric(dat$GradeG3),
  GradeG4 = as.numeric(dat$GradeG4)
)

if (any(!is.finite(bridge_data$time)) || any(bridge_data$time < 0)) {
  stop("Invalid survival time encountered in bridge-sampling data.")
}

# ------------------------------------------------------------
# Normalized log posterior
#
# The marginal likelihood requires NORMALIZED prior densities.
# dnorm(..., log = TRUE) includes the normalizing constants.
#
# Weibull AFT model:
#
#   eta_i =
#     beta0
#     + beta2 StageII_i
#     + beta3 StageIII_i
#     + beta4 StageIV_i
#     + betaA Age10_i
#     + betaS Male_i
#     + gamma3 GradeG3_i
#     + gamma4 GradeG4_i
#
#   lambda_i = exp(eta_i)
#   k        = exp(log_k)
#
# Death:
#   Weibull density
#
# Censored:
#   Weibull survival probability
# ------------------------------------------------------------

log_posterior_weibull <- function(theta, data) {

  # bridge_sampler may pass an UNNAMED numeric vector.
  # bridge_samples was explicitly ordered as bridge_par_names,
  # so use fixed positions rather than name lookup.
  theta <- as.numeric(theta)

  if (length(theta) != 9 || any(!is.finite(theta))) {
    return(-Inf)
  }

  beta0  <- theta[1]
  beta2  <- theta[2]
  beta3  <- theta[3]
  beta4  <- theta[4]
  betaA  <- theta[5]
  betaS  <- theta[6]
  gamma3 <- theta[7]
  gamma4 <- theta[8]
  log_k  <- theta[9]

  k <- exp(log_k)

  if (!is.finite(k) || k <= 0) {
    return(-Inf)
  }

  eta <-
    beta0 +
    beta2  * data$StageII +
    beta3  * data$StageIII +
    beta4  * data$StageIV +
    betaA  * data$Age10 +
    betaS  * data$Male +
    gamma3 * data$GradeG3 +
    gamma4 * data$GradeG4

  if (any(!is.finite(eta))) {
    return(-Inf)
  }

  # --------------------------
  # Normalized prior density
  # --------------------------
  log_prior <-
    dnorm(beta0,  mean =  3.219, sd = 0.200, log = TRUE) +
    dnorm(log_k,  mean =  0.182, sd = 0.150, log = TRUE) +
    dnorm(beta2,  mean = -0.343, sd = 0.278, log = TRUE) +
    dnorm(beta3,  mean = -0.773, sd = 0.307, log = TRUE) +
    dnorm(beta4,  mean = -1.939, sd = 0.496, log = TRUE) +
    dnorm(betaA,  mean = -0.330, sd = 0.098, log = TRUE) +
    dnorm(betaS,  mean =  0.000, sd = 0.175, log = TRUE) +
    dnorm(gamma3, mean = -0.283, sd = 0.198, log = TRUE) +
    dnorm(gamma4, mean = -0.773, sd = 0.279, log = TRUE)

  if (!is.finite(log_prior)) {
    return(-Inf)
  }

  if (any(data$time <= 0) || any(!is.finite(data$time))) {
    return(-Inf)
  }

  # Stable Weibull AFT likelihood.
  # z = log cumulative hazard = k * (log(t) - eta)
  log_time <- log(data$time)
  log_t_over_lambda <- log_time - eta
  z <- k * log_t_over_lambda

  if (any(!is.finite(z))) {
    return(-Inf)
  }

  # If z exceeds log(max double), the cumulative hazard is
  # effectively infinite and the true log posterior is -Inf.
  if (any(z > log(.Machine$double.xmax))) {
    return(-Inf)
  }

  H <- exp(z)
  died <- data$event == 1L

  ll <- -H

  if (any(died)) {
    ll[died] <-
      log_k - eta[died] +
      (k - 1) * log_t_over_lambda[died] -
      H[died]
  }

  out <- log_prior + sum(ll)

  if (!is.finite(out)) {
    return(-Inf)
  }

  as.numeric(out)
}

# All nine sampled parameters are unconstrained real values.
bridge_lb <- rep(-Inf, length(bridge_par_names))
bridge_ub <- rep( Inf, length(bridge_par_names))

names(bridge_lb) <- bridge_par_names
names(bridge_ub) <- bridge_par_names

# ------------------------------------------------------------
# VALIDATE THE CUSTOM LIKELIHOOD AGAINST STAN'S SAVED log_lik
# ------------------------------------------------------------
# This is a critical check.  At an actual posterior draw, the
# custom R likelihood used by bridge sampling should reproduce
# the pointwise log_lik values generated by Stan.

if (any(bridge_data$time <= 0)) {
  zero_rows <- which(bridge_data$time <= 0)
  cat("\nRows with non-positive survival time:\n")
  print(zero_rows)
  stop(
    "Bridge sampling cannot proceed because the recreated model ",
    "data contain non-positive survival times. This does not match ",
    "the final Weibull fitting script, which required OS_years > 0."
  )
}

stored_log_lik_matrix <- fit$draws(
  variables = "log_lik",
  format = "matrix"
)

if (nrow(stored_log_lik_matrix) != nrow(bridge_samples)) {
  stop("Posterior-draw order mismatch between parameters and log_lik.")
}

weibull_pointwise_loglik <- function(theta, data) {

  # Use the SAME fixed parameter order as bridge_samples.
  theta <- as.numeric(theta)

  if (length(theta) != 9 || any(!is.finite(theta))) {
    return(rep(-Inf, length(data$time)))
  }

  beta0  <- theta[1]
  beta2  <- theta[2]
  beta3  <- theta[3]
  beta4  <- theta[4]
  betaA  <- theta[5]
  betaS  <- theta[6]
  gamma3 <- theta[7]
  gamma4 <- theta[8]
  log_k  <- theta[9]

  k <- exp(log_k)

  eta <-
    beta0 +
    beta2  * data$StageII +
    beta3  * data$StageIII +
    beta4  * data$StageIV +
    betaA  * data$Age10 +
    betaS  * data$Male +
    gamma3 * data$GradeG3 +
    gamma4 * data$GradeG4

  log_t_over_lambda <- log(data$time) - eta
  z <- k * log_t_over_lambda

  # Do not mask all 521 observations if only one contribution is
  # numerically invalid; flag only the affected observations.
  ll <- rep(-Inf, length(data$time))
  good <- is.finite(z) & z <= log(.Machine$double.xmax)

  if (any(good)) {
    H <- exp(z[good])
    good_idx <- which(good)
    ll[good_idx] <- -H

    died_good <- data$event[good_idx] == 1L

    if (any(died_good)) {
      idx <- good_idx[died_good]
      ll[idx] <-
        log_k - eta[idx] +
        (k - 1) * log_t_over_lambda[idx] -
        exp(z[idx])
    }
  }

  as.numeric(ll)
}

# Diagnostic printout for the first posterior draw.
cat("\nBridge parameter order used:\n")
print(bridge_par_names)
cat("\nFirst posterior draw used for bridge validation:\n")
print(setNames(as.numeric(bridge_samples[1, ]), bridge_par_names))

custom_ll_1 <- weibull_pointwise_loglik(
  bridge_samples[1, ],
  bridge_data
)

stan_ll_1 <- as.numeric(stored_log_lik_matrix[1, ])

ll_diff_1 <- custom_ll_1 - stan_ll_1

cat("\nStan log-likelihood sum, draw 1   =", sum(stan_ll_1), "\n")
cat("Custom log-likelihood sum, draw 1 =", sum(custom_ll_1), "\n")
cat("Maximum absolute pointwise difference =",
    max(abs(ll_diff_1)), "\n")

if (
  any(!is.finite(custom_ll_1)) ||
  max(abs(ll_diff_1)) > 1e-6
) {

  bad <- which(
    !is.finite(custom_ll_1) |
    !is.finite(stan_ll_1) |
    abs(ll_diff_1) > 1e-6
  )

  cat("\nFirst rows where custom likelihood differs from Stan:\n")
  print(head(bad, 20))

  if (length(bad) > 0) {
    diagnostic_table <- data.frame(
      Row = bad,
      Time_years = bridge_data$time[bad],
      Event = bridge_data$event[bad],
      Stan_log_lik = stan_ll_1[bad],
      Custom_log_lik = custom_ll_1[bad],
      Difference = ll_diff_1[bad]
    )
    print(head(diagnostic_table, 20))
  }

  stop(
    "The custom Weibull likelihood does not reproduce Stan's saved ",
    "log_lik. Do not compute a Bayes factor until this matches."
  )
}

cat("Custom Weibull likelihood matches Stan log_lik for draw 1.\n")

# Also test several posterior draws, not just the first one.
test_draw_ids <- unique(round(seq(
  1, nrow(bridge_samples), length.out = min(10, nrow(bridge_samples))
)))

for (ii in test_draw_ids) {
  custom_ll_i <- weibull_pointwise_loglik(
    bridge_samples[ii, ], bridge_data
  )
  stan_ll_i <- as.numeric(stored_log_lik_matrix[ii, ])
  max_diff_i <- max(abs(custom_ll_i - stan_ll_i))

  if (!is.finite(max_diff_i) || max_diff_i > 1e-6) {
    stop(
      "Custom likelihood failed validation at posterior draw ", ii,
      ". Maximum pointwise difference = ", max_diff_i
    )
  }
}

cat("Custom Weibull likelihood also matches Stan at ",
    length(test_draw_ids), " posterior draws.\n", sep = "")

# Verify the complete normalized log posterior is finite at an
# actual posterior draw.
bridge_test_lp <- log_posterior_weibull(
  bridge_samples[1, ],
  bridge_data
)

cat(
  "\nBridge log-posterior test value =",
  bridge_test_lp,
  "\n"
)

if (!is.finite(bridge_test_lp)) {
  stop(
    "The stable custom bridge log posterior is still not finite at ",
    "an observed posterior draw. The diagnostic comparison above ",
    "should identify the mismatch."
  )
}

# Save priors used, so all team members can document the
# Bayes-factor calculation consistently.
weibull_prior_table <- data.frame(
  Parameter = c(
    "beta0", "log_k",
    "beta2", "beta3", "beta4",
    "betaA", "betaS",
    "gamma3", "gamma4"
  ),
  Prior_Mean = c(
     3.219,  0.182,
    -0.343, -0.773, -1.939,
    -0.330,  0.000,
    -0.283, -0.773
  ),
  Prior_SD = c(
    0.200, 0.150,
    0.278, 0.307, 0.496,
    0.098, 0.175,
    0.198, 0.279
  )
)

write.csv(
  weibull_prior_table,
  file.path(output_dir, "Weibull_Priors_For_Bayes_Factor.csv"),
  row.names = FALSE
)

cat("\n============================================\n")
cat("BRIDGE SAMPLING: WEIBULL LOG MARGINAL LIKELIHOOD\n")
cat("============================================\n")

set.seed(638)

bridge_weibull <- bridgesampling::bridge_sampler(
  samples = bridge_samples,
  log_posterior = log_posterior_weibull,
  data = bridge_data,
  lb = bridge_lb,
  ub = bridge_ub,
  repetitions = 5,
  method = "warp3",
  cores = 1,
  use_neff = TRUE,
  silent = FALSE
)

print(bridge_weibull)

# Repeated warp3 estimates.
bridge_logml_values <- as.numeric(bridge_weibull$logml)

if (!any(is.finite(bridge_logml_values))) {
  stop("Bridge sampling did not produce any finite log marginal likelihood.")
}

log_marginal_likelihood <- median(
  bridge_logml_values,
  na.rm = TRUE
)

bridge_logml_mean <- mean(
  bridge_logml_values,
  na.rm = TRUE
)

bridge_logml_sd <- sd(
  bridge_logml_values,
  na.rm = TRUE
)

bridge_logml_min <- min(
  bridge_logml_values,
  na.rm = TRUE
)

bridge_logml_max <- max(
  bridge_logml_values,
  na.rm = TRUE
)

bridge_logml_iqr <- IQR(
  bridge_logml_values,
  na.rm = TRUE
)

# Newer bridgesampling versions also report MCSE on the log scale.
if (!is.null(bridge_weibull$mcse_logml)) {
  bridge_mcse_median <- median(
    as.numeric(bridge_weibull$mcse_logml),
    na.rm = TRUE
  )
} else {
  bridge_mcse_median <- NA_real_
}

cat("\nBayes-factor component\n")
cat(
  "Median log p(y | M_Weibull) =",
  log_marginal_likelihood,
  "\n"
)
cat(
  "Mean log p(y | M_Weibull)   =",
  bridge_logml_mean,
  "\n"
)
cat(
  "SD across repetitions       =",
  bridge_logml_sd,
  "\n"
)
cat(
  "IQR across repetitions      =",
  bridge_logml_iqr,
  "\n"
)
cat(
  "Median bridge MCSE          =",
  bridge_mcse_median,
  "\n"
)

# Save the full object and repetition-level estimates.
saveRDS(
  bridge_weibull,
  file.path(output_dir, "BridgeSampling_Weibull.rds")
)

bridge_repetition_table <- data.frame(
  Repetition = seq_along(bridge_logml_values),
  Log_Marginal_Likelihood = bridge_logml_values
)

if (!is.null(bridge_weibull$mcse_logml)) {
  bridge_repetition_table$MCSE_LogML <-
    as.numeric(bridge_weibull$mcse_logml)
}

write.csv(
  bridge_repetition_table,
  file.path(output_dir, "Weibull_Bridge_LogML_Repetitions.csv"),
  row.names = FALSE
)

# eta matrix:
# rows = posterior draws
# columns = patients
eta <- outer(post$beta0, rep(1, N)) +
  outer(post$beta2,  dat$StageII) +
  outer(post$beta3,  dat$StageIII) +
  outer(post$beta4,  dat$StageIV) +
  outer(post$betaA,  dat$Age10) +
  outer(post$betaS,  dat$Male) +
  outer(post$gamma3, dat$GradeG3) +
  outer(post$gamma4, dat$GradeG4)

# Posterior mean survival probability for every patient at time t
posterior_mean_survival <- function(t) {

  if (t == 0) return(rep(1, N))

  # log[(t/lambda)^k] = k * (log(t) - eta)
  z <- sweep(
    log(t) - eta,
    MARGIN = 1,
    STATS = post$k,
    FUN = "*"
  )

  # protect against numerical overflow
  z <- pmin(z, 700)

  S_draws <- exp(-exp(z))

  colMeans(S_draws)
}

# ------------------------------------------------------------
# 7. IPCW Brier score
#
# Lower = better.
#
# IMPORTANT:
# These are supporting apparent-prediction metrics.
# PSIS-LOO remains the primary out-of-sample comparison.
# ------------------------------------------------------------

# KM estimate of censoring survival G(t)
cens_fit <- survfit(
  Surv(OS_years, 1 - event) ~ 1,
  data = dat
)

# G(t), right-continuous
G_right <- function(x) {
  vapply(
    x,
    function(z) {
      idx <- which(cens_fit$time <= z)
      if (length(idx) == 0) return(1)
      cens_fit$surv[max(idx)]
    },
    numeric(1)
  )
}

# G(t-), left limit
G_left <- function(x) {
  vapply(
    x,
    function(z) {
      idx <- which(cens_fit$time < z)
      if (length(idx) == 0) return(1)
      cens_fit$surv[max(idx)]
    },
    numeric(1)
  )
}

G_Tminus <- G_left(dat$OS_years)

ipcw_brier <- function(t, S_hat) {

  G_t <- G_right(t)

  if (G_t <= 1e-8) {
    warning("G(t) is near zero at t = ", t)
  }

  contribution <- rep(0, N)

  # observed death by t
  died <- dat$OS_years <= t & dat$event == 1

  contribution[died] <-
    S_hat[died]^2 /
    pmax(G_Tminus[died], 1e-8)

  # known alive beyond t
  alive_beyond <- dat$OS_years > t

  contribution[alive_beyond] <-
    (1 - S_hat[alive_beyond])^2 /
    pmax(G_t, 1e-8)

  # censored before/equal t contributes 0
  mean(contribution)
}

# ------------------------------------------------------------
# 8. Brier scores at 1, 5, 10 years
# ------------------------------------------------------------

times_report <- c(1, 5, 10)

brier_report <- sapply(
  times_report,
  function(t) {
    S_hat <- posterior_mean_survival(t)
    ipcw_brier(t, S_hat)
  }
)

names(brier_report) <- paste0("Brier_", times_report, "yr")

cat("\nBrier scores\n")
print(brier_report)

# ------------------------------------------------------------
# 9. Integrated Brier Score, 0-10 years
#
# All teammates should use EXACTLY this same time grid.
# ------------------------------------------------------------

ibs_grid <- seq(
  0,
  10,
  by = 0.25
)

bs_grid <- numeric(length(ibs_grid))

for (j in seq_along(ibs_grid)) {

  t <- ibs_grid[j]

  cat(
    "Calculating Brier score at t =",
    sprintf("%.2f", t),
    "\n"
  )

  S_hat <- posterior_mean_survival(t)

  bs_grid[j] <- ipcw_brier(
    t = t,
    S_hat = S_hat
  )
}

brier_curve <- data.frame(
  Time_years = ibs_grid,
  Brier = bs_grid
)

# trapezoidal numerical integration
IBS_0_10 <- sum(
  diff(ibs_grid) *
    (head(bs_grid, -1) + tail(bs_grid, -1)) / 2
) / 10

cat("\nIBS 0-10 years =", IBS_0_10, "\n")

write.csv(
  brier_curve,
  file.path(output_dir, "Weibull_Brier_Curve.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. One-row table for your team comparison
# ------------------------------------------------------------

metrics <- data.frame(
  Model = "Weibull",
  N = N,

  ELPD_LOO = elpd_loo,
  SE_ELPD_LOO = se_elpd,
  p_LOO = p_loo,
  LOOIC = looic,
  SE_LOOIC = se_looic,

  Max_Pareto_k = max_pk,
  N_Pareto_k_over_0_7 = n_pk_07,
  N_Pareto_k_over_1_0 = n_pk_10,

  Brier_1yr = unname(brier_report["Brier_1yr"]),
  Brier_5yr = unname(brier_report["Brier_5yr"]),
  Brier_10yr = unname(brier_report["Brier_10yr"]),
  IBS_0_10yr = IBS_0_10,

  Max_Rhat = max_rhat,
  Min_ESS_Bulk = min_ess_bulk,
  Min_ESS_Tail = min_ess_tail,
  Divergences = n_divergent,
  Max_Treedepth_Hits = n_treedepth,

  # Bayes-factor component:
  # pairwise log BF(A,B) = LogML_A - LogML_B
  Log_Marginal_Likelihood = log_marginal_likelihood,
  Bridge_LogML_Mean = bridge_logml_mean,
  Bridge_LogML_SD = bridge_logml_sd,
  Bridge_LogML_IQR = bridge_logml_iqr,
  Bridge_LogML_Min = bridge_logml_min,
  Bridge_LogML_Max = bridge_logml_max,
  Bridge_MCSE_Median = bridge_mcse_median
)

cat("\n============================================\n")
cat("FINAL WEIBULL COMPARISON METRICS\n")
cat("============================================\n")
print(metrics, digits = 5)

write.csv(
  metrics,
  file.path(
    output_dir,
    "Weibull_Model_Comparison_Metrics.csv"
  ),
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  file.path(
    output_dir,
    "Weibull_Model_Comparison_Metrics.csv"
  ),
  "\n"
)
cat(
  file.path(
    output_dir,
    "LOO_Weibull.rds"
  ),
  "\n"
)
cat(
  file.path(
    output_dir,
    "BridgeSampling_Weibull.rds"
  ),
  "\n"
)
cat(
  file.path(
    output_dir,
    "Weibull_Bridge_LogML_Repetitions.csv"
  ),
  "\n"
)
cat(
  "\nFor Bayes-factor comparison, report:\n",
  "log p(y | M_Weibull) = ",
  log_marginal_likelihood,
  "\n"
)

# ============================================================
# 11. WHEN ALL THREE TEAMMATES ARE DONE
#
# Each teammate should save the LOO object:
#
#   LOO_Weibull.rds
#   LOO_Lognormal.rds
#   LOO_Loglogistic.rds
#
# Then use:
#
# loo_w  <- readRDS("LOO_Weibull.rds")
# loo_ln <- readRDS("LOO_Lognormal.rds")
# loo_ll <- readRDS("LOO_Loglogistic.rds")
#
# comparison <- loo::loo_compare(
#   Weibull = loo_w,
#   Lognormal = loo_ln,
#   Loglogistic = loo_ll
# )
#
# print(comparison)
#
# The best model appears first.
# loo_compare gives ELPD differences and their standard errors.
#
# ALL THREE MODELS MUST USE:
#   - the same patients
#   - the same row order
#   - the same covariates
#   - the same coding
#   - the same Brier/IBS time grid
# ============================================================
