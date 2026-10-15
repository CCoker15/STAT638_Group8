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
  Max_Treedepth_Hits = n_treedepth
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
