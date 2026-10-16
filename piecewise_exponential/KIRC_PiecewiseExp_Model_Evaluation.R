###############################################################################
# STAT 638 – Group 6 – Piecewise Exponential Model
# SECTION 10: Standardized survival curves + model-comparison metrics
#
# Run AFTER piecewise_exponential_KIRC_GibbsSampler.R, which must
#   (1) use time = days / 365   (same as the Weibull model), and
#   (2) save its chains with saveRDS(...) to PWE_Fit_<prior_type>.rds
#
# Produces (same definitions as KIRC_Weibull_Model_Evaluation.R):
#   10.1  data check (same 521 patients, same order, same time unit)
#   10.2  standardized survival by STAGE  (averaged over age, sex, grade)
#   10.3  standardized survival by AGE    (45/55/65/75; averaged over stage, sex, grade)
#   10.4  standardized survival by GRADE  (averaged over stage, age, sex)
#   10.5  PSIS-LOO  (ELPD, LOOIC, Pareto-k) + pointwise elpd for all 521 patients
#   10.6  IPCW Brier score at 1, 5, 10 yrs and Integrated Brier Score (0-10 yrs)
#   10.7  one-row comparison-metrics table (same columns as the Weibull CSV)
#   10.8  posterior predictive check: model vs Kaplan-Meier (overall and by stage)
###############################################################################

# install.packages(c("readxl","dplyr","survival","loo","posterior","ggplot2"))
library(readxl)
library(dplyr)
library(survival)
library(loo)
library(posterior)
library(ggplot2)

set.seed(638)

prior_type <- "informative"                          # or "vague" for the sensitivity run
fit_file   <- paste0("PWE_Fit_", prior_type, ".rds")
tag        <- paste0("PiecewiseExp_", prior_type)    # prefix for every output file

n_draws_pred <- 8000      # posterior draws used for survival predictions (Weibull used 8000)
n_draws_loo  <- 16800     # posterior draws used for PSIS-LOO (every 5th retained draw)

# =============================================================================
# 10.1  DATA – recreate the EXACT comparison dataset
#       Same filter, same order, same coding and same time unit as the Weibull script.
# =============================================================================
raw <- read_excel("KIRC - Filtered Data.xlsx")

dat <- raw %>%
  filter(toupper(trimws(histologic_grade)) %in% c("G1", "G2", "G3", "G4")) %>%   # drops GX + missing
  arrange(bcr_patient_barcode) %>%          # <- ADD: same order as the original file
  mutate(
    OS_years = as.numeric(overall_survival_time) / 365,      # SAME as Weibull (days / 365)
    event    = as.integer(overall_survival),
    Age10    = (as.numeric(age_at_diagnosis) - 60) / 10,
    Male     = as.integer(toupper(trimws(gender)) == "MALE"),
    stage    = factor(sub("Stage ", "", trimws(pathologic_stage)), levels = c("I", "II", "III", "IV")),
    grade    = factor(ifelse(histologic_grade %in% c("G1", "G2"), "G1&2", histologic_grade),
                      levels = c("G1&2", "G3", "G4"))
  )

N <- nrow(dat)
stopifnot(N == 521, sum(dat$event) == 173, !anyNA(dat$stage))
cat("Comparison N =", N, " deaths =", sum(dat$event), "\n")

# Patient-level design matrix: same column names as the posterior draws
X_pat <- cbind(
  stageII  = as.numeric(dat$stage == "II"),
  stageIII = as.numeric(dat$stage == "III"),
  stageIV  = as.numeric(dat$stage == "IV"),
  age10    = dat$Age10,
  male     = dat$Male,
  gradeG3  = as.numeric(dat$grade == "G3"),
  gradeG4  = as.numeric(dat$grade == "G4")
)
beta_names <- colnames(X_pat)
lam_names  <- paste0("lambda", 1:5)

# =============================================================================
# Load the saved posterior chains
# =============================================================================
fit <- readRDS(fit_file)
stopifnot(fit$time_divisor == 365)               # must match the Weibull time unit
cut_all <- fit$cut_all                           # 0, 1, 2, 3, 5, Inf
lo <- head(cut_all, -1); hi <- tail(cut_all, -1)
time_in <- function(t) pmax(0, pmin(t, hi) - lo) # years spent in each interval by time t

# Retained draws per chain (burn-in removed)
chains  <- lapply(fit$fits, function(f) f$draws[-(1:fit$burn), c(lam_names, beta_names)])
n_chain <- length(chains)
n_iter  <- nrow(chains[[1]])

# Evenly spaced draws from every chain, keeping track of which chain each came from
take_draws <- function(n_total) {
  idx <- round(seq(1, n_iter, length.out = ceiling(n_total / n_chain)))
  list(draws = do.call(rbind, lapply(chains, function(d) d[idx, , drop = FALSE])),
       chain = rep(seq_len(n_chain), each = length(idx)))
}

# =============================================================================
# MCMC diagnostics (same quantities as the Weibull table)
# =============================================================================
arr <- array(NA_real_, dim = c(n_iter, n_chain, ncol(chains[[1]])),
             dimnames = list(NULL, NULL, colnames(chains[[1]])))
for (ch in seq_len(n_chain)) arr[, ch, ] <- chains[[ch]]
diag_tab <- summarise_draws(as_draws_array(arr), "rhat", "ess_bulk", "ess_tail")
print(diag_tab, n = Inf)
max_rhat     <- max(diag_tab$rhat)
min_ess_bulk <- min(diag_tab$ess_bulk)
min_ess_tail <- min(diag_tab$ess_tail)

# =============================================================================
# Prediction pieces shared by 10.2-10.4, 10.6 and 10.8
#   S_i(t) = exp( -exp(x_i'beta) * sum_k lambda_k * time_in_k(t) )   for each draw
# =============================================================================
D     <- take_draws(n_draws_pred)
lamD  <- D$draws[, lam_names]                    # M x 5
betaD <- D$draws[, beta_names]                   # M x 7
M     <- nrow(lamD)

# Survival at time t for every draw (rows) and patient (cols), given a design matrix
surv_matrix <- function(t, X) {
  if (t == 0) return(matrix(1, M, nrow(X)))
  H0 <- drop(lamD %*% time_in(t))                # baseline cumulative hazard, one per draw
  exp(-exp(betaD %*% t(X)) * H0)                 # (M x N) * length-M vector -> row-wise
}

time_grid      <- round(seq(0, 10, by = 0.1), 2)
times_interest <- c(1, 5, 10)

# =============================================================================
# 10.2-10.4  STANDARDIZED (MARGINALIZED) SURVIVAL CURVES
#   For each posterior draw m:
#     Sbar^(m)(t) = (1/N) * sum_i S^(m)(t | fixed value, patient i's other covariates)
#   then median and 95% CrI across draws.
# =============================================================================
set_fixed <- function(X, fix, value) {
  if (fix == "stage") {
    X[, c("stageII", "stageIII", "stageIV")] <- 0
    if (value != "I") X[, paste0("stage", value)] <- 1
  }
  if (fix == "age") X[, "age10"] <- (as.numeric(value) - 60) / 10
  if (fix == "grade") {
    X[, c("gradeG3", "gradeG4")] <- 0
    if (value != "G1&2") X[, paste0("grade", value)] <- 1
  }
  X
}

standardized_curve <- function(fix, value) {
  Xf  <- set_fixed(X_pat, fix, value)
  elp <- exp(betaD %*% t(Xf))                    # M x N, computed once per level
  bind_rows(lapply(time_grid, function(t) {
    Sbar <- if (t == 0) rep(1, M) else rowMeans(exp(-elp * drop(lamD %*% time_in(t))))
    data.frame(Time = t,
               Posterior_mean   = mean(Sbar),
               Posterior_median = median(Sbar),
               Lower_95_CrI     = quantile(Sbar, 0.025, names = FALSE),
               Upper_95_CrI     = quantile(Sbar, 0.975, names = FALSE))
  })) %>% mutate(Level = as.character(value))
}

run_study <- function(fix, values, labels) {
  cat("Standardized curves by", fix, "...\n")
  out <- bind_rows(lapply(values, function(v) standardized_curve(fix, v)))
  out$Level <- factor(out$Level, levels = as.character(values), labels = labels)
  out
}

stage_curves <- run_study("stage", c("I", "II", "III", "IV"),
                          c("Stage I", "Stage II", "Stage III", "Stage IV"))
age_curves   <- run_study("age",   c(45, 55, 65, 75),
                          c("45 years", "55 years", "65 years", "75 years"))
grade_curves <- run_study("grade", c("G1&2", "G3", "G4"),
                          c("G1-G2", "G3", "G4"))
# (a few minutes: 11 levels x 101 time points x 8,000 draws x 521 patients)

summary_1_5_10 <- function(curves) {
  curves %>%
    filter(Time %in% times_interest) %>%
    transmute(Level, Time_years = Time, Posterior_mean, Posterior_median,
              Lower_95_CrI, Upper_95_CrI,
              Median_percent   = 100 * Posterior_median,
              Lower_95_percent = 100 * Lower_95_CrI,
              Upper_95_percent = 100 * Upper_95_CrI)
}
stage_tab <- summary_1_5_10(stage_curves); print(stage_tab, digits = 3)
age_tab   <- summary_1_5_10(age_curves);   print(age_tab,   digits = 3)
grade_tab <- summary_1_5_10(grade_curves); print(grade_tab, digits = 3)

plot_std <- function(curves, title, subtitle, legend_title) {
  ggplot(curves, aes(Time, Posterior_median, colour = Level, fill = Level)) +
    geom_ribbon(aes(ymin = Lower_95_CrI, ymax = Upper_95_CrI), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 1) +
    scale_x_continuous(limits = c(0, 10), breaks = 0:10) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       labels = function(x) paste0(round(100 * x), "%")) +
    labs(title = title, subtitle = subtitle, x = "Years after diagnosis",
         y = "Survival probability", colour = legend_title, fill = legend_title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
}

p_stage <- plot_std(stage_curves, "Standardized Posterior Survival by Pathologic Stage",
                    "Piecewise exponential; averaged over observed age, sex and grade", "Pathologic stage")
p_age   <- plot_std(age_curves, "Standardized Posterior Survival by Age",
                    "Piecewise exponential; averaged over observed stage, sex and grade", "Age")
p_grade <- plot_std(grade_curves, "Standardized Posterior Survival by Histologic Grade",
                    "Piecewise exponential; averaged over observed stage, age and sex", "Grade")
print(p_stage); print(p_age); print(p_grade)

#ggsave(paste0(tag, "_Standardized_Stage.png"), p_stage, width = 8, height = 6, dpi = 300)
#ggsave(paste0(tag, "_Standardized_Age.png"),   p_age,   width = 8, height = 6, dpi = 300)
#ggsave(paste0(tag, "_Standardized_Grade.png"), p_grade, width = 8, height = 6, dpi = 300)
#write.csv(stage_curves, paste0(tag, "_Standardized_Curves_Stage.csv"), row.names = FALSE)
#write.csv(age_curves,   paste0(tag, "_Standardized_Curves_Age.csv"),   row.names = FALSE)
#write.csv(grade_curves, paste0(tag, "_Standardized_Curves_Grade.csv"), row.names = FALSE)
#write.csv(stage_tab, paste0(tag, "_Standardized_Stage_1_5_10.csv"), row.names = FALSE)
#write.csv(age_tab,   paste0(tag, "_Standardized_Age_1_5_10.csv"),   row.names = FALSE)
#write.csv(grade_tab, paste0(tag, "_Standardized_Grade_1_5_10.csv"), row.names = FALSE)

# =============================================================================
# 10.5  PSIS-LOO
#   Pointwise log-likelihood per PATIENT (not per long-format row), using the
#   survival likelihood so it is on the same scale as the Weibull log_lik:
#     death:    log f(t_i) = log h(t_i) - H(t_i) = log lambda_{k(t_i)} + x_i'beta - H(t_i)
#     censored: log S(t_i) = -H(t_i)
#   (the Poisson form's extra log(exposure) term must NOT be included)
# =============================================================================
L     <- take_draws(n_draws_loo)
lamL  <- L$draws[, lam_names]
lpL   <- L$draws[, beta_names] %*% t(X_pat)                  # M_loo x N : x_i'beta
E_pat <- t(sapply(dat$OS_years, time_in))                     # N x 5 : time in each interval
k_pat <- findInterval(dat$OS_years, cut_all[2:5], left.open = TRUE) + 1   # interval of t_i

H_pat   <- (lamL %*% t(E_pat)) * exp(lpL)                     # cumulative hazard at t_i
log_haz <- log(lamL[, k_pat]) + lpL                           # log hazard at t_i
log_lik <- sweep(log_haz, 2, dat$event, "*") - H_pat          # M_loo x N
stopifnot(ncol(log_lik) == N, all(is.finite(log_lik)))
rm(H_pat, log_haz, lpL); invisible(gc())

r_eff   <- relative_eff(exp(log_lik), chain_id = L$chain)
loo_pwe <- loo(log_lik, r_eff = r_eff, cores = 1)
print(loo_pwe)

saveRDS(loo_pwe, paste0("LOO_", tag, ".rds"))                 # contains pointwise elpd_loo

elpd_loo <- loo_pwe$estimates["elpd_loo", "Estimate"]
se_elpd  <- loo_pwe$estimates["elpd_loo", "SE"]
p_loo    <- loo_pwe$estimates["p_loo", "Estimate"]
looic    <- loo_pwe$estimates["looic", "Estimate"]
se_looic <- loo_pwe$estimates["looic", "SE"]
pk       <- pareto_k_values(loo_pwe)

id_col <- intersect(c("sample", "bcr_patient_barcode", "patient_id"), names(dat))[1]

pointwise <- data.frame(row      = seq_len(N),
                        id       = if (!is.na(id_col)) dat[[id_col]] else NA,
                        event    = dat$event,
                        OS_years = dat$OS_years,
                        elpd_loo = loo_pwe$pointwise[, "elpd_loo"],
                        pareto_k = pk)
write.csv(pointwise, paste0(tag, "_Pointwise_ELPD.csv"), row.names = FALSE)

cat("\nELPD-LOO =", round(elpd_loo, 2), " SE =", round(se_elpd, 2),
    "\nLOOIC =", round(looic, 2), "  max Pareto-k =", round(max(pk), 3), "\n")

# =============================================================================
# 10.6  IPCW BRIER SCORE AND INTEGRATED BRIER SCORE  (same code as the Weibull script)
# =============================================================================
post_mean_survival <- function(t) colMeans(surv_matrix(t, X_pat))   # length N

cens_fit <- survfit(Surv(OS_years, 1 - event) ~ 1, data = dat)      # KM of censoring
G_right <- function(x) vapply(x, function(z) {
  idx <- which(cens_fit$time <= z); if (length(idx) == 0) 1 else cens_fit$surv[max(idx)] }, numeric(1))
G_left  <- function(x) vapply(x, function(z) {
  idx <- which(cens_fit$time <  z); if (length(idx) == 0) 1 else cens_fit$surv[max(idx)] }, numeric(1))
G_Tminus <- G_left(dat$OS_years)

ipcw_brier <- function(t, S_hat) {
  G_t <- G_right(t)
  contribution <- rep(0, N)
  died <- dat$OS_years <= t & dat$event == 1                       # died by t
  contribution[died] <- S_hat[died]^2 / pmax(G_Tminus[died], 1e-8)
  alive <- dat$OS_years > t                                         # known alive at t
  contribution[alive] <- (1 - S_hat[alive])^2 / pmax(G_t, 1e-8)
  mean(contribution)                                                # censored before t -> 0
}

brier_report <- sapply(times_interest, function(t) ipcw_brier(t, post_mean_survival(t)))
names(brier_report) <- paste0("Brier_", times_interest, "yr")
print(round(brier_report, 4))

ibs_grid <- seq(0, 10, by = 0.25)                                   # SAME grid as Weibull
bs_grid  <- sapply(ibs_grid, function(t) ipcw_brier(t, post_mean_survival(t)))
IBS_0_10 <- sum(diff(ibs_grid) * (head(bs_grid, -1) + tail(bs_grid, -1)) / 2) / 10
cat("IBS 0-10 years =", round(IBS_0_10, 4), "\n")

write.csv(data.frame(Time_years = ibs_grid, Brier = bs_grid),
          paste0(tag, "_Brier_Curve.csv"), row.names = FALSE)

# =============================================================================
# 10.7  ONE-ROW METRICS TABLE (same columns as Weibull_Model_Comparison_Metrics.csv)
#   Divergences / tree-depth hits exist only for Stan's sampler -> NA here.
# =============================================================================
metrics <- data.frame(
  Model = "Piecewise exponential", N = N,
  ELPD_LOO = elpd_loo, SE_ELPD_LOO = se_elpd, p_LOO = p_loo,
  LOOIC = looic, SE_LOOIC = se_looic,
  Max_Pareto_k = max(pk), N_Pareto_k_over_0_7 = sum(pk > 0.7), N_Pareto_k_over_1_0 = sum(pk > 1.0),
  Brier_1yr = unname(brier_report["Brier_1yr"]),
  Brier_5yr = unname(brier_report["Brier_5yr"]),
  Brier_10yr = unname(brier_report["Brier_10yr"]),
  IBS_0_10yr = IBS_0_10,
  Max_Rhat = max_rhat, Min_ESS_Bulk = min_ess_bulk, Min_ESS_Tail = min_ess_tail,
  Divergences = NA, Max_Treedepth_Hits = NA
)
print(metrics, digits = 5)
write.csv(metrics, paste0(tag, "_Model_Comparison_Metrics.csv"), row.names = FALSE)

# Side by side with the Weibull, if its files are in this folder
if (file.exists("Weibull_Model_Comparison_Metrics.csv")) {
  print(bind_rows(read.csv("Weibull_Model_Comparison_Metrics.csv"), metrics), digits = 5)
}
if (file.exists("LOO_Weibull.rds")) {
  print(loo_compare(list(Weibull = readRDS("LOO_Weibull.rds"),
                         PiecewiseExp = loo_pwe)))
}

# =============================================================================
# 10.8  POSTERIOR PREDICTIVE CHECK: model vs Kaplan-Meier
#   Model curve for a group = average of the posterior survival of the patients
#   actually in that group (their own covariates), compared with the group's KM.
# =============================================================================
pp_group <- function(rows, label) {
  bind_rows(lapply(time_grid, function(t) {
    Sbar <- rowMeans(surv_matrix(t, X_pat[rows, , drop = FALSE]))
    data.frame(Time = t, Median = median(Sbar),
               Lo = quantile(Sbar, 0.025, names = FALSE), Hi = quantile(Sbar, 0.975, names = FALSE),
               Group = label)
  }))
}
groups   <- c(list(All = seq_len(N)), split(seq_len(N), paste("Stage", dat$stage)))
pp_model <- bind_rows(lapply(names(groups), function(g) pp_group(groups[[g]], g)))

km_df <- bind_rows(lapply(names(groups), function(g) {
  km <- survfit(Surv(OS_years, event) ~ 1, data = dat[groups[[g]], ])
  data.frame(Time = c(0, km$time), KM = c(1, km$surv), Group = g)
}))

p_pp <- ggplot() +
  geom_ribbon(data = pp_model, aes(Time, ymin = Lo, ymax = Hi), fill = "steelblue", alpha = 0.2) +
  geom_line(data = pp_model, aes(Time, Median), colour = "steelblue", linewidth = 1) +
  geom_step(data = km_df, aes(Time, KM), colour = "black") +
  facet_wrap(~ Group) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 1)) +
  labs(title = "Posterior predictive check: piecewise exponential vs Kaplan-Meier",
       subtitle = "Blue = model (median, 95% CrI); black = Kaplan-Meier",
       x = "Years after diagnosis", y = "Survival probability") +
  theme_classic(base_size = 11)
print(p_pp)
ggsave(paste0(tag, "_PPC_vs_KM.png"), p_pp, width = 10, height = 7, dpi = 300)

cat("\nDone. Files written with prefix", tag, "\n")
