###############################################################################
# STAT 638 – Group 6 – Prior derivation for the Piecewise Exponential model
#
# Purpose: translate the group's shared priors (built for the Weibull model)
#          into the parameters of the piecewise exponential model.
#   Part 1: covariate priors  -> same literature HRs, written on the log-HR scale
#   Part 2: baseline priors   -> Weibull beta0 / k prior translated into a
#                                Gamma prior for each interval hazard lambda_k
#   Part 3: prior predictive check for the reference patient
#           (Stage I, age 60, female, grade G1-G2)
#
# Run this ONCE. Copy the printed values into STEP 4 of
# piecewise_exponential_KIRC.R. (Plain R – no Stan.)
###############################################################################

library(ggplot2)

set.seed(638)
N <- 200000

# ---------------------------------------------------------------------------
# PART 1. COVARIATE PRIORS (log hazard ratio scale)
#
# Inputs = the HR centers and plausible ranges in the group write-up
# (same as Parameter_Prior_2_3_4_A_S.r, plus the two grade rows).
#
# The Weibull needs beta = -log(HR)/k because its coefficients are on the
# survival-time scale. In the piecewise model h(t) = lambda_k * exp(x'beta),
# so beta = log(HR) directly: no k, no minus sign.
#
#   mean = log(center)
#   SD   = (log(upper) - log(lower)) / (2 * 1.96)   (range = approx. 95% interval)
#   Age is given per year; x10 because the model uses Age10 = (age - 60)/10
# ---------------------------------------------------------------------------
HR <- data.frame(
  param  = c("stageII", "stageIII", "stageIV", "age10", "male", "gradeG3", "gradeG4"),
  label  = c("Stage II vs I", "Stage III vs I", "Stage IV vs I", "Age (per year)",
             "Male vs female", "Grade G3 vs G1-G2", "Grade G4 vs G1-G2"),
  center = c(1.50, 2.50, 10.00, 1.04, 1.00, 1.40, 2.50),
  lower  = c(0.80, 1.30,  4.00, 1.02, 0.67, 0.90, 1.40),
  upper  = c(2.80, 4.80, 25.00, 1.06, 1.50, 2.20, 4.50),
  mult   = c(1,    1,     1,    10,   1,    1,    1)
)

HR$mu0 <- HR$mult * log(HR$center)
HR$sd0 <- HR$mult * (log(HR$upper) - log(HR$lower)) / (2 * 1.96)

# Implied 95% prior range for the HR actually used in the model (Age = per 10 yrs)
HR$HR_lo95 <- exp(HR$mu0 - 1.96 * HR$sd0)
HR$HR_hi95 <- exp(HR$mu0 + 1.96 * HR$sd0)

cat("\n=========== PART 1: COVARIATE PRIORS (beta = log HR) ===========\n")
print(HR[, c("label", "center", "lower", "upper", "mu0", "sd0", "HR_lo95", "HR_hi95")],
      digits = 3, row.names = FALSE)
# Expect: stageII 0.405/0.320, stageIII 0.916/0.333, stageIV 2.303/0.467,
#         age10 0.392/0.098, male 0/0.206, gradeG3 0.336/0.228, gradeG4 0.916/0.298

# Cross-check against the Weibull priors: Weibull beta ~ -(log HR)/k, median k ~ 1.2
cat("\nCheck: -mu0/1.2 should be close to the Weibull prior means",
    "(-0.343, -0.773, -1.939, -0.330, 0, -0.283, -0.773):\n")
print(round(-HR$mu0 / exp(0.182), 3))

# ---------------------------------------------------------------------------
# PART 2. BASELINE PRIORS (one Gamma per interval hazard)
#
# The Weibull baseline prior (from Parameter_Prior_1.r):
#   beta0 = log(lambda0) ~ N(3.219, 0.2^2),  log(k) ~ N(0.182, 0.15^2)
# Each draw (beta0, k) is one possible survival curve for the reference patient.
#
# For each draw, the AVERAGE Weibull hazard over interval [a, b) is
#   [H(b) - H(a)] / (b - a),   with H(t) = (t / lambda0)^k
# which is exactly what a constant lambda_k represents.
# A Gamma is then matched to each interval's draws by mean and variance:
#   shape = mean^2 / var,   rate = mean / var
# Gamma keeps the lambda_k step of the Gibbs sampler conjugate (exact draws).
#
# The open-ended 5+ interval is matched over 5-8 years (where follow-up exists).
# ---------------------------------------------------------------------------
beta0  <- rnorm(N, mean = 3.219, sd = 0.20)
k      <- exp(rnorm(N, mean = 0.182, sd = 0.15))
lambda0 <- exp(beta0)
H <- function(t) (t / lambda0)^k               # Weibull cumulative hazard

lo <- c(0, 1, 2, 3, 5)
hi <- c(1, 2, 3, 5, 8)

base_tab <- do.call(rbind, lapply(1:5, function(j) {
  h <- (H(hi[j]) - H(lo[j])) / (hi[j] - lo[j])  # avg hazard per year in interval j
  m <- mean(h); v <- var(h)
  data.frame(interval    = c("0-1", "1-2", "2-3", "3-5", "5+")[j],
             hazard_mean = m,
             hazard_sd   = sqrt(v),
             hazard_lo95 = quantile(h, 0.025, names = FALSE),
             hazard_hi95 = quantile(h, 0.975, names = FALSE),
             shape       = m^2 / v,
             rate        = m / v)
}))

cat("\n=========== PART 2: BASELINE PRIORS lambda_k ~ Gamma(shape, rate) ===========\n")
print(base_tab, digits = 3, row.names = FALSE)
# Expect approx: Gamma(2.74, 113), (5.99, 209), (8.52, 273), (11.6, 343), (14.9, 402)

# ---------------------------------------------------------------------------
# PART 3. PRIOR PREDICTIVE CHECK (reference patient)
# Draw lambda_1..5 from the Gamma priors and compute S(t). The medians should
# fall inside the group's targets: S(1) 0.90-0.99, S(5) 0.80-0.95, S(10) 0.55-0.90
# ---------------------------------------------------------------------------
lam_draws <- sapply(1:5, function(j) rgamma(N, base_tab$shape[j], base_tab$rate[j]))
cut_lo <- c(0, 1, 2, 3, 5); cut_hi <- c(1, 2, 3, 5, Inf)
time_in <- function(t) pmax(0, pmin(t, cut_hi) - cut_lo)    # years in each interval by t

targets <- data.frame(Time = c(1, 5, 10), Lower = c(0.90, 0.80, 0.55), Upper = c(0.99, 0.95, 0.90))

check <- do.call(rbind, lapply(seq_len(nrow(targets)), function(i) {
  t <- targets$Time[i]
  S <- exp(-drop(lam_draws %*% time_in(t)))
  data.frame(Time = t,
             S_median = median(S),
             S_lo95   = quantile(S, 0.025, names = FALSE),
             S_hi95   = quantile(S, 0.975, names = FALSE),
             Target   = paste0(targets$Lower[i], "-", targets$Upper[i]),
             Coverage = mean(S >= targets$Lower[i] & S <= targets$Upper[i]))
}))

cat("\n=========== PART 3: PRIOR PREDICTIVE CHECK (Stage I, 60, female, G1-G2) ===========\n")
print(check, digits = 3, row.names = FALSE)
# Expect medians about 0.98, 0.86, 0.72 (all inside the targets)

# Prior predictive curve (20,000 draws to keep memory small)
tgrid <- seq(0, 10, by = 0.05)
Emat  <- t(sapply(tgrid, time_in))                      # time x 5
S_mat <- exp(-lam_draws[1:20000, ] %*% t(Emat))         # draws x time
curve <- data.frame(Time   = tgrid,
                    Median = apply(S_mat, 2, median),
                    Lo     = apply(S_mat, 2, quantile, 0.025),
                    Hi     = apply(S_mat, 2, quantile, 0.975))

p_prior <- ggplot(curve, aes(Time, Median)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = 0.2) +
  geom_line(linewidth = 1.1) +
  geom_segment(data = targets, aes(x = Time, xend = Time, y = Lower, yend = Upper),
               inherit.aes = FALSE, linewidth = 6, alpha = 0.25) +
  scale_x_continuous(breaks = 0:10) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Prior predictive survival – piecewise exponential",
       subtitle = "Reference patient: Stage I, age 60, female, grade G1–G2 (grey bars = target ranges)",
       x = "Years after diagnosis", y = "Survival probability") +
  theme_classic(base_size = 12)
print(p_prior)
#ggsave("PWE_Prior_Predictive.png", p_prior, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------------------
# PASTE-READY VALUES FOR STEP 4 OF piecewise_exponential_KIRC.R
# ---------------------------------------------------------------------------
cat("\n=========== COPY INTO STEP 4 ===========\n")
cat("a0  <- c(", paste(sprintf("%.2f", base_tab$shape), collapse = ", "), ")\n", sep = "")
cat("b0  <- c(", paste(sprintf("%.1f", base_tab$rate),  collapse = ", "), ")\n", sep = "")
cat("mu0 <- c(", paste0(HR$param, " = ", sprintf("%.3f", HR$mu0), collapse = ", "), ")\n", sep = "")
cat("sd0 <- c(", paste(sprintf("%.3f", HR$sd0), collapse = ", "), ")\n", sep = "")

#write.csv(HR,       "PWE_Covariate_Priors.csv", row.names = FALSE)
#write.csv(base_tab, "PWE_Baseline_Priors.csv",  row.names = FALSE)
