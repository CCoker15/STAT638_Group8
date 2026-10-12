# ============================================================
# PRIORS FOR STAGE, AGE, AND SEX COEFFICIENTS
#
# Weibull AFT model:
#
# log(lambda_i) =
#   beta0
# + beta2 * StageII
# + beta3 * StageIII
# + beta4 * StageIV
# + betaA * Age10
# + betaS * Male
#
# For a Weibull model with common shape k:
#
# HR = exp(-k * beta)
#
# Therefore:
#
# beta = -log(HR) / k
#
# ============================================================

set.seed(12345)

N <- 100000


# ============================================================
# 1. PRIOR FOR WEIBULL SHAPE k
#
# Current selected prior:
#
# log(k) ~ Normal(0, 0.2^2)
#
# If your final k prior changes, change these two values.
# ============================================================

mu_logk <- 0.182
sd_logk <- 0.15


logk_draws <- rnorm(
  N,
  mean = mu_logk,
  sd = sd_logk
)

k_draws <- exp(logk_draws)


# ============================================================
# 2. LITERATURE-INFORMED HR PRIOR ASSUMPTIONS
#
# IMPORTANT:
#
# HR_Lower and HR_Upper below are broad PRIOR PLAUSIBLE
# ranges. They are deliberately wider than published
# statistical confidence intervals.
#
# Age HR is specified PER YEAR.
# It will be converted to a 10-year effect because:
#
# Age10 = (Age - 60)/10
# ============================================================

HR_data <- data.frame(

  Parameter = c(
    "Stage II vs Stage I",
    "Stage III vs Stage I",
    "Stage IV vs Stage I",
    "Age: per year",
    "Male vs Female"
  ),

  HR_center = c(
    1.50,
    2.50,
    10.00,
    1.04,
    1.00
  ),

  HR_Lower = c(
    0.80,
    1.30,
    4.00,
    1.02,
    0.67
  ),

  HR_Upper = c(
    2.80,
    4.80,
    25.00,
    1.06,
    1.50
  ),

  # Effect multiplier:
  #
  # Stage = 1
  # Sex   = 1
  #
  # Age = 10 because beta_A represents
  # a 10-year increase in age.

  Multiplier = c(
    1,
    1,
    1,
    10,
    1
  )
)


print(HR_data)


# ============================================================
# 3. FUNCTION TO CONVERT HR PRIOR INTO beta PRIOR
# ============================================================

get_beta_prior <- function(
    HR_center,
    HR_Lower,
    HR_Upper,
    multiplier,
    k_draws) {

  # ----------------------------------------------------------
  # Put HR on log scale
  # ----------------------------------------------------------

  mu_logHR <- log(HR_center)

  # Treat the supplied broad HR range as an
  # approximate 95% prior interval

  sd_logHR <-
    (
      log(HR_Upper) -
      log(HR_Lower)
    ) / (2 * 1.96)


  # ----------------------------------------------------------
  # Draw possible log(HR) values
  # ----------------------------------------------------------

  logHR_draws <- rnorm(
    length(k_draws),
    mean = mu_logHR,
    sd = sd_logHR
  )


  # ----------------------------------------------------------
  # Convert age effect from 1 year to 10 years
  #
  # For age:
  #
  # log(HR_10) = 10 * log(HR_1)
  #
  # Other coefficients have multiplier = 1.
  # ----------------------------------------------------------

  logHR_model <- multiplier * logHR_draws


  # ----------------------------------------------------------
  # Convert to Weibull AFT coefficient
  #
  # beta = -log(HR) / k
  # ----------------------------------------------------------

  beta_draws <-
    -logHR_model / k_draws


  # ----------------------------------------------------------
  # Summary
  # ----------------------------------------------------------

  data.frame(

    Beta_mean =
      mean(beta_draws),

    Beta_SD =
      sd(beta_draws),

    Beta_median =
      median(beta_draws),

    Beta_lower95 =
      quantile(
        beta_draws,
        0.025,
        names = FALSE
      ),

    Beta_upper95 =
      quantile(
        beta_draws,
        0.975,
        names = FALSE
      ),

    ScaleRatio_median =
      median(
        exp(beta_draws)
      ),

    ScaleRatio_lower95 =
      quantile(
        exp(beta_draws),
        0.025,
        names = FALSE
      ),

    ScaleRatio_upper95 =
      quantile(
        exp(beta_draws),
        0.975,
        names = FALSE
      )
  )
}


# ============================================================
# 4. CALCULATE PRIORS
# ============================================================

results <- vector(
  "list",
  nrow(HR_data)
)


for (i in seq_len(nrow(HR_data))) {

  results[[i]] <-
    get_beta_prior(

      HR_center =
        HR_data$HR_center[i],

      HR_Lower =
        HR_data$HR_Lower[i],

      HR_Upper =
        HR_data$HR_Upper[i],

      multiplier =
        HR_data$Multiplier[i],

      k_draws =
        k_draws
    )
}


# ============================================================
# 5. CREATE FINAL TABLE
# ============================================================

prior_table <- do.call(
  rbind,
  results
)


prior_table <- cbind(

  Parameter =
    HR_data$Parameter,

  HR_center =
    HR_data$HR_center,

  HR_Lower =
    HR_data$HR_Lower,

  HR_Upper =
    HR_data$HR_Upper,

  prior_table
)


row.names(prior_table) <- NULL


cat(
  "\n============================================\n"
)

cat(
  "REGRESSION COEFFICIENT PRIORS\n"
)

cat(
  "============================================\n\n"
)


print(
  prior_table,
  digits = 4
)


# ============================================================
# 6. PRINT PRIORS IN MODEL NOTATION
# ============================================================

beta_names <- c(
  "beta_2",
  "beta_3",
  "beta_4",
  "beta_A",
  "beta_S"
)


cat(
  "\n============================================\n"
)

cat(
  "NORMAL PRIOR SPECIFICATIONS\n"
)

cat(
  "============================================\n\n"
)


for (i in seq_len(nrow(prior_table))) {

  cat(

    beta_names[i],
    " ~ Normal(",
    round(
      prior_table$Beta_mean[i],
      3
    ),
    ", ",
    round(
      prior_table$Beta_SD[i],
      3
    ),
    "^2)\n",

    sep = ""
  )
}


# ============================================================
# 7. SAVE RESULTS
# ============================================================

project_dir <-
  "D:/OneDrive - Texas A&M University/PhD/Courses/STAT 638/Project"


write.csv(

  prior_table,

  file.path(
    project_dir,
    "Final_Regression_Coefficient_Priors.csv"
  ),

  row.names = FALSE
)


cat(
  "\nResults saved successfully.\n"
)