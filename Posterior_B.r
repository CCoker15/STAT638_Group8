# ============================================================
# BAYESIAN EMPIRICAL STANDARDIZATION - KIRC
#
# IMPORTANT:
#   This script DOES NOT run MCMC again.
#   It reuses the exact retained posterior draws from the
#   Bayesian Weibull model.
#
# Analyses:
#   A. Stage-specific standardized survival
#      Fix stage; marginalize over observed age + sex
#
#   B. Age-specific standardized survival
#      Fix age at 45, 55, 65, 75;
#      marginalize over observed stage + sex
#
#   C. Sex-specific standardized survival
#      Fix sex; marginalize over observed age + stage
#
# Outputs:
#   - posterior median survival curves
#   - 95% credible intervals
#   - survival at 1, 5, and 10 years
#   - CSV files and PNG figures
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. PROJECT DIRECTORY
# ============================================================

project_dir <- "D:/OneDrive - Texas A&M University/PhD/Courses/STAT 638/Project"

data_file <- file.path(
  project_dir,
  "KIRC - Filtered Data.xlsx"
)

posterior_file <- file.path(
  project_dir,
  "KIRC_Weibull_Posterior_Draws.csv"
)

if (!file.exists(data_file)) {
  stop("KIRC data file not found: ", data_file)
}

if (!file.exists(posterior_file)) {
  stop(
    "Posterior draws file not found: ", posterior_file,
    "\nRun the Bayesian Weibull MCMC model first."
  )
}

# ============================================================
# 2. READ KIRC DATA
# ============================================================

dat <- read_excel(data_file)

cat("\n============================================\n")
cat("KIRC DATA\n")
cat("============================================\n")
cat("Number of patients =", nrow(dat), "\n")

required_data_columns <- c(
  "age_at_diagnosis",
  "gender",
  "pathologic_stage"
)

missing_data_columns <- setdiff(
  required_data_columns,
  names(dat)
)

if (length(missing_data_columns) > 0) {
  stop(
    "Missing data columns: ",
    paste(missing_data_columns, collapse = ", ")
  )
}

# ============================================================
# 3. PREPARE AGE
#
# age_at_diagnosis is in years.
# Age10 = (Age - 60) / 10
# ============================================================

dat$Age_years <- as.numeric(dat$age_at_diagnosis)
dat$Age10 <- (dat$Age_years - 60) / 10

# ============================================================
# 4. PREPARE SEX
#
# Female = 0
# Male   = 1
# ============================================================

sex_text <- toupper(
  trimws(
    as.character(dat$gender)
  )
)

if (!all(sex_text %in% c("FEMALE", "MALE"))) {
  stop(
    "Unexpected gender values found: ",
    paste(
      unique(sex_text[!sex_text %in% c("FEMALE", "MALE")]),
      collapse = ", "
    )
  )
}

dat$Male <- ifelse(
  sex_text == "MALE",
  1L,
  0L
)

cat("\nSex distribution:\n")
print(table(dat$gender))

# ============================================================
# 5. PREPARE STAGE INDICATORS
#
# Stage I is the reference category.
# ============================================================

stage_text <- trimws(
  as.character(dat$pathologic_stage)
)

valid_stages <- c(
  "Stage I",
  "Stage II",
  "Stage III",
  "Stage IV"
)

if (!all(stage_text %in% valid_stages)) {
  stop(
    "Unexpected stage values found: ",
    paste(
      unique(stage_text[!stage_text %in% valid_stages]),
      collapse = ", "
    )
  )
}

dat$StageII  <- as.integer(stage_text == "Stage II")
dat$StageIII <- as.integer(stage_text == "Stage III")
dat$StageIV  <- as.integer(stage_text == "Stage IV")

cat("\nStage distribution:\n")
print(table(stage_text))

# Check for missing covariates used in standardization
needed_complete <- complete.cases(
  dat$Age10,
  dat$Male,
  dat$StageII,
  dat$StageIII,
  dat$StageIV
)

if (!all(needed_complete)) {
  stop(
    sum(!needed_complete),
    " patients have missing age, sex, or stage."
  )
}

N <- nrow(dat)

# ============================================================
# 6. READ EXACT RETAINED POSTERIOR DRAWS
# ============================================================

post <- read.csv(posterior_file)

# If k was not saved but log_k was, reconstruct it.
if (!"k" %in% names(post)) {
  if (!"log_k" %in% names(post)) {
    stop("Neither k nor log_k was found in the posterior draws.")
  }

  post$k <- exp(post$log_k)
}

required_parameters <- c(
  "beta0",
  "beta2",
  "beta3",
  "beta4",
  "betaA",
  "betaS",
  "k"
)

missing_parameters <- setdiff(
  required_parameters,
  names(post)
)

if (length(missing_parameters) > 0) {
  stop(
    "Missing posterior parameters: ",
    paste(missing_parameters, collapse = ", ")
  )
}

# Keep only finite rows for the required parameters.
finite_rows <- apply(
  post[, required_parameters, drop = FALSE],
  1,
  function(x) all(is.finite(x))
)

if (!all(finite_rows)) {
  stop(
    sum(!finite_rows),
    " posterior rows contain non-finite parameter values."
  )
}

M <- nrow(post)

cat("\n============================================\n")
cat("POSTERIOR DRAWS\n")
cat("============================================\n")
cat("Posterior draws being used =", M, "\n")
cat("Patients used for standardization =", N, "\n")

if (M != 8000) {
  cat(
    "\nNOTE: Expected 8000 draws if using 4 chains x 2000 retained,\n",
    "but found ", M, ". The script will use all available draws.\n",
    sep = ""
  )
}

# ============================================================
# 7. TIME GRID
# ============================================================

time_grid <- seq(
  0,
  10,
  by = 0.1
)

times_interest <- c(
  1,
  5,
  10
)

# ============================================================
# 8. HELPER FUNCTION
#
# eta_matrix:
#   rows    = posterior draws
#   columns = empirical covariate groups
#
# weights:
#   empirical frequency of each covariate group
#
# For each posterior draw m:
#
#   Sbar^(m)(t)
#   = sum_g w_g * exp[-(t / lambda_g^(m))^k^(m)]
#
# where lambda = exp(eta).
#
# Grouping identical covariate patterns is EXACTLY equivalent
# to predicting all 529 patients individually and averaging.
# ============================================================

standardized_survival_from_eta <- function(
    eta_matrix,
    weights,
    k_draws,
    times) {

  if (nrow(eta_matrix) != length(k_draws)) {
    stop("eta_matrix rows must equal number of k draws.")
  }

  if (ncol(eta_matrix) != length(weights)) {
    stop("eta_matrix columns must equal number of weights.")
  }

  if (abs(sum(weights) - 1) > 1e-10) {
    stop("Empirical weights do not sum to 1.")
  }

  M_local <- nrow(eta_matrix)

  # For each posterior draw and empirical group:
  #
  # S(t) = exp[-(t / exp(eta))^k]
  #      = exp[-t^k * exp(-k*eta)]
  #
  # k_draws has length M, so R recycles it down each column,
  # preserving the same k for every empirical group within a
  # given posterior draw.

  k_eta <- eta_matrix * k_draws
  A <- exp(-k_eta)

  results <- vector(
    "list",
    length(times)
  )

  for (j in seq_along(times)) {

    t <- times[j]

    if (t == 0) {

      Sbar_draws <- rep(
        1,
        M_local
      )

    } else {

      t_k <- t^k_draws

      # M x G matrix:
      # one survival probability for every posterior draw
      # and every empirical covariate group.
      S_matrix <- exp(
        -A * t_k
      )

      # Weighted average over the empirical covariate
      # distribution, separately for every posterior draw.
      Sbar_draws <- drop(
        S_matrix %*% weights
      )
    }

    results[[j]] <- data.frame(
      Time = t,
      Posterior_mean = mean(Sbar_draws),
      Posterior_median = median(Sbar_draws),
      Lower_95_CrI = quantile(
        Sbar_draws,
        0.025,
        names = FALSE
      ),
      Upper_95_CrI = quantile(
        Sbar_draws,
        0.975,
        names = FALSE
      )
    )
  }

  bind_rows(results)
}

# ============================================================
# ANALYSIS A:
# STANDARDIZED POSTERIOR SURVIVAL BY STAGE
#
# FIX:
#   Stage = I, II, III, or IV
#
# MARGINALIZE OVER:
#   observed age + sex distribution
#
# Sbar_s^(m)(t)
# = 1/N sum_i S^(m)(t | Stage=s, Age_i, Sex_i)
# ============================================================

# Empirical joint distribution of age and sex
age_sex_distribution <- dat %>%
  count(
    Age10,
    Male,
    name = "n"
  ) %>%
  mutate(
    weight = n / N
  )

cat(
  "\nUnique age-sex combinations =",
  nrow(age_sex_distribution),
  "\n"
)

stage_levels <- c(
  "Stage I",
  "Stage II",
  "Stage III",
  "Stage IV"
)

stage_curve_list <- list()

for (stage in stage_levels) {

  cat(
    "\nCalculating standardized curve for",
    stage,
    "...\n"
  )

  stage_effect <- switch(
    stage,
    "Stage I"   = rep(0, M),
    "Stage II"  = post$beta2,
    "Stage III" = post$beta3,
    "Stage IV"  = post$beta4
  )

  G <- nrow(age_sex_distribution)

  # Everyone is assigned the SAME stage.
  # Age and sex follow the observed empirical distribution.
  eta <- outer(
    post$beta0 + stage_effect,
    rep(1, G)
  ) +
    outer(
      post$betaA,
      age_sex_distribution$Age10
    ) +
    outer(
      post$betaS,
      age_sex_distribution$Male
    )

  tmp <- standardized_survival_from_eta(
    eta_matrix = eta,
    weights = age_sex_distribution$weight,
    k_draws = post$k,
    times = time_grid
  )

  tmp$Stage <- stage

  # CORRECT LIST:
  stage_curve_list[[stage]] <- tmp

  rm(eta)
  gc()
}

stage_curves <- bind_rows(
  stage_curve_list
)

stage_curves$Stage <- factor(
  stage_curves$Stage,
  levels = stage_levels
)

# ============================================================
# ANALYSIS B:
# STANDARDIZED POSTERIOR SURVIVAL BY AGE
#
# FIX:
#   Age = 45, 55, 65, 75
#
# MARGINALIZE OVER:
#   observed stage + sex distribution
#
# Sbar_a^(m)(t)
# = 1/N sum_i S^(m)(t | Age=a, Stage_i, Sex_i)
# ============================================================

ages_to_compare <- c(
  45,
  55,
  65,
  75
)

# Empirical joint distribution of stage and sex
stage_sex_distribution <- dat %>%
  count(
    StageII,
    StageIII,
    StageIV,
    Male,
    name = "n"
  ) %>%
  mutate(
    weight = n / N
  )

cat(
  "\nUnique stage-sex combinations =",
  nrow(stage_sex_distribution),
  "\n"
)

age_curve_list <- list()

for (age_value in ages_to_compare) {

  cat(
    "\nCalculating standardized curve for age",
    age_value,
    "...\n"
  )

  age10_fixed <- (
    age_value - 60
  ) / 10

  G <- nrow(stage_sex_distribution)

  # Everyone is assigned the SAME age.
  # Stage and sex follow their observed joint distribution.
  eta <- outer(
    post$beta0 +
      post$betaA * age10_fixed,
    rep(1, G)
  ) +
    outer(
      post$beta2,
      stage_sex_distribution$StageII
    ) +
    outer(
      post$beta3,
      stage_sex_distribution$StageIII
    ) +
    outer(
      post$beta4,
      stage_sex_distribution$StageIV
    ) +
    outer(
      post$betaS,
      stage_sex_distribution$Male
    )

  tmp <- standardized_survival_from_eta(
    eta_matrix = eta,
    weights = stage_sex_distribution$weight,
    k_draws = post$k,
    times = time_grid
  )

  tmp$Age <- age_value

  # CORRECT DOUBLE-BRACKET LIST ASSIGNMENT:
  age_curve_list[[as.character(age_value)]] <- tmp

  rm(eta)
  gc()
}

age_curves <- bind_rows(
  age_curve_list
)

# ============================================================
# ANALYSIS C:
# STANDARDIZED POSTERIOR SURVIVAL BY SEX
#
# FIX:
#   Female or Male
#
# MARGINALIZE OVER:
#   observed age + stage distribution
#
# Sbar_x^(m)(t)
# = 1/N sum_i S^(m)(t | Sex=x, Age_i, Stage_i)
# ============================================================

# Empirical joint distribution of age and stage
age_stage_distribution <- dat %>%
  count(
    Age10,
    StageII,
    StageIII,
    StageIV,
    name = "n"
  ) %>%
  mutate(
    weight = n / N
  )

cat(
  "\nUnique age-stage combinations =",
  nrow(age_stage_distribution),
  "\n"
)

sex_levels <- c(
  "Female",
  "Male"
)

sex_curve_list <- list()

for (sex_value in sex_levels) {

  cat(
    "\nCalculating standardized curve for",
    sex_value,
    "...\n"
  )

  male_fixed <- ifelse(
    sex_value == "Male",
    1,
    0
  )

  G <- nrow(age_stage_distribution)

  # Everyone is assigned the SAME sex.
  # Age and stage follow their observed joint distribution.
  eta <- outer(
    post$beta0 +
      post$betaS * male_fixed,
    rep(1, G)
  ) +
    outer(
      post$betaA,
      age_stage_distribution$Age10
    ) +
    outer(
      post$beta2,
      age_stage_distribution$StageII
    ) +
    outer(
      post$beta3,
      age_stage_distribution$StageIII
    ) +
    outer(
      post$beta4,
      age_stage_distribution$StageIV
    )

  tmp <- standardized_survival_from_eta(
    eta_matrix = eta,
    weights = age_stage_distribution$weight,
    k_draws = post$k,
    times = time_grid
  )

  tmp$Sex <- sex_value

  sex_curve_list[[sex_value]] <- tmp

  rm(eta)
  gc()
}

sex_curves <- bind_rows(
  sex_curve_list
)

sex_curves$Sex <- factor(
  sex_curves$Sex,
  levels = sex_levels
)

# ============================================================
# 9. EXTRACT 1-, 5-, AND 10-YEAR RESULTS
# ============================================================

stage_summary <- stage_curves %>%
  filter(
    Time %in% times_interest
  ) %>%
  select(
    Stage,
    Time,
    Posterior_mean,
    Posterior_median,
    Lower_95_CrI,
    Upper_95_CrI
  )

age_summary <- age_curves %>%
  filter(
    Time %in% times_interest
  ) %>%
  select(
    Age,
    Time,
    Posterior_mean,
    Posterior_median,
    Lower_95_CrI,
    Upper_95_CrI
  )

sex_summary <- sex_curves %>%
  filter(
    Time %in% times_interest
  ) %>%
  select(
    Sex,
    Time,
    Posterior_mean,
    Posterior_median,
    Lower_95_CrI,
    Upper_95_CrI
  )

# Add percentage columns
add_percent_columns <- function(x) {
  x %>%
    mutate(
      Mean_percent = 100 * Posterior_mean,
      Median_percent = 100 * Posterior_median,
      Lower_95_percent = 100 * Lower_95_CrI,
      Upper_95_percent = 100 * Upper_95_CrI
    )
}

stage_summary <- add_percent_columns(
  stage_summary
)

age_summary <- add_percent_columns(
  age_summary
)

sex_summary <- add_percent_columns(
  sex_summary
)

# ============================================================
# 10. PRINT RESULTS
# ============================================================

cat("\n\n============================================\n")
cat("STANDARDIZED SURVIVAL BY STAGE\n")
cat("Age and sex marginalized out\n")
cat("============================================\n\n")
print(stage_summary, digits = 4)

cat("\n\n============================================\n")
cat("STANDARDIZED SURVIVAL BY AGE\n")
cat("Stage and sex marginalized out\n")
cat("============================================\n\n")
print(age_summary, digits = 4)

cat("\n\n============================================\n")
cat("STANDARDIZED SURVIVAL BY SEX\n")
cat("Age and stage marginalized out\n")
cat("============================================\n\n")
print(sex_summary, digits = 4)

# ============================================================
# 11. SAVE NUMERICAL RESULTS
# ============================================================

write.csv(
  stage_curves,
  file.path(
    project_dir,
    "Standardized_Posterior_Curves_Stage.csv"
  ),
  row.names = FALSE
)

write.csv(
  age_curves,
  file.path(
    project_dir,
    "Standardized_Posterior_Curves_Age.csv"
  ),
  row.names = FALSE
)

write.csv(
  sex_curves,
  file.path(
    project_dir,
    "Standardized_Posterior_Curves_Sex.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage_summary,
  file.path(
    project_dir,
    "Standardized_Stage_Survival_1_5_10.csv"
  ),
  row.names = FALSE
)

write.csv(
  age_summary,
  file.path(
    project_dir,
    "Standardized_Age_Survival_1_5_10.csv"
  ),
  row.names = FALSE
)

write.csv(
  sex_summary,
  file.path(
    project_dir,
    "Standardized_Sex_Survival_1_5_10.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 12. PLOT SETTINGS
# ============================================================

font_name <- "Times New Roman"

base_theme <- theme_classic(
  base_family = font_name
) +
  theme(
    text = element_text(
      family = font_name,
      size = 12
    ),
    plot.title = element_text(
      family = font_name,
      face = "bold",
      size = 14,
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      family = font_name,
      size = 11,
      hjust = 0.5
    ),
    axis.title = element_text(
      family = font_name,
      size = 12
    ),
    axis.text = element_text(
      family = font_name,
      size = 11
    ),
    legend.title = element_text(
      family = font_name,
      size = 11
    ),
    legend.text = element_text(
      family = font_name,
      size = 10
    )
  )

# ============================================================
# 13. STAGE PLOT
# ============================================================

p_stage <- ggplot(
  stage_curves,
  aes(
    x = Time,
    y = Posterior_median,
    colour = Stage,
    fill = Stage
  )
) +
  geom_ribbon(
    aes(
      ymin = Lower_95_CrI,
      ymax = Upper_95_CrI
    ),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = 0:10
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) {
      paste0(
        round(100 * x),
        "%"
      )
    }
  ) +
  labs(
    title = "Standardized Posterior Survival by Pathologic Stage",
    subtitle = "Marginalized over the observed age and sex distribution",
    x = "Years after diagnosis",
    y = "Survival probability",
    colour = "Pathologic Stage",
    fill = "Pathologic Stage"
  ) +
  base_theme

print(p_stage)

# ============================================================
# 14. AGE PLOT
# ============================================================

age_curves$Age_label <- factor(
  paste0(
    age_curves$Age,
    " years"
  ),
  levels = paste0(
    ages_to_compare,
    " years"
  )
)

p_age <- ggplot(
  age_curves,
  aes(
    x = Time,
    y = Posterior_median,
    colour = Age_label,
    fill = Age_label
  )
) +
  geom_ribbon(
    aes(
      ymin = Lower_95_CrI,
      ymax = Upper_95_CrI
    ),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = 0:10
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) {
      paste0(
        round(100 * x),
        "%"
      )
    }
  ) +
  labs(
    title = "Standardized Posterior Survival by Age",
    subtitle = "Marginalized over the observed stage and sex distribution",
    x = "Years after diagnosis",
    y = "Survival probability",
    colour = "Age",
    fill = "Age"
  ) +
  base_theme

print(p_age)

# ============================================================
# 15. SEX PLOT
# ============================================================

p_sex <- ggplot(
  sex_curves,
  aes(
    x = Time,
    y = Posterior_median,
    colour = Sex,
    fill = Sex
  )
) +
  geom_ribbon(
    aes(
      ymin = Lower_95_CrI,
      ymax = Upper_95_CrI
    ),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(
    linewidth = 1
  ) +
  scale_x_continuous(
    limits = c(0, 10),
    breaks = 0:10
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) {
      paste0(
        round(100 * x),
        "%"
      )
    }
  ) +
  labs(
    title = "Standardized Posterior Survival by Sex",
    subtitle = "Marginalized over the observed age and stage distribution",
    x = "Years after diagnosis",
    y = "Survival probability",
    colour = "Sex",
    fill = "Sex"
  ) +
  base_theme

print(p_sex)

# ============================================================
# 16. SAVE FIGURES
# ============================================================

ggsave(
  filename = file.path(
    project_dir,
    "Standardized_Posterior_Survival_Stage.png"
  ),
  plot = p_stage,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(
    project_dir,
    "Standardized_Posterior_Survival_Age.png"
  ),
  plot = p_age,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = file.path(
    project_dir,
    "Standardized_Posterior_Survival_Sex.png"
  ),
  plot = p_sex,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 17. FINISHED
# ============================================================

cat("\n\n============================================\n")
cat("EMPIRICAL STANDARDIZATION COMPLETED\n")
cat("============================================\n")

cat("\nStage analysis:\n")
cat("  Stage fixed at I, II, III, IV\n")
cat("  Age + sex marginalized\n")

cat("\nAge analysis:\n")
cat("  Age fixed at 45, 55, 65, 75\n")
cat("  Stage + sex marginalized\n")

cat("\nSex analysis:\n")
cat("  Sex fixed at Female or Male\n")
cat("  Age + stage marginalized\n")

cat("\nSame exact posterior draws used for ALL analyses.\n")
