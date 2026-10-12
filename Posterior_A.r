# ============================================================
# BAYESIAN WEIBULL SURVIVAL MODEL - KIRC
#
# POSTERIOR STUDY 1
#
# Fixed prediction profile:
#   Age = 60 years
#   Sex = Female
#
# Compare:
#   Stage I
#   Stage II
#   Stage III
#   Stage IV
#
# MCMC:
#   4 chains
#   2000 warm-up per chain
#   2000 retained per chain
#
# Total retained posterior draws = 8000
# ============================================================


library(readxl)
library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)


set.seed(638)


# ============================================================
# 1. PROJECT DIRECTORY
# ============================================================

project_dir <-
  "D:/OneDrive - Texas A&M University/PhD/Courses/STAT 638/Project"


data_file <-
  file.path(
    project_dir,
    "KIRC - Filtered Data.xlsx"
  )


if (!file.exists(data_file)) {
  stop(
    "KIRC data file not found:\n",
    data_file
  )
}


# ============================================================
# 2. READ KIRC DATA
# ============================================================

dat <- read_excel(data_file)


cat("\n========================================\n")
cat("KIRC DATA\n")
cat("========================================\n")

cat("\nNumber of rows =", nrow(dat), "\n")

cat("\nColumn names:\n")
print(names(dat))


# ============================================================
# 3. IDENTIFY AGE AND SEX COLUMNS
#
# The survival/stage columns are already known:
#
# overall_survival_time
# overall_survival
# pathologic_stage
#
# This block searches common age/sex names.
# ============================================================

find_first_column <- function(data, candidates) {

  found <- candidates[candidates %in% names(data)]

  if (length(found) == 0) {
    return(NA_character_)
  }

  found[1]
}


AGE_COL <- find_first_column(
  dat,
  c(
    "age_at_diagnosis",
    "age_at_index",
    "age",
    "Age",
    "AGE"
  )
)


SEX_COL <- find_first_column(
  dat,
  c(
    "gender",
    "sex",
    "Gender",
    "Sex",
    "SEX"
  )
)


if (is.na(AGE_COL)) {

  stop(
    "\nAge column could not be detected.\n",
    "Available columns are:\n",
    paste(names(dat), collapse = ", ")
  )
}


if (is.na(SEX_COL)) {

  stop(
    "\nSex column could not be detected.\n",
    "Available columns are:\n",
    paste(names(dat), collapse = ", ")
  )
}


cat("\nAge column detected:", AGE_COL, "\n")
cat("Sex column detected:", SEX_COL, "\n")


# ============================================================
# 4. PREPARE SURVIVAL TIME
#
# overall_survival_time is in DAYS.
#
# Priors for beta0/lambda were calibrated using YEARS,
# therefore the likelihood must also use years.
# ============================================================

dat$OS_years <-
  as.numeric(dat$overall_survival_time) / 365


# Event:
# 1 = death
# 0 = censored

dat$event <-
  as.integer(dat$overall_survival)


# ============================================================
# 5. PREPARE AGE
# ============================================================

age_raw <-
  as.numeric(dat[[AGE_COL]])


# TCGA age_at_diagnosis is sometimes stored in days.
# If values are clearly > human age in years,
# convert to years.

if (median(age_raw, na.rm = TRUE) > 150) {

  cat(
    "\nAge appears to be stored in days.",
    "\nConverting age to years.\n"
  )

  dat$Age_years <-
    age_raw / 365.25

} else {

  dat$Age_years <-
    age_raw
}


# Age centered at 60 and scaled per 10 years

dat$Age10 <-
  (dat$Age_years - 60) / 10


cat(
  "\nAge range in years:",
  round(min(dat$Age_years, na.rm = TRUE), 1),
  "to",
  round(max(dat$Age_years, na.rm = TRUE), 1),
  "\n"
)


# ============================================================
# 6. PREPARE SEX
#
# Female = 0
# Male   = 1
# ============================================================

sex_text <-
  tolower(
    trimws(
      as.character(dat[[SEX_COL]])
    )
  )


dat$Male <- NA_integer_


dat$Male[
  sex_text %in% c(
    "female",
    "f"
  )
] <- 0


dat$Male[
  sex_text %in% c(
    "male",
    "m"
  )
] <- 1


if (any(is.na(dat$Male))) {

  cat(
    "\nUnrecognized sex values:\n"
  )

  print(
    unique(
      dat[[SEX_COL]][is.na(dat$Male)]
    )
  )

  stop(
    "Sex variable contains values that were not recognized."
  )
}


cat("\nSex distribution:\n")
print(table(dat$Male))

cat("\n0 = Female, 1 = Male\n")


# ============================================================
# 7. PREPARE PATHOLOGIC STAGE
#
# Stage I is reference.
# ============================================================

stage_text <-
  trimws(
    as.character(dat$pathologic_stage)
  )


dat$StageII <-
  as.integer(stage_text == "Stage II")


dat$StageIII <-
  as.integer(stage_text == "Stage III")


dat$StageIV <-
  as.integer(stage_text == "Stage IV")


valid_stage <-
  stage_text %in%
  c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )


if (!all(valid_stage)) {

  cat(
    "\nUnexpected pathologic stage values:\n"
  )

  print(
    unique(
      stage_text[!valid_stage]
    )
  )

  stop(
    "Unexpected stage categories found."
  )
}


cat("\nStage distribution:\n")
print(table(stage_text))


# ============================================================
# 8. CHECK DATA BEFORE MODEL FITTING
# ============================================================

if (any(dat$OS_years <= 0, na.rm = TRUE)) {

  stop(
    "Survival time must be > 0 for Weibull likelihood."
  )
}


if (!all(dat$event %in% c(0, 1))) {

  stop(
    "overall_survival must contain only 0 and 1."
  )
}


required_complete <- complete.cases(
  dat$OS_years,
  dat$event,
  dat$StageII,
  dat$StageIII,
  dat$StageIV,
  dat$Age10,
  dat$Male
)


if (!all(required_complete)) {

  stop(
    sum(!required_complete),
    " rows contain missing model variables."
  )
}


N <- nrow(dat)


cat("\nPatients entering Bayesian model =", N, "\n")
cat("Deaths =", sum(dat$event == 1), "\n")
cat("Censored =", sum(dat$event == 0), "\n")


# ============================================================
# 9. STAN DATA
# ============================================================

stan_data <- list(

  N = N,

  time =
    as.numeric(dat$OS_years),

  event =
    as.integer(dat$event),

  stage2 =
    as.numeric(dat$StageII),

  stage3 =
    as.numeric(dat$StageIII),

  stage4 =
    as.numeric(dat$StageIV),

  age10 =
    as.numeric(dat$Age10),

  male =
    as.numeric(dat$Male)
)


# ============================================================
# 10. STAN MODEL
#
# Model:
#
# log(lambda_i) =
#   beta0
# + beta2 * StageII
# + beta3 * StageIII
# + beta4 * StageIV
# + betaA * Age10
# + betaS * Male
#
# T_i ~ Weibull(k, lambda_i)
#
# Death:
#   f(t_i)
#
# Censored:
#   S(t_i)
# ============================================================

stan_code <- '

data {

  int<lower=1> N;

  vector<lower=0>[N] time;

  array[N] int<lower=0, upper=1> event;

  vector[N] stage2;
  vector[N] stage3;
  vector[N] stage4;

  vector[N] age10;

  vector[N] male;
}


parameters {

  real beta0;

  real beta2;
  real beta3;
  real beta4;

  real betaA;
  real betaS;

  real log_k;
}


transformed parameters {

  real<lower=0> k;

  vector[N] eta;

  vector<lower=0>[N] lambda;


  k = exp(log_k);


  eta =
      beta0
    + beta2 * stage2
    + beta3 * stage3
    + beta4 * stage4
    + betaA * age10
    + betaS * male;


  lambda = exp(eta);
}


model {

  // ------------------------------------------
  // PRIORS
  // ------------------------------------------

  beta0 ~ normal(3.219, 0.200);

  beta2 ~ normal(-0.343, 0.278);

  beta3 ~ normal(-0.773, 0.307);

  beta4 ~ normal(-1.939, 0.496);

  betaA ~ normal(-0.330, 0.098);

  betaS ~ normal(0.000, 0.175);

  log_k ~ normal(0.182, 0.150);


  // ------------------------------------------
  // CENSORED WEIBULL LIKELIHOOD
  // ------------------------------------------

  for (i in 1:N) {

    if (event[i] == 1) {

      // Death observed
      target +=
        weibull_lpdf(
          time[i] |
          k,
          lambda[i]
        );

    } else {

      // Right censored:
      // probability of surviving beyond time[i]

      target +=
        weibull_lccdf(
          time[i] |
          k,
          lambda[i]
        );
    }
  }
}
'


# ============================================================
# 11. WRITE STAN MODEL TO FILE
# ============================================================

stan_file <-
  file.path(
    project_dir,
    "KIRC_Bayesian_Weibull_Study1.stan"
  )


writeLines(
  stan_code,
  stan_file
)


cat(
  "\nStan model written to:\n",
  stan_file,
  "\n"
)


# ============================================================
# 12. COMPILE STAN MODEL
# ============================================================

model <-
  cmdstan_model(
    stan_file
  )


# ============================================================
# 13. RUN MCMC
#
# 4 chains
#
# EACH chain:
#   2000 warm-up
#   2000 retained
#
# Total retained:
#
#   4 x 2000 = 8000 posterior draws
# ============================================================

fit <-
  model$sample(

    data =
      stan_data,

    seed =
      638,

    chains =
      4,

    parallel_chains =
      min(
        4,
        parallel::detectCores()
      ),

    iter_warmup =
      2000,

    iter_sampling =
      2000,

    adapt_delta =
      0.95,

    max_treedepth =
      12,

    refresh =
      100
  )


# ============================================================
# 14. SAVE FIT
# ============================================================

fit$save_object(

  file =
    file.path(
      project_dir,
      "KIRC_Bayesian_Weibull_Fit.rds"
    )
)


# ============================================================
# 15. MCMC DIAGNOSTICS
# ============================================================

parameters <-
  c(
    "beta0",
    "beta2",
    "beta3",
    "beta4",
    "betaA",
    "betaS",
    "log_k",
    "k"
  )


diagnostics <-
  fit$summary(
    variables =
      parameters
  )


cat(
  "\n========================================\n"
)

cat(
  "POSTERIOR + MCMC DIAGNOSTICS\n"
)

cat(
  "========================================\n\n"
)


print(
  diagnostics
)


cat(
  "\nSampler diagnostic summary:\n"
)

print(
  fit$diagnostic_summary()
)


# ============================================================
# 16. EXTRACT ALL 8000 POSTERIOR DRAWS
# ============================================================

draws <-
  as_draws_df(
    fit$draws(
      variables =
        parameters
    )
  )


posterior_draws <-
  as.data.frame(draws)


cat(
  "\nNumber of retained posterior draws =",
  nrow(posterior_draws),
  "\n"
)


# Save exact posterior draws

write.csv(

  posterior_draws,

  file.path(
    project_dir,
    "KIRC_Weibull_Posterior_Draws.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 17. POSTERIOR PARAMETER SUMMARY
# ============================================================

make_parameter_summary <- function(
    posterior_draws,
    parameter_names) {

  output <- lapply(
    parameter_names,
    function(v) {

      x <- posterior_draws[[v]]

      data.frame(

        Parameter =
          v,

        Mean =
          mean(x),

        SD =
          sd(x),

        Median =
          median(x),

        Lower_95_CrI =
          quantile(
            x,
            0.025,
            names = FALSE
          ),

        Upper_95_CrI =
          quantile(
            x,
            0.975,
            names = FALSE
          )
      )
    }
  )

  do.call(
    rbind,
    output
  )
}


posterior_parameter_summary <-
  make_parameter_summary(
    posterior_draws,
    parameters
  )


cat(
  "\n========================================\n"
)

cat(
  "POSTERIOR PARAMETER SUMMARY\n"
)

cat(
  "========================================\n\n"
)


print(
  posterior_parameter_summary,
  digits = 4
)


write.csv(

  posterior_parameter_summary,

  file.path(
    project_dir,
    "KIRC_Weibull_Posterior_Parameter_Summary.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 18. STUDY 1
#
# FIXED REFERENCE PROFILE:
#
# Age = 60
# Female
#
# Therefore:
#
# Age10 = 0
# Male  = 0
#
# We now reuse the EXACT SAME posterior rows.
# ============================================================


# Time grid for complete survival curves

time_grid <-
  seq(
    0,
    10,
    by = 0.05
  )


stage_names <-
  c(
    "Stage I",
    "Stage II",
    "Stage III",
    "Stage IV"
  )


# ============================================================
# 19. FUNCTION TO CALCULATE SURVIVAL FOR ONE STAGE
# ============================================================

get_eta_for_stage <- function(
    posterior_draws,
    stage) {

  if (stage == "Stage I") {

    eta <-
      posterior_draws$beta0

  } else if (stage == "Stage II") {

    eta <-
      posterior_draws$beta0 +
      posterior_draws$beta2

  } else if (stage == "Stage III") {

    eta <-
      posterior_draws$beta0 +
      posterior_draws$beta3

  } else if (stage == "Stage IV") {

    eta <-
      posterior_draws$beta0 +
      posterior_draws$beta4

  } else {

    stop(
      "Unknown stage."
    )
  }


  return(eta)
}


# ============================================================
# 20. CALCULATE POSTERIOR SURVIVAL CURVES
# ============================================================

curve_list <- list()


for (stage in stage_names) {


  # ------------------------------------------
  # For every posterior row:
  #
  # lambda = exp(eta)
  #
  # IMPORTANT:
  # k and beta values come from SAME row.
  # ------------------------------------------

  eta_draws <-
    get_eta_for_stage(
      posterior_draws,
      stage
    )


  lambda_draws <-
    exp(
      eta_draws
    )


  k_draws <-
    posterior_draws$k


  # ------------------------------------------
  # 8000 x number-of-time-points matrix
  # ------------------------------------------

  S_matrix <-
    sapply(
      time_grid,
      function(t) {

        exp(
          -(
            t /
            lambda_draws
          )^k_draws
        )
      }
    )


  # ------------------------------------------
  # Summarize the 8000 survival estimates
  # at every time point
  # ------------------------------------------

  stage_curve <-
    data.frame(

      Stage =
        stage,

      Time =
        time_grid,

      Median =
        apply(
          S_matrix,
          2,
          median
        ),

      Lower95 =
        apply(
          S_matrix,
          2,
          quantile,
          probs = 0.025
        ),

      Upper95 =
        apply(
          S_matrix,
          2,
          quantile,
          probs = 0.975
        )
    )


  curve_list[[stage]] <-
    stage_curve
}


posterior_curves <-
  bind_rows(
    curve_list
  )


# ============================================================
# 21. SAVE COMPLETE POSTERIOR CURVES
# ============================================================

write.csv(

  posterior_curves,

  file.path(
    project_dir,
    "Study1_Posterior_Survival_Curves.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 22. POSTERIOR SURVIVAL AT 1, 5, 10 YEARS
# ============================================================

times_interest <-
  c(
    1,
    5,
    10
  )


survival_table_list <- list()

counter <- 1


for (stage in stage_names) {


  eta_draws <-
    get_eta_for_stage(
      posterior_draws,
      stage
    )


  lambda_draws <-
    exp(
      eta_draws
    )


  k_draws <-
    posterior_draws$k


  for (t in times_interest) {


    S_draws <-
      exp(
        -(
          t /
          lambda_draws
        )^k_draws
      )


    survival_table_list[[counter]] <-
      data.frame(

        Stage =
          stage,

        Time_years =
          t,

        Posterior_mean =
          mean(
            S_draws
          ),

        Posterior_median =
          median(
            S_draws
          ),

        Lower_95_CrI =
          quantile(
            S_draws,
            0.025,
            names = FALSE
          ),

        Upper_95_CrI =
          quantile(
            S_draws,
            0.975,
            names = FALSE
          )
      )


    counter <-
      counter + 1
  }
}


survival_table <-
  bind_rows(
    survival_table_list
  )


# Convert to percentage

survival_table <-
  survival_table %>%

  mutate(

    Median_percent =
      100 *
      Posterior_median,

    Lower_95_percent =
      100 *
      Lower_95_CrI,

    Upper_95_percent =
      100 *
      Upper_95_CrI
  )


cat(
  "\n========================================\n"
)

cat(
  "STUDY 1: POSTERIOR SURVIVAL\n"
)

cat(
  "Reference patient: Age 60, Female\n"
)

cat(
  "========================================\n\n"
)


print(
  survival_table,
  digits = 4
)


write.csv(

  survival_table,

  file.path(
    project_dir,
    "Study1_Posterior_Survival_1_5_10_years.csv"
  ),

  row.names = FALSE
)


# ============================================================
# 23. POSTERIOR SURVIVAL CURVE FIGURE
# ============================================================

font_name <-
  "Times New Roman"


p_stage <-
  ggplot(

    posterior_curves,

    aes(
      x = Time,
      y = Median,
      group = Stage
    )

  ) +


  # 95% credible bands
  geom_ribbon(

    aes(
      ymin = Lower95,
      ymax = Upper95,
      fill = Stage
    ),

    alpha =
      0.12,

    colour =
      NA
  ) +


  # Posterior median curves
  geom_line(

    aes(
      colour = Stage
    ),

    linewidth =
      1
  ) +


  scale_y_continuous(

    limits =
      c(0, 1),

    breaks =
      seq(
        0,
        1,
        by = 0.2
      ),

    labels =
      function(x) {
        paste0(
          round(
            100 * x
          ),
          "%"
        )
      }
  ) +


  scale_x_continuous(

    limits =
      c(0, 10),

    breaks =
      0:10
  ) +


  labs(

    title =
      "Posterior Survival by Pathologic Stage",

    subtitle =
      "Reference patient: 60-year-old female",

    x =
      "Years after diagnosis",

    y =
      "Posterior survival probability",

    colour =
      "Pathologic Stage",

    fill =
      "Pathologic Stage"
  ) +


  theme_classic(
    base_family =
      font_name
  ) +


  theme(

    text =
      element_text(
        family =
          font_name,
        size =
          12
      ),

    plot.title =
      element_text(
        family =
          font_name,
        size =
          15,
        face =
          "bold",
        hjust =
          0.5
      ),

    plot.subtitle =
      element_text(
        family =
          font_name,
        size =
          12,
        hjust =
          0.5
      ),

    axis.title =
      element_text(
        family =
          font_name,
        size =
          12
      ),

    axis.text =
      element_text(
        family =
          font_name,
        size =
          11
      ),

    legend.title =
      element_text(
        family =
          font_name,
        size =
          11
      ),

    legend.text =
      element_text(
        family =
          font_name,
        size =
          10
      )
  )


print(
  p_stage
)


# ============================================================
# 24. SAVE FIGURE
# ============================================================

ggsave(

  filename =
    file.path(
      project_dir,
      "Study1_Posterior_Survival_by_Stage.png"
    ),

  plot =
    p_stage,

  width =
    8,

  height =
    6,

  dpi =
    300
)


# ============================================================
# 25. FINISHED
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "ANALYSIS COMPLETED\n"
)

cat(
  "========================================\n"
)


cat(
  "\nImportant outputs:\n"
)


cat(
  "\n1. KIRC_Weibull_Posterior_Draws.csv"
)

cat(
  "\n   -> all 8000 exact joint posterior draws\n"
)


cat(
  "\n2. KIRC_Weibull_Posterior_Parameter_Summary.csv"
)

cat(
  "\n   -> posterior estimates of beta and k\n"
)


cat(
  "\n3. Study1_Posterior_Survival_Curves.csv"
)

cat(
  "\n   -> complete Stage I-IV posterior curves\n"
)


cat(
  "\n4. Study1_Posterior_Survival_1_5_10_years.csv"
)

cat(
  "\n   -> numerical survival estimates at 1, 5, and 10 years\n"
)


cat(
  "\n5. Study1_Posterior_Survival_by_Stage.png"
)

cat(
  "\n   -> final Study 1 figure\n"
)