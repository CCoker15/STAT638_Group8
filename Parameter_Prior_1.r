# ============================================================
# JOINT PRIOR CALIBRATION FOR beta_0 AND k
#
# Weibull model:
#   S(t) = exp(-(t/lambda_0)^k)
#
# Reference patient:
#   Stage I, age 60, female
#
# Goal:
#   Choose beta_0 and k priors jointly so that
#   most prior-predictive survival probabilities fall
#   inside externally plausible survival ranges.
# ============================================================


library(ggplot2)

set.seed(12345)


# ============================================================
# 1. PROJECT DIRECTORY
# ============================================================

project_dir <-
  "D:/OneDrive - Texas A&M University/PhD/Courses/STAT 638/Project"


# ============================================================
# 2. NUMBER OF PRIOR DRAWS
# ============================================================

N <- 10000


# ============================================================
# 3. EXTERNALLY PLAUSIBLE SURVIVAL RANGES
#
# These apply to the reference profile:
# Stage I, age 60, female
# ============================================================

targets <- data.frame(

  Time = c(
    1,
    5,
    10
  ),

  Lower = c(
    0.90,
    0.80,
    0.55
  ),

  Upper = c(
    0.99,
    0.95,
    0.90
  )
)


print(targets)


# ============================================================
# 4. REQUIRED PRIOR-PREDICTIVE COVERAGE
#
# Require at least 80% of prior draws at each time point
# to fall inside its plausible survival range.
# ============================================================

required_coverage <- 0.80


# ============================================================
# 5. CANDIDATE beta_0 PRIORS
#
# beta_0 = log(lambda_0)
#
# Lambda_center is the median prior value of lambda_0
# ============================================================

lambda_centers <- c(
  10,
  12,
  15,
  18,
  20,
  25,
  30
)


beta0_sds <- c(
  0.20,
  0.30,
  0.40,
  0.50
)


# ============================================================
# 6. CANDIDATE k PRIORS
#
# log(k) ~ Normal(log(k_center), logk_sd^2)
# ============================================================

k_centers <- c(
  0.8,
  1.0,
  1.2,
  1.5,
  1.8,
  2.0
)


logk_sds <- c(
  0.15,
  0.20,
  0.30,
  0.40
)


# ============================================================
# 7. FUNCTION TO EVALUATE ONE PRIOR COMBINATION
# ============================================================

evaluate_prior <- function(
    lambda_center,
    beta0_sd,
    k_center,
    logk_sd) {


  # ----------------------------------------------------------
  # Draw beta_0
  # ----------------------------------------------------------

  mu_beta0 <- log(lambda_center)

  beta0 <- rnorm(
    N,
    mean = mu_beta0,
    sd = beta0_sd
  )


  lambda0 <- exp(beta0)


  # ----------------------------------------------------------
  # Draw k
  # ----------------------------------------------------------

  mu_logk <- log(k_center)

  logk <- rnorm(
    N,
    mean = mu_logk,
    sd = logk_sd
  )


  k <- exp(logk)


  # ----------------------------------------------------------
  # Containers
  # ----------------------------------------------------------

  medians <- numeric(nrow(targets))

  lower95 <- numeric(nrow(targets))

  upper95 <- numeric(nrow(targets))

  coverage <- numeric(nrow(targets))


  # ----------------------------------------------------------
  # Prior predictive calculation
  # ----------------------------------------------------------

  for (j in seq_len(nrow(targets))) {


    t <- targets$Time[j]


    S <- exp(
      -(t / lambda0)^k
    )


    # Median survival
    medians[j] <- median(S)


    # 95% prior predictive interval
    lower95[j] <- quantile(
      S,
      probs = 0.025,
      names = FALSE
    )


    upper95[j] <- quantile(
      S,
      probs = 0.975,
      names = FALSE
    )


    # --------------------------------------------------------
    # MAIN QUANTITY:
    #
    # Fraction of prior predictions falling inside
    # the plausible external survival range
    # --------------------------------------------------------

    coverage[j] <- mean(
      S >= targets$Lower[j] &
      S <= targets$Upper[j]
    )
  }


  # ==========================================================
  # 8. SCORE
  #
  # No penalty if coverage >= 80%.
  #
  # If coverage < 80%, penalize the shortfall.
  # ==========================================================

  coverage_deficit <-
    pmax(
      required_coverage - coverage,
      0
    )


  coverage_penalty <-
    sum(
      coverage_deficit^2
    )


  # ----------------------------------------------------------
  # Also require median survival to lie inside target range.
  # ----------------------------------------------------------

  median_penalty <- 0


  for (j in seq_len(nrow(targets))) {


    if (medians[j] < targets$Lower[j]) {

      median_penalty <-
        median_penalty +
        (targets$Lower[j] - medians[j])^2
    }


    if (medians[j] > targets$Upper[j]) {

      median_penalty <-
        median_penalty +
        (medians[j] - targets$Upper[j])^2
    }
  }


  # ----------------------------------------------------------
  # Total score
  # ----------------------------------------------------------

  total_score <-
    coverage_penalty +
    10 * median_penalty


  # ----------------------------------------------------------
  # Did this candidate satisfy all three coverage requirements?
  # ----------------------------------------------------------

  passes_all <-
    all(
      coverage >= required_coverage
    )


  # ----------------------------------------------------------
  # Return result
  # ----------------------------------------------------------

  data.frame(

    Lambda_center =
      lambda_center,

    Beta0_mean =
      mu_beta0,

    Beta0_SD =
      beta0_sd,

    K_center =
      k_center,

    LogK_mean =
      mu_logk,

    LogK_SD =
      logk_sd,


    S1_median =
      medians[1],

    S5_median =
      medians[2],

    S10_median =
      medians[3],


    S1_lower95 =
      lower95[1],

    S1_upper95 =
      upper95[1],


    S5_lower95 =
      lower95[2],

    S5_upper95 =
      upper95[2],


    S10_lower95 =
      lower95[3],

    S10_upper95 =
      upper95[3],


    Coverage_1yr =
      coverage[1],

    Coverage_5yr =
      coverage[2],

    Coverage_10yr =
      coverage[3],


    Minimum_coverage =
      min(coverage),


    Mean_coverage =
      mean(coverage),


    Passes_all =
      passes_all,


    Score =
      total_score
  )
}



# ============================================================
# 9. RUN GRID SEARCH
# ============================================================

results_list <- list()

counter <- 1


for (lc in lambda_centers) {

  for (b_sd in beta0_sds) {

    for (kc in k_centers) {

      for (k_sd in logk_sds) {


        results_list[[counter]] <-
          evaluate_prior(

            lambda_center =
              lc,

            beta0_sd =
              b_sd,

            k_center =
              kc,

            logk_sd =
              k_sd
          )


        counter <-
          counter + 1
      }
    }
  }
}


prior_results <-
  do.call(
    rbind,
    results_list
  )



# ============================================================
# 10. RANK RESULTS
#
# Priority:
#
# 1. Pass all 80% requirements
# 2. Highest minimum coverage
# 3. Highest mean coverage
# 4. Lowest score
# ============================================================

prior_results <-
  prior_results[
    order(
      -prior_results$Passes_all,
      -prior_results$Minimum_coverage,
      -prior_results$Mean_coverage,
      prior_results$Score
    ),
  ]


row.names(prior_results) <- NULL



# ============================================================
# 11. DISPLAY TOP RESULTS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "TOP 20 JOINT PRIOR CANDIDATES\n"
)

cat(
  "============================================\n\n"
)


print(
  head(
    prior_results,
    20
  ),
  digits = 4
)



# ============================================================
# 12. HOW MANY CANDIDATES PASSED?
# ============================================================

n_pass <-
  sum(
    prior_results$Passes_all
  )


cat(
  "\nNumber of candidates meeting >= ",
  required_coverage * 100,
  "% coverage at ALL three times: ",
  n_pass,
  "\n",
  sep = ""
)



# ============================================================
# 13. SELECT BEST CANDIDATE
# ============================================================

best <-
  prior_results[1, ]


cat(
  "\n============================================\n"
)

cat(
  "BEST CANDIDATE\n"
)

cat(
  "============================================\n\n"
)


print(
  best,
  digits = 4
)



# ============================================================
# 14. PRINT PRIORS CLEARLY
# ============================================================

cat(
  "\nSelected beta_0 prior:\n"
)

cat(
  "beta_0 ~ Normal(",
  round(best$Beta0_mean, 3),
  ", ",
  round(best$Beta0_SD, 3),
  "^2)\n",
  sep = ""
)



cat(
  "\nSelected log(k) prior:\n"
)

cat(
  "log(k) ~ Normal(",
  round(best$LogK_mean, 3),
  ", ",
  round(best$LogK_SD, 3),
  "^2)\n",
  sep = ""
)



cat(
  "\nMedian prior lambda_0 = ",
  best$Lambda_center,
  " years\n",
  sep = ""
)


cat(
  "Median prior k = ",
  best$K_center,
  "\n",
  sep = ""
)



# ============================================================
# 15. PRINT COVERAGE IN PERCENTAGES
# ============================================================

cat(
  "\nPrior-predictive coverage:\n"
)


cat(
  "1 year : ",
  round(
    100 * best$Coverage_1yr,
    1
  ),
  "% inside [0.90, 0.99]\n",
  sep = ""
)


cat(
  "5 years: ",
  round(
    100 * best$Coverage_5yr,
    1
  ),
  "% inside [0.80, 0.95]\n",
  sep = ""
)


cat(
  "10 years: ",
  round(
    100 * best$Coverage_10yr,
    1
  ),
  "% inside [0.55, 0.90]\n",
  sep = ""
)



# ============================================================
# 16. GENERATE FINAL PRIOR PREDICTIVE CURVE
# ============================================================

set.seed(9876)


beta0_final <- rnorm(

  N,

  mean =
    best$Beta0_mean,

  sd =
    best$Beta0_SD
)


lambda_final <-
  exp(beta0_final)



logk_final <- rnorm(

  N,

  mean =
    best$LogK_mean,

  sd =
    best$LogK_SD
)


k_final <-
  exp(logk_final)



# ============================================================
# 17. TIME GRID
# ============================================================

time_grid <-
  seq(
    0,
    10,
    by = 0.05
  )



# ============================================================
# 18. PRIOR PREDICTIVE SURVIVAL MATRIX
# ============================================================

S_matrix <- sapply(

  time_grid,

  function(t) {

    exp(
      -(t / lambda_final)^k_final
    )
  }
)



# ============================================================
# 19. SUMMARIZE CURVES
# ============================================================

curve_data <- data.frame(

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



# ============================================================
# 20. PLOT
# ============================================================

p <- ggplot(
  curve_data,
  aes(
    x = Time,
    y = Median
  )
) +


  # 95% prior predictive interval
  geom_ribbon(

    aes(
      ymin = Lower95,
      ymax = Upper95
    ),

    alpha = 0.20
  ) +


  # Median prior predictive curve
  geom_line(
    linewidth = 1.1
  ) +


  # Plausible survival ranges
  geom_segment(

    data = targets,

    aes(
      x = Time,
      xend = Time,
      y = Lower,
      yend = Upper
    ),

    inherit.aes = FALSE,

    linewidth = 6,

    alpha = 0.25
  ) +


  scale_y_continuous(

    limits = c(
      0,
      1
    ),

    breaks =
      seq(
        0,
        1,
        by = 0.2
      )
  ) +


  scale_x_continuous(

    breaks =
      0:10
  ) +


  labs(

    title =
      "Joint Prior Predictive Check",

    subtitle =
      paste0(
        "Reference patient: Stage I, age 60, female",
        "\n",
        "Required coverage = ",
        required_coverage * 100,
        "%"
      ),

    x =
      "Years after diagnosis",

    y =
      "Survival probability"
  ) +


  theme_classic(
    base_family =
      "Times New Roman"
  ) +


  theme(

    text =
      element_text(
        family =
          "Times New Roman",
        size = 12
      ),

    plot.title =
      element_text(
        family =
          "Times New Roman",
        face = "bold",
        hjust = 0.5,
        size = 17
      ),

    plot.subtitle =
      element_text(
        family =
          "Times New Roman",
        hjust = 0.5,
        size = 12
      ),

    axis.title =
      element_text(
        family =
          "Times New Roman",
        size = 13
      ),

    axis.text =
      element_text(
        family =
          "Times New Roman",
        size = 11
      )
  )



print(p)



# ============================================================
# 21. SAVE PLOT
# ============================================================

ggsave(

  filename =
    file.path(
      project_dir,
      "Final_Beta0_K_Prior_80pct_Coverage.png"
    ),

  plot =
    p,

  width =
    8,

  height =
    6,

  dpi =
    300
)



# ============================================================
# 22. SAVE FULL GRID SEARCH
# ============================================================

write.csv(

  prior_results,

  file.path(
    project_dir,
    "Beta0_K_Prior_80pct_Grid_Search.csv"
  ),

  row.names =
    FALSE
)



cat(
  "\nFiles saved successfully.\n"
)