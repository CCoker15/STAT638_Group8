###############################################################################
# STAT 638 | TCGA-KIRC | Bayesian log-normal AFT survival analysis
# Complete, self-contained analysis aligned with Weibull / piecewise-exponential
# 521 patients; stage, age per decade, sex, grade; informative vs vague priors
# Requires: KIRC_filtered_521_for_models.xlsx (select interactively if needed)
# Outputs are saved to Lognormal_Results/ beside the selected workbook.
# NOTE: This script has not been executed against an R/Stan runtime here.
###############################################################################

# STEP 0: Packages (install missing packages once)
packages <- c('readxl','dplyr','survival','brms','posterior','loo','ggplot2')
missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install missing packages first: ',paste(missing,collapse=', '),
                        '\nRun: install.packages(c(',paste(sprintf('"%s"',missing),collapse=', '),'))')
suppressPackageStartupMessages({
  library(readxl);library(dplyr);library(survival);library(brms)
  library(posterior);library(loo);library(ggplot2)
})
options(mc.cores = min(4L, parallel::detectCores()))
set.seed(638)

# STEP 1: Load data and harmonize cohort and order
input_file <- 'KIRC_filtered_521_for_models.xlsx'
if(!file.exists(input_file)) {
  message('Select KIRC_filtered_521_for_models.xlsx in the file chooser')
  input_file <- file.choose()
}
input_file <- normalizePath(input_file)
out_dir <- file.path(dirname(input_file),'Lognormal_Results')
dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
raw <- readxl::read_excel(input_file)
needed <- c('bcr_patient_barcode','overall_survival_time','overall_survival',
            'pathologic_stage','age_at_diagnosis','gender','histologic_grade')
stopifnot(all(needed %in% names(raw)))
dat <- raw %>%
  mutate(grade_raw=toupper(trimws(as.character(histologic_grade)))) %>%
  filter(grade_raw %in% c('G1','G2','G3','G4')) %>%
  arrange(bcr_patient_barcode) %>%
  mutate(OS_years=as.numeric(overall_survival_time)/365,
         event=as.integer(overall_survival),
         censored=1L-event, # brms: 0=death observed; 1=right censored
         Age10=(as.numeric(age_at_diagnosis)-60)/10,
         Male=as.integer(toupper(trimws(as.character(gender)))=='MALE'),
         stage=factor(sub('^Stage\\s+','',trimws(as.character(pathologic_stage))),
                      levels=c('I','II','III','IV')),
         grade=factor(ifelse(grade_raw %in% c('G1','G2'),'G1-G2',grade_raw),
                      levels=c('G1-G2','G3','G4')),
         StageII=as.integer(stage=='II'),StageIII=as.integer(stage=='III'),
         StageIV=as.integer(stage=='IV'),
         GradeG3=as.integer(grade=='G3'),GradeG4=as.integer(grade=='G4'))
stopifnot(nrow(dat)==521L,sum(dat$event)==173L,sum(dat$censored)==348L,
          all(is.finite(dat$OS_years)),all(dat$OS_years>0),
          all(!is.na(dat$stage)),all(!is.na(dat$grade)),
          all(complete.cases(dat[,c('Age10','Male')])),
          all(dat$event %in% c(0L,1L)))
write.csv(dat[,c('bcr_patient_barcode','OS_years','event','censored','stage','Age10','Male','grade')],
          file.path(out_dir,'Lognormal_Analysis_Cohort.csv'),row.names=FALSE)
cat('N =',nrow(dat),' deaths =',sum(dat$event),' censored =',sum(dat$censored),'\n')
print(table(dat$stage,dat$event));print(table(dat$grade,dat$event))

# STEP 2: Log-normal AFT model
# log(T_i) ~ Normal(mu_i,sigma^2), mu_i = b0 + X_i*b
# Reference: Stage I, age 60, female, grade G1-G2
f_ln <- brms::bf(OS_years | cens(censored) ~ StageII+StageIII+StageIV+
                  Age10+Male+GradeG3+GradeG4)
# Death likelihood f(t); censored likelihood S(t).

# STEP 3: Shared literature-based HR evidence (from PWE_Prior_Derivation.R)
# HR centers and plausible ranges are shared with the group's other models.
hr <- data.frame(
  term=c('StageII','StageIII','StageIV','Age10','Male','GradeG3','GradeG4'),
  center=c(1.5,2.5,10,1.04,1,1.4,2.5),
  lower=c(.8,1.3,4,1.02,.67,.9,1.4),
  upper=c(2.8,4.8,25,1.06,1.5,2.2,4.5),
  multiplier=c(1,1,1,10,1,1,1))
hr$logHR_mu <- hr$multiplier*log(hr$center)
hr$logHR_sd <- hr$multiplier*(log(hr$upper)-log(hr$lower))/(2*1.96)
# IMPORTANT: A log-normal AFT model does NOT have constant hazard ratios.
# To transport shared prior INFORMATION, use Weibull prior median k=exp(.182)
# as an approximate AFT calibration: beta_AFT = -log(HR)/k.
# This is an APPROXIMATION, not an exact log-normal HR conversion.
# We propagate uncertainty in k and HR with Monte Carlo and moment matching.
set.seed(638)
n_prior <- 200000L
k_draw <- exp(rnorm(n_prior,mean=.182,sd=.15))
prior_cov <- do.call(rbind,lapply(seq_len(nrow(hr)),function(j){
  h <- rnorm(n_prior,hr$logHR_mu[j],hr$logHR_sd[j])
  b <- -h/k_draw
  data.frame(term=hr$term[j],mean=mean(b),sd=sd(b),
             q025=unname(quantile(b,.025)),q975=unname(quantile(b,.975)))
}))
print(prior_cov)
write.csv(prior_cov,file.path(out_dir,'Lognormal_Prior_Coefficients.csv'),row.names=FALSE)

# ============================================================
# STEP 4: BASELINE LOG-NORMAL PRIOR CALIBRATION
# ============================================================

# Purpose:
# Calibrate the log-normal baseline priors using the same
# Weibull reference-patient prior assumptions as the group.
#
# Weibull baseline priors:
# log(scale) ~ Normal(3.219, 0.20^2)
# log(shape) ~ Normal(0.182, 0.15^2)
#
# Calibration times: 1, 5, and 10 years

set.seed(639)

# Number of prior simulations
n_prior <- 200000

# ------------------------------------------------------------
# 4.1 Generate Weibull prior draws
# ------------------------------------------------------------

weib_mu <- rnorm(
  n_prior,
  mean = 3.219,
  sd = 0.20
)

weib_k <- exp(
  rnorm(
    n_prior,
    mean = 0.182,
    sd = 0.15
  )
)

# ------------------------------------------------------------
# 4.2 Calculate Weibull cumulative probabilities
# ------------------------------------------------------------

cal_times <- c(1, 5, 10)

x <- log(cal_times)

Fmat <- sapply(
  cal_times,
  function(t) {
    -expm1(
      -exp(weib_k * (log(t) - weib_mu))
    )
  }
)

# ------------------------------------------------------------
# 4.3 Convert Weibull probabilities to normal quantiles
# ------------------------------------------------------------

# IMPORTANT:
# qnorm() can remove matrix dimensions.
# matrix() restores the correct dimensions.

z <- matrix(
  qnorm(
    pmin(
      1 - 1e-9,
      pmax(1e-9, Fmat)
    )
  ),
  nrow = nrow(Fmat),
  ncol = ncol(Fmat)
)

# Verify dimensions
cat("Dimensions of Fmat:", dim(Fmat), "\n")
cat("Dimensions of z:", dim(z), "\n")

# ------------------------------------------------------------
# 4.4 Fit log-normal approximation for each prior draw
# ------------------------------------------------------------

# Log-normal CDF:
# F(t) = Phi((log(t) - mu0) / sigma)
#
# Therefore:
# z(t) = a + b * log(t)
#
# b = 1/sigma
# a = -mu0/sigma

xb <- mean(x)

xss <- sum((x - xb)^2)

# Calculate regression slope for each prior draw
slope <- as.vector(
  z %*% (x - xb)
) / xss

# Calculate regression intercept
intercept <- rowMeans(z) - slope * xb

# Convert regression parameters into log-normal parameters
sigma_draw <- 1 / slope

mu_draw <- -intercept / slope

# ------------------------------------------------------------
# 4.5 Remove invalid draws
# ------------------------------------------------------------

valid <- is.finite(mu_draw) &
  is.finite(sigma_draw) &
  sigma_draw > 0

mu_draw <- mu_draw[valid]

sigma_draw <- sigma_draw[valid]

cat("Valid prior draws:", length(mu_draw), "\n")

# ------------------------------------------------------------
# 4.6 Moment-match log-normal prior distributions
# ------------------------------------------------------------

# Intercept prior: mu0 ~ Normal(mean, sd^2)

mu_prior_mean <- mean(mu_draw)

mu_prior_sd <- sd(mu_draw)

# Residual scale prior:
# log(sigma) ~ Normal(mean, sd^2)

logsigma_prior_mean <- mean(log(sigma_draw))

logsigma_prior_sd <- sd(log(sigma_draw))

# ------------------------------------------------------------
# 4.7 Display calibrated prior parameters
# ------------------------------------------------------------

cat("\n====================================\n")
cat("CALIBRATED LOG-NORMAL PRIORS\n")
cat("====================================\n")

cat("\nIntercept (mu0):\n")

print(
  c(
    mean = mu_prior_mean,
    sd = mu_prior_sd
  )
)

cat("\nLog residual scale (log sigma):\n")

print(
  c(
    mean = logsigma_prior_mean,
    sd = logsigma_prior_sd
  )
)

# ------------------------------------------------------------
# 4.8 Compare prior survival probabilities
# ------------------------------------------------------------

# Weibull survival probabilities
S_weib <- sapply(
  cal_times,
  function(t) {
    exp(
      -exp(weib_k * (log(t) - weib_mu))
    )
  }
)

# Log-normal survival probabilities
S_logn <- sapply(
  cal_times,
  function(t) {
    plnorm(
      t,
      meanlog = mu_draw,
      sdlog = sigma_draw,
      lower.tail = FALSE
    )
  }
)

comparison <- data.frame(
  Time_Years = cal_times,
  Weibull_Median = apply(S_weib, 2, median),
  Lognormal_Median = apply(S_logn, 2, median)
)

cat("\nPrior survival comparison:\n")

print(comparison)

cat("\nSTEP 4 COMPLETED SUCCESSFULLY\n")

# Variable names required by subsequent steps

base_mu_mean <- mu_prior_mean
base_mu_sd   <- mu_prior_sd

logsig_mean <- logsigma_prior_mean
logsig_sd   <- logsigma_prior_sd

base_lsigma_mean <- logsigma_prior_mean
base_lsigma_sd   <- logsigma_prior_sd

# STEP 5: Prior predictive calibration (reference patient)
set.seed(640)
mu_test <- rnorm(50000, base_mu_mean, base_mu_sd)
sig_test <- rlnorm(50000, logsig_mean, logsig_sd)
prior_check <- data.frame(Time_years=cal_times,
  Weibull_Prior_Median=vapply(cal_times,function(t)
    median(exp(-exp(weib_k*(log(t)-weib_mu)))),numeric(1)),
  Lognormal_Prior_Median=vapply(cal_times,function(t)
    median(plnorm(t,meanlog=mu_test,sdlog=sig_test,lower.tail=FALSE)),numeric(1)),
  Lognormal_Lower=vapply(cal_times,function(t)
    unname(quantile(plnorm(t,mu_test,sig_test,lower.tail=FALSE),.025)),numeric(1)),
  Lognormal_Upper=vapply(cal_times,function(t)
    unname(quantile(plnorm(t,mu_test,sig_test,lower.tail=FALSE),.975)),numeric(1)))
print(prior_check)
write.csv(prior_check,file.path(out_dir,'Lognormal_Prior_Predictive_Check.csv'),row.names=FALSE)
# Inspect this table. If the two prior predictive curves disagree materially,
# recalibration is needed before treating models as harmonized.

# STEP 6: Specify informative and vague priors
# brms prior on class='Intercept' refers to centered design-matrix intercept,
# NOT the reference patient's intercept. To preserve baseline calibration,
# fit an explicit intercept via '0 + Intercept' and class='b', coef='Intercept'.
f_ln <- brms::bf(OS_years | cens(censored) ~ 0+Intercept+
                  StageII+StageIII+StageIV+Age10+Male+GradeG3+GradeG4)
dat$Intercept <- 1
p_inf <- c(
  set_prior(sprintf('normal(%.8f, %.8f)',base_mu_mean,base_mu_sd),
            class='b',coef='Intercept'),
  set_prior(sprintf('lognormal(%.8f, %.8f)',logsig_mean,logsig_sd),class='sigma'))
for(j in seq_len(nrow(prior_cov)))
  p_inf <- c(p_inf,set_prior(sprintf('normal(%.8f, %.8f)',prior_cov$mean[j],prior_cov$sd[j]),
                              class='b',coef=prior_cov$term[j]))
# Vague priors: not identical numerical priors to PWE (different parameters).
# Broad but proper priors, with positive sigma; check prior predictive plausibility.
p_vague <- c(set_prior('normal(2, 5)',class='b',coef='Intercept'),
             set_prior('normal(0, 5)',class='b'),
             set_prior('lognormal(0, 1)',class='sigma'))
print(p_inf);print(p_vague)
print(brms::get_prior(f_ln,data=dat,family=lognormal()))

# STEP 7: Fit informative and vague models (can take several minutes)
fit_one <- function(priors,tag){
  path <- file.path(out_dir,paste0('Lognormal_Fit_',tag,'.rds'))
  if(file.exists(path)) return(readRDS(path))
  fit <- brms::brm(f_ln,data=dat,family=lognormal(),prior=priors,
                   chains=4,iter=4000,warmup=2000,cores=min(4L,parallel::detectCores()),
                   seed=638,control=list(adapt_delta=.97,max_treedepth=12),
                   save_pars=brms::save_pars(all=TRUE),refresh=200)
  saveRDS(fit,path)
  fit
}
fit_inf <- fit_one(p_inf,'informative')
fit_vague <- fit_one(p_vague,'vague')

# STEP 8: Diagnostics and coefficient summaries
summary(fit_inf);summary(fit_vague)
for(tag in c('informative','vague')){
  fit <- if(tag=='informative')fit_inf else fit_vague
  sm <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(fit)))
  write.csv(sm,file.path(out_dir,paste0('Lognormal_MCMC_Diagnostics_',tag,'.csv')),row.names=FALSE)
  pdf(file.path(out_dir,paste0('Lognormal_Trace_',tag,'.pdf')),
      width=11,height=8,onefile=TRUE)
  plot(fit, ask=FALSE)
  dev.off()
  pdf(file.path(out_dir,paste0('Lognormal_PPCheck_',tag,'.pdf')),width=9,height=6)
  print(pp_check(fit,ndraws=40));dev.off()
  # brms pp_check for censored models may exclude censored responses: interpret accordingly.
}
# Log-normal AFT coefficient exp(beta) = median survival-time ratio, NOT HR.
coef_table <- function(fit,tag){
  d <- posterior::as_draws_df(fit)
  nm <- paste0('b_',prior_cov$term)
  tab <- do.call(rbind,lapply(nm,function(v){
    z <- exp(as.numeric(d[[v]]))
    data.frame(term=sub('^b_','',v),median=median(z),
               lower=unname(quantile(z,.025)),upper=unname(quantile(z,.975)))
  }))
  write.csv(tab,file.path(out_dir,paste0('Lognormal_Time_Ratios_',tag,'.csv')),row.names=FALSE)
  tab
}
print(coef_table(fit_inf,'informative'))
print(coef_table(fit_vague,'vague'))

# STEP 9: Survival prediction helper (posterior draws x patients)
# Using brms posterior_linpred returns mu; sigma drawn separately.
# Restrict draws to control memory / runtime.
get_draw_indices <- function(fit,n=800L){
  total <- posterior::ndraws(fit)
  unique(as.integer(round(seq(1,total,length.out=min(n,total)))))
}
pred_surv <- function(fit,newdata,t,draw_ids=NULL,return_draws=FALSE){
  if(is.null(draw_ids)) draw_ids <- get_draw_indices(fit)
  lp <- brms::posterior_linpred(fit,newdata=newdata,transform=FALSE,
                                draw_ids=draw_ids)
  sig <- as.numeric(posterior::as_draws_df(fit)$sigma)[draw_ids]
  # Log-normal S(t) = Phi((mu-log(t))/sigma)
  if(t<=0) S <- matrix(1,nrow=nrow(lp),ncol=ncol(lp)) else S <- pnorm(sweep(lp-log(t),1,sig,'/'))
  if(return_draws) return(S)
  colMeans(S)
}

# STEP 10: Study A, reference patient stage survival at 1,5,10 years
ref <- data.frame(Intercept=1,StageII=0,StageIII=0,StageIV=0,
                  Age10=0,Male=0,GradeG3=0,GradeG4=0)
make_stage <- function(df,s){
  df$StageII <- as.integer(s=='II');df$StageIII <- as.integer(s=='III')
  df$StageIV <- as.integer(s=='IV');df
}
make_age <- function(df,a){df$Age10 <- (a-60)/10;df}
make_grade <- function(df,g){
  df$GradeG3 <- as.integer(g=='G3');df$GradeG4 <- as.integer(g=='G4');df
}
ref_table <- function(fit,tag){
  ids <- get_draw_indices(fit,2000)
  rows <- list()
  for(sex in c(0,1))for(s in c('I','II','III','IV'))for(t in c(1,5,10)){
    nd <- make_stage(ref,s);nd$Male <- sex
    v <- as.vector(pred_surv(fit,nd,t,ids,TRUE))
    rows[[length(rows)+1]] <- data.frame(sex=ifelse(sex==1,'Male','Female'),
      stage=s,Time_years=t,median=median(v),
      lower=unname(quantile(v,.025)),upper=unname(quantile(v,.975)))
  }
  tab <- dplyr::bind_rows(rows)
  write.csv(tab,file.path(out_dir,paste0('Lognormal_Reference_Survival_',tag,'.csv')),row.names=FALSE)
  tab
}
print(ref_table(fit_inf,'informative'))
print(ref_table(fit_vague,'vague'))

# STEP 11: Study B, empirically standardized survival curves
# For each posterior draw, average predictions over the SAME 521 patients.
# Then summarize distribution of that cohort-average across draws.
standardized <- function(fit,tag,n_draws=800L){
  ids <- get_draw_indices(fit,n_draws)
  grid <- seq(0,10,by=.1)
  groups <- list(
    stage=c('I','II','III','IV'),
    age=c('45','55','65','75'),
    grade=c('G1-G2','G3','G4'))
  rows <- list()
  for(kind in names(groups))for(value in groups[[kind]]){
    nd <- dat
    if(kind=='stage')nd <- make_stage(nd,value)
    if(kind=='age')nd <- make_age(nd,as.numeric(value))
    if(kind=='grade')nd <- make_grade(nd,value)
    mu <- brms::posterior_linpred(fit,newdata=nd,transform=FALSE,draw_ids=ids)
    sig <- as.numeric(posterior::as_draws_df(fit)$sigma)[ids]
    for(t in grid){
      avg <- if(t==0)rep(1,length(ids)) else
        rowMeans(pnorm(sweep(mu-log(t),1,sig,'/')))
      rows[[length(rows)+1]] <- data.frame(group=kind,level=value,Time_years=t,
        median=median(avg),lower=unname(quantile(avg,.025)),
        upper=unname(quantile(avg,.975)))
    }
  }
  res <- dplyr::bind_rows(rows)
  write.csv(res,file.path(out_dir,paste0('Lognormal_Standardized_Survival_',tag,'.csv')),
            row.names=FALSE)
  for(kind in names(groups)){
    sub <- subset(res,group==kind)
    p <- ggplot(sub,aes(Time_years,median,color=level,fill=level))+
      geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.12,colour=NA)+
      geom_line(linewidth=.85)+coord_cartesian(ylim=c(0,1))+
      labs(title=paste('Log-normal standardized survival by',kind),
           x='Years since diagnosis',y='Survival probability',color=kind,fill=kind)+
      theme_minimal()
    ggsave(file.path(out_dir,paste0('Lognormal_Standardized_',kind,'_',tag,'.png')),
           p,width=8,height=5,dpi=250)
  }
  res
}
std_inf <- standardized(fit_inf,'informative')
std_vague <- standardized(fit_vague,'vague')

# STEP 12: PSIS-LOO (patient-level censored log-likelihood, ordered by barcode)
# Important: use same patient order for ALL models' pointwise comparisons.
loo_ln <- loo::loo(fit_inf,save_psis=TRUE)
print(loo_ln)
saveRDS(loo_ln,file.path(out_dir,'LOO_Lognormal.rds'))
write.csv(data.frame(bcr_patient_barcode=dat$bcr_patient_barcode,
                     elpd_loo=loo_ln$pointwise[,'elpd_loo']),
          file.path(out_dir,'Lognormal_Pointwise_ELPD.csv'),row.names=FALSE)

# STEP 13: IPCW Brier and IBS (match Weibull evaluation script)
N <- nrow(dat)
cens_fit <- survival::survfit(survival::Surv(OS_years,1-event)~1,data=dat)
G_right <- function(x) vapply(x,function(z){
  idx <- which(cens_fit$time<=z)
  if(!length(idx)) 1 else cens_fit$surv[max(idx)]},numeric(1))
G_left <- function(x) vapply(x,function(z){
  idx <- which(cens_fit$time<z)
  if(!length(idx)) 1 else cens_fit$surv[max(idx)]},numeric(1))
G_Tminus <- G_left(dat$OS_years)
ipcw_brier <- function(t,S_hat){
  gt <- G_right(t)
  if(gt<=1e-8) warning('Censoring survival near zero at ',t)
  contribution <- numeric(N)
  died <- dat$OS_years<=t & dat$event==1
  alive_beyond <- dat$OS_years>t
  contribution[died] <- S_hat[died]^2/pmax(G_Tminus[died],1e-8)
  contribution[alive_beyond] <- (1-S_hat[alive_beyond])^2/pmax(gt,1e-8)
  mean(contribution)
}
ids_eval <- get_draw_indices(fit_inf,800)
mu_eval <- brms::posterior_linpred(fit_inf,newdata=dat,transform=FALSE,draw_ids=ids_eval)
sig_eval <- as.numeric(posterior::as_draws_df(fit_inf)$sigma)[ids_eval]
mean_surv <- function(t){
  if(t<=0)return(rep(1,N))
  colMeans(pnorm(sweep(mu_eval-log(t),1,sig_eval,'/')))
}
ibs_grid <- seq(0,10,by=.25)
bs <- vapply(ibs_grid,function(t)ipcw_brier(t,mean_surv(t)),numeric(1))
curve <- data.frame(Time_years=ibs_grid,Brier=bs)
write.csv(curve,file.path(out_dir,'Lognormal_Brier_Curve.csv'),row.names=FALSE)
IBS <- sum(diff(ibs_grid)*(head(bs,-1)+tail(bs,-1))/2)/10
brier_at <- function(t) bs[which.min(abs(ibs_grid-t))]

cat("\n===== LOG-NORMAL MODEL EVALUATION =====\n")

cat("Brier score at 1 year:", brier_at(1), "\n")
cat("Brier score at 5 years:", brier_at(5), "\n")
cat("Brier score at 10 years:", brier_at(10), "\n")
cat("Integrated Brier Score (0–10 years):", IBS, "\n")

# STEP 14: Summarize diagnostics and comparison metrics
sm <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(fit_inf)))
max_rhat <- max(sm$rhat,na.rm=TRUE)
min_bulk <- min(sm$ess_bulk,na.rm=TRUE)
min_tail <- min(sm$ess_tail,na.rm=TRUE)
np <- brms::nuts_params(fit_inf)
divergences <- sum(np$Parameter=='divergent__' & np$Value==1)
pk <- loo::pareto_k_values(loo_ln)
metrics <- data.frame(Model='Lognormal',N=N,Deaths=sum(dat$event),
  ELPD_LOO=unname(loo_ln$estimates['elpd_loo','Estimate']),
  SE_ELPD_LOO=unname(loo_ln$estimates['elpd_loo','SE']),
  p_LOO=unname(loo_ln$estimates['p_loo','Estimate']),
  LOOIC=unname(loo_ln$estimates['looic','Estimate']),
  SE_LOOIC=unname(loo_ln$estimates['looic','SE']),
  Max_Pareto_k=max(pk,na.rm=TRUE),Pareto_k_above_0_7=sum(pk>.7),
  Brier_1yr=brier_at(1),Brier_5yr=brier_at(5),Brier_10yr=brier_at(10),
  IBS_0_10yr=IBS,Max_Rhat=max_rhat,Min_ESS_bulk=min_bulk,
  Min_ESS_tail=min_tail,Divergences=divergences)
print(metrics)
write.csv(metrics,file.path(out_dir,'Lognormal_Model_Comparison_Metrics.csv'),row.names=FALSE)

# STEP 15: COMPARISON OF ALL THREE BAYESIAN MODELS
# ============================================================

library(loo)

other <- c(
  Weibull = "LOO_Weibull.rds",
  PiecewiseExp = "LOO_PiecewiseExp.rds"
)

other_paths <- file.path(dirname(input_file), unname(other))

# Load LOO objects
loo_w <- readRDS(other_paths[1])
loo_p <- readRDS(other_paths[2])

# Verify observation counts
print(c(
  Lognormal = nrow(loo_ln$pointwise),
  Weibull = nrow(loo_w$pointwise),
  PiecewiseExp = nrow(loo_p$pointwise)
))

# IMPORTANT:
# Confirm all three objects use the same 521 patients,
# in identical order, with comparable likelihood definitions.

# Perform model comparison
comparison <- loo::loo_compare(
  loo_ln,
  loo_w,
  loo_p
)

print(comparison)

# Save results
write.csv(
  as.data.frame(comparison),
  file.path(out_dir, "Three_Model_LOO_Comparison.csv")
)

cat("\nSTEP 15 COMPLETED SUCCESSFULLY\n")
