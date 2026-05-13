#---------------------------------------- ========================
######### Fit OASIS-3 Benchmark and Bayesian Spatial Models
####### Bayesian Longitudinal Spatial Normative Modeling
#=   ================----------------- ============================

library(dplyr)
library(lme4)
library(cmdstanr)
library(posterior)
library(readr)

source("R/utilities.R")

ensure_dir("oasis_application/outputs/model_objects")
ensure_dir("oasis_application/outputs/model_objects/cmdstan_output_clean")

dat_model <- readRDS("oasis_application/outputs/model_objects/dat_model.rds")
stan_data <- readRDS("oasis_application/outputs/model_objects/stan_data_oasis.rds")

### MODEL A: INDEPENDENT CROSS-SECTIONAL MODEL

fit_A <- lm(
  value_z ~ 0 + region + region:age_z + fs_version,
  data = dat_model
)

saveRDS(
  fit_A,
  "oasis_application/outputs/model_objects/fit_A_independent_cross_sectional.rds"
)

### MODEL B: LONGITUDINAL NON-SPATIAL MIXED MODEL

fit_B <- lmer(
  value_z ~ 0 + region + region:age_z + fs_version + (1 | subject),
  data = dat_model,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

saveRDS(
  fit_B,
  "oasis_application/outputs/model_objects/fit_B_longitudinal_nonspatial.rds"
)

### MODEL C: BAYESIAN SUBJECT-SPECIFIC SPATIAL MODEL

mod_C <- cmdstan_model(
  "stan/blsnm_oasis_model.stan",
  force_recompile = TRUE
)

fit_C_stan <- mod_C$sample(
  data = stan_data,
  seed = 89,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 100,
  adapt_delta = 0.99,
  max_treedepth = 13,
  output_dir = "oasis_application/outputs/
  model_objects/cmdstan_output_clean"
)

saveRDS(
  fit_C_stan,
  "oasis_application/outputs/model_objects/
  fit_C_subject_specific_spatial_stan.rds"
)

fit_C_summary <- fit_C_stan$summary()

write_csv(
  fit_C_summary,
  "oasis_application/outputs/model_objects/
  fit_C_subject_specific_spatial_summary.csv"
)

key_summary <- fit_C_stan$summary(
  c("sigma", "sigma_b", "tau_u", "sigma_version")
)

print(key_summary)

write_csv(
  key_summary,
  "oasis_application/outputs/model_objects/
  fit_C_subject_specific_spatial_key_summary.csv"
)

#"Worst Rhat:",
max(fit_C_summary$rhat, na.rm = TRUE)
# "Minimum bulk ESS:"
min(fit_C_summary$ess_bulk, na.rm = TRUE)
## "Minimum tail ESS:"
fit_C_summary$ess_tail, na.rm = TRUE)

