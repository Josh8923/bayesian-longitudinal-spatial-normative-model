#---------------------------------------- ========================
######### Extract Posterior Deviations and Model Comparison Metrics
####### Bayesian Longitudinal Spatial Normative Modeling
#=   ================----------------- ============================

library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(posterior)

source("R/utilities.R")

ensure_dir("oasis_application/outputs/model_objects")
ensure_dir("oasis_application/outputs/tables")

dat_stan <- readRDS("oasis_application/outputs/model_objects/dat_stan_input.rds")
stan_data <- readRDS("oasis_application/outputs/model_objects/stan_data_oasis.rds")
fit_A <- readRDS("oasis_application/outputs/model_objects/fit_A_independent_cross_sectional.rds")
fit_B <- readRDS("oasis_application/outputs/model_objects/fit_B_longitudinal_nonspatial.rds")
fit_C_stan <- readRDS("oasis_application/outputs/model_objects/fit_C_subject_specific_spatial_stan.rds")

fit_C_summary <- read_csv(
  "oasis_application/outputs/model_objects/fit_C_subject_specific_spatial_summary.csv",
  show_col_types = FALSE
)

key_summary <- read_csv(
  "oasis_application/outputs/model_objects/fit_C_subject_specific_spatial_key_summary.csv",
  show_col_types = FALSE
)

R <- stan_data$R
I <- stan_data$I
K_version <- stan_data$K_version

### EXTRACT POSTERIOR MEAN PARAMETERS

draws_C <- fit_C_stan$draws(format = "draws_matrix")

saveRDS(
  draws_C,
  "oasis_application/outputs/model_objects/fit_C_subject_specific_spatial_draws_matrix.rds"
)

alpha_mean <- get_mean_vector(draws_C, "alpha", R)
beta_age_mean <- get_mean_vector(draws_C, "beta_age", R)
gamma_version_mean <- get_mean_vector(draws_C, "gamma_version", K_version)
b_mean <- get_mean_vector(draws_C, "b", I)

u_cols <- grep("^u\\[", colnames(draws_C), value = TRUE)

u_index <- stringr::str_match(u_cols, "^u\\[(\\d+),(\\d+)\\]$")
u_index <- as.data.frame(u_index)
names(u_index) <- c("full", "region_id", "subject_id")

u_index$region_id <- as.integer(u_index$region_id)
u_index$subject_id <- as.integer(u_index$subject_id)
u_index$mean <- colMeans(draws_C[, u_cols, drop = FALSE])

u_mean_mat <- matrix(NA_real_, nrow = R, ncol = I)

for (k in seq_len(nrow(u_index))) {
  u_mean_mat[u_index$region_id[k], u_index$subject_id[k]] <- u_index$mean[k]
}

dat_stan$mu_C <-
  alpha_mean[dat_stan$region_id] +
  beta_age_mean[dat_stan$region_id] * dat_stan$age_z +
  gamma_version_mean[dat_stan$fs_version_id] +
  b_mean[dat_stan$subject_id] +
  u_mean_mat[cbind(dat_stan$region_id, dat_stan$subject_id)]

sigma_mean <- mean(draws_C[, "sigma"])

dat_stan$dev_C <- dat_stan$value_z - dat_stan$mu_C
dat_stan$z_C <- dat_stan$dev_C / sigma_mean

### BENCHMARK PREDICTIONS

dat_stan$mu_A <- predict(fit_A, newdata = dat_stan)
dat_stan$mu_B <- predict(fit_B, newdata = dat_stan, allow.new.levels = TRUE)

dat_stan <- dat_stan %>%
  mutate(
    dev_A = value_z - mu_A,
    dev_B = value_z - mu_B,
    z_A = dev_A / sd(dev_A, na.rm = TRUE),
    z_B = dev_B / sd(dev_B, na.rm = TRUE)
  )

write_csv(
  dat_stan,
  "oasis_application/outputs/model_objects/oasis3_final_model_deviation_scores.csv"
)

saveRDS(
  dat_stan,
  "oasis_application/outputs/model_objects/oasis3_final_model_deviation_scores.rds"
)

### MODEL COMPARISON METRICS

compare_tbl <- tibble(
  model = c(
    "Independent cross-sectional model",
    "Longitudinal non-spatial model",
    "Bayesian subject-specific spatial model"
  ),
  residual_sd = c(
    sd(dat_stan$dev_A, na.rm = TRUE),
    sd(dat_stan$dev_B, na.rm = TRUE),
    sd(dat_stan$dev_C, na.rm = TRUE)
  ),
  mean_abs_deviation = c(
    mean(abs(dat_stan$dev_A), na.rm = TRUE),
    mean(abs(dat_stan$dev_B), na.rm = TRUE),
    mean(abs(dat_stan$dev_C), na.rm = TRUE)
  ),
  rmse = c(
    sqrt(mean(dat_stan$dev_A^2, na.rm = TRUE)),
    sqrt(mean(dat_stan$dev_B^2, na.rm = TRUE)),
    sqrt(mean(dat_stan$dev_C^2, na.rm = TRUE))
  ),
  z_mean = c(
    mean(dat_stan$z_A, na.rm = TRUE),
    mean(dat_stan$z_B, na.rm = TRUE),
    mean(dat_stan$z_C, na.rm = TRUE)
  ),
  z_var = c(
    var(dat_stan$z_A, na.rm = TRUE),
    var(dat_stan$z_B, na.rm = TRUE),
    var(dat_stan$z_C, na.rm = TRUE)
  ),
  tail_probability = c(
    mean(abs(dat_stan$z_A) > 1.96, na.rm = TRUE),
    mean(abs(dat_stan$z_B) > 1.96, na.rm = TRUE),
    mean(abs(dat_stan$z_C) > 1.96, na.rm = TRUE)
  )
)

print(compare_tbl)

write_csv(
  compare_tbl,
  "oasis_application/outputs/tables/oasis3_model_comparison_metrics.csv"
)

### REGION AND SUBJECT DEVIATION SUMMARIES

region_dev_C <- dat_stan %>%
  group_by(region) %>%
  summarise(
    n_obs = n(),
    mean_z = mean(z_C, na.rm = TRUE),
    sd_z = sd(z_C, na.rm = TRUE),
    tail_prob = mean(abs(z_C) > 1.96, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_z)

write_csv(
  region_dev_C,
  "oasis_application/outputs/tables/oasis3_region_deviation_summary_model_C.csv"
)

subject_dev_C <- dat_stan %>%
  group_by(subject) %>%
  summarise(
    n_obs = n(),
    n_visits = n_distinct(mr_session),
    mean_abs_z = mean(abs(z_C), na.rm = TRUE),
    max_abs_z = max(abs(z_C), na.rm = TRUE),
    prop_extreme = mean(abs(z_C) > 1.96, na.rm = TRUE),
    mean_mmse = mean(mmse, na.rm = TRUE),
    mean_cdrtot = mean(cdrtot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_z))

write_csv(
  subject_dev_C,
  "oasis_application/outputs/tables/oasis3_subject_deviation_summary_model_C.csv"
)

analysis_objects <- list(
  dat_stan = dat_stan,
  compare_tbl = compare_tbl,
  region_dev_C = region_dev_C,
  subject_dev_C = subject_dev_C,
  fit_C_summary = fit_C_summary,
  key_summary = key_summary
)

saveRDS(
  analysis_objects,
  "oasis_application/outputs/model_objects/oasis_analysis_objects.rds"
)

