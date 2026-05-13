#---------------------------------------- ========================
######### Build OASIS-3 Modeling Data and Spatial Adjacency
####### Bayesian Longitudinal Spatial Normative Modeling
#=   ================----------------- ============================

library(dplyr)
library(stringr)
library(readr)

source("R/utilities.R")

ensure_dir("oasis_application/outputs/model_objects")

dat_long_2plus <- readRDS("oasis_application/outputs/oasis3_model_ready_long_2plus.rds")

### MODELING DATA

dat_model <- dat_long_2plus %>%
  filter(
    !is.na(value),
    !is.na(age_z),
    !is.na(region),
    !is.na(subject),
    !is.na(fs_version)
  ) %>%
  mutate(
    region = factor(region),
    subject = factor(subject),
    fs_version = factor(fs_version)
  )

dat_model <- dat_model %>%
  group_by(region) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup() %>%
  filter(!is.na(value_z))

cat("Subjects:", n_distinct(dat_model$subject), "\n")
cat("Regions:", n_distinct(dat_model$region), "\n")
cat("Observations:", nrow(dat_model), "\n")

region_levels <- levels(dat_model$region)

### BUILD SPATIAL ADJACENCY MATRIX

region_info <- tibble(
  region = region_levels,
  base_region = get_base_region(region_levels),
  hemisphere = case_when(
    str_detect(region, "^lh_|^left_") ~ "left",
    str_detect(region, "^rh_|^right_") ~ "right",
    TRUE ~ "midline"
  ),
  family = case_when(
    str_detect(base_region, "hippocampus|amygdala|entorhinal|parahippocampal|temporalpole") ~ "medial_temporal",
    str_detect(base_region, "inferiortemporal|middletemporal") ~ "lateral_temporal",
    str_detect(base_region, "precuneus|posteriorcingulate") ~ "posterior_medial",
    TRUE ~ "other"
  )
)

R <- length(region_levels)

W <- matrix(0, R, R)
rownames(W) <- colnames(W) <- region_levels

for (i in seq_len(R)) {
  for (j in seq_len(R)) {
    if (i != j) {
      same_base <- region_info$base_region[i] == region_info$base_region[j]
      same_family <- region_info$family[i] == region_info$family[j]

      if (same_base) {
        W[i, j] <- 1.00
      } else if (same_family) {
        W[i, j] <- 0.35
      }
    }
  }
}

diag(W) <- 0

D <- diag(rowSums(W))
rho <- 0.90

Q_stan <- D - rho * W
Q_stan <- Q_stan + diag(1e-4, R)

eigen_check <- eigen(Q_stan, symmetric = TRUE, only.values = TRUE)$values
cat("Minimum eigenvalue of Q_stan:", min(eigen_check), "\n")

if (min(eigen_check) <= 0) {
  stop("Q_stan is not positive definite. Increase ridge term.")
}

write.csv(W, "oasis_application/outputs/model_objects/oasis3_region_adjacency_matrix_W.csv")
write.csv(Q_stan, "oasis_application/outputs/model_objects/oasis3_region_precision_matrix_Q.csv")
write_csv(region_info, "oasis_application/outputs/model_objects/oasis3_region_info.csv")

### STAN DATA

dat_stan <- dat_model %>%
  mutate(
    subject_id = as.integer(factor(subject)),
    region_id = as.integer(factor(region, levels = region_levels)),
    fs_version_id = as.integer(factor(fs_version))
  ) %>%
  arrange(subject_id, region_id, age_z)

N <- nrow(dat_stan)
I <- n_distinct(dat_stan$subject_id)
R <- length(region_levels)
K_version <- n_distinct(dat_stan$fs_version_id)

stan_data <- list(
  N = N,
  I = I,
  R = R,
  K_version = K_version,
  subject_id = dat_stan$subject_id,
  region_id = dat_stan$region_id,
  fs_version_id = dat_stan$fs_version_id,
  age_z = dat_stan$age_z,
  y = dat_stan$value_z,
  Q = Q_stan
)

cat("Stan data dimensions:\n")
cat("N observations:", N, "\n")
cat("Subjects:", I, "\n")
cat("Regions:", R, "\n")
cat("FreeSurfer versions:", K_version, "\n")

saveRDS(dat_model, "oasis_application/outputs/model_objects/dat_model.rds")
saveRDS(dat_stan, "oasis_application/outputs/model_objects/dat_stan_input.rds")
saveRDS(stan_data, "oasis_application/outputs/model_objects/stan_data_oasis.rds")
saveRDS(region_levels, "oasis_application/outputs/model_objects/region_levels.rds")
saveRDS(Q_stan, "oasis_application/outputs/model_objects/Q_stan.rds")

