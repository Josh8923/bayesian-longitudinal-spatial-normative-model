library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(Matrix)
library(mvtnorm)
library(lme4)
library(cmdstanr)
library(posterior)

source("simulation_function.R")
scenarios <- get_simulation_scenarios()
sim_args <- scenarios$Scenario_1_NoSpatial
stopifnot(!is.null(sim_args))

sim_results <- run_simulation_study(
  M = 1000,
  scenario_name = "Scenario_1_NoSpatial",
  sim_args = sim_args,
  stan_file = "blsnm_simulation_model.stan",
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 2025,
  adapt_delta = 0.95,
  max_treedepth = 12
)

saveRDS(sim_results, "sim_results_Scenario_1_NoSpatial.rds")
