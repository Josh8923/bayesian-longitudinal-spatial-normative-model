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
sim_args <- scenarios$Scenario_5_MissingFollowup
stopifnot(!is.null(sim_args))

sim_results <- run_simulation_study(
  M = 1000,
  scenario_name = "Scenario_5_MissingFollowup",
  sim_args = sim_args,
  stan_file = "blsnm_simulation_model.stan",
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 2029,
  adapt_delta = 0.99,
  max_treedepth = 13
)

saveRDS(sim_results, "sim_results_Scenario_5_MissingFollowup.rds")
