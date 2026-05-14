# ##------------------------------------------------------------
### Bayesian longitudinal spatial normative model simulation study
####. Core utilities for data generation, model fitting, and evaluation.
# ##------------------------------------------------------------

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



`%||%` <- function(x, y) if (is.null(x)) y else x  ## # function to clean results later


## Spatial adjacency and precision matrix utilities
make_chain_adjacency <- function(R) {
  if (!is.numeric(R) || length(R) != 1 || R < 2) {
    stop("R must be a single integer >= 2.", call. = FALSE)
  }

  W <- matrix(0, nrow = R, ncol = R)
  for (r in seq_len(R - 1)) {
    W[r, r + 1] <- 1
    W[r + 1, r] <- 1
  }

  W
}

make_precision_matrix <- function(W, rho, jitter = 1e-5) {
  if (!is.matrix(W) || nrow(W) != ncol(W)) {
    stop("W must be a square matrix.", call. = FALSE)
  }

  if (rho < 0 || rho >= 1) {
    stop("rho must lie in [0, 1).", call. = FALSE)
  }

  D <- diag(rowSums(W))
  Q <- D - rho * W + diag(jitter, nrow(W))  ### Proper CAR precision
  eigvals <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values

  if (min(eigvals) <= 0) {
    stop(
      "Precision matrix is not positive definite. Reduce rho or increase jitter.",
      call. = FALSE
    )
  }

  Q
}

## Design matrix utilities
make_design_matrix <- function(age, sex, nonlinear = FALSE) {
  if (!nonlinear) {
    X <- cbind(1, age, sex)
    colnames(X) <- c("Intercept", "Age", "Sex")
  } else {
    X <- cbind(1, age, age^2, sex)
    colnames(X) <- c("Intercept", "Age", "Age2", "Sex")
  }

  X
}

standardize_design_matrix <- function(X) {
  X_std <- X
  means <- rep(0, ncol(X))
  sds <- rep(1, ncol(X))

  if (ncol(X) >= 2) {
    for (j in 2:ncol(X)) {
      means[j] <- mean(X[, j], na.rm = TRUE)
      sds[j] <- sd(X[, j], na.rm = TRUE)
      if (is.na(sds[j]) || sds[j] == 0) sds[j] <- 1
      X_std[, j] <- (X[, j] - means[j]) / sds[j]
    }
  }

  list(X = X_std, means = means, sds = sds)
}


# #Abnormality signal generation
make_abnormal_shift <- function(R,
                                pattern = c("localized", "diffuse"),
                                effect_size = 0.8,
                                cluster_size = max(3, floor(0.15 * R)),
                                diffuse_prop = 0.30) {
  pattern <- match.arg(pattern)
  delta <- rep(0, R)

  if (pattern == "localized") {
    start_idx <- sample.int(max(1, R - cluster_size + 1), 1)
    abnormal_regions <- start_idx:(start_idx + cluster_size - 1)
    delta[abnormal_regions] <- effect_size
  } else {
    abnormal_regions <- sort(
      sample.int(R, size = max(1, floor(diffuse_prop * R)), replace = FALSE)
    )
    delta[abnormal_regions] <- effect_size
  }

  list(delta = delta, abnormal_regions = abnormal_regions)
}


## Longitudinal neuroimaging data generation
simulate_brain_longitudinal <- function(
    n_subj = 120,
    R = 20,
    visits = 3,
    visits_range = c(2L, 5L),
    variable_visits = FALSE,
    nonlinear = FALSE,
    rho = 0.5,
    sigma_b = 0.5,
    tau_u = 0.8,
    sigma = 0.6,
    beta_intercept_mean = 5.0,
    beta_intercept_sd = 0.5,
    beta_age_mean = -0.03,
    beta_age_sd = 0.01,
    beta_sex_mean = 0.15,
    beta_sex_sd = 0.05,
    beta_age2_mean = 0.00015,
    beta_age2_sd = 0.00005,
    proportion_abnormal = 0.25,
    abnormal_pattern = c("localized", "diffuse"),
    abnormal_effect_size = 0.8,
    missing_followup = FALSE,
    missing_strength = 0.15,
    W = NULL,
    seed = NULL
) {
  abnormal_pattern <- match.arg(abnormal_pattern)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (is.null(W)) {
    W <- make_chain_adjacency(R)
  }

  Q <- make_precision_matrix(W, rho = rho)
  Sigma_u <- solve(Q + diag(1e-5, nrow(Q))) * tau_u^2  # Spatial covariance

  beta0 <- rnorm(R, mean = beta_intercept_mean, sd = beta_intercept_sd)
  beta1 <- rnorm(R, mean = beta_age_mean, sd = beta_age_sd)
  beta2 <- rnorm(R, mean = beta_sex_mean, sd = beta_sex_sd)

  if (nonlinear) {
    beta3 <- rnorm(R, mean = beta_age2_mean, sd = beta_age2_sd)
    beta_mat <- rbind(beta0, beta1, beta3, beta2)
    rownames(beta_mat) <- c("Intercept", "Age", "Age2", "Sex")
  } else {
    beta_mat <- rbind(beta0, beta1, beta2)
    rownames(beta_mat) <- c("Intercept", "Age", "Sex")
  }

  n_abnormal <- floor(proportion_abnormal * n_subj)
  abnormal_ids <- if (n_abnormal > 0) {
    sort(sample.int(n_subj, size = n_abnormal, replace = FALSE))
  } else {
    integer(0)
  }

  abnormal_object <- make_abnormal_shift(
    R = R,
    pattern = abnormal_pattern,
    effect_size = abnormal_effect_size
  )

  delta_vec <- abnormal_object$delta
  long_data <- vector("list", n_subj)
  truth_u <- matrix(NA_real_, nrow = n_subj, ncol = R)
  truth_b <- rep(NA_real_, n_subj)

  for (i in seq_len(n_subj)) {
    Ti <- if (variable_visits) {
      sample(seq.int(visits_range[1], visits_range[2]), size = 1)
    } else {
      visits
    }

    age0 <- runif(1, min = 60, max = 85)
    sex_i <- rbinom(1, size = 1, prob = 0.5)
    b_i <- rnorm(1, mean = 0, sd = sigma_b)
    u_i <- as.numeric(
      mvtnorm::rmvnorm(1, mean = rep(0, R), sigma = Sigma_u)
    )

    truth_b[i] <- b_i
    truth_u[i, ] <- u_i

    if (Ti == 1) {
      age_vec <- age0
    } else {
      increments <- c(0, cumsum(runif(Ti - 1, min = 0.8, max = 2.0)))
      age_vec <- age0 + increments
    }

    X_i <- make_design_matrix(age = age_vec, sex = sex_i, nonlinear = nonlinear)
    mean_mat <- X_i %*% beta_mat  # ##Region-specific normative mean

    subject_df <- expand.grid(
      subject = i,
      visit = seq_len(Ti),
      region = seq_len(R),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      mutate(
        age = age_vec[visit],
        sex = sex_i,
        abnormal_subject = as.integer(i %in% abnormal_ids)
      )

    y_vec <- numeric(nrow(subject_df))
    mu_vec <- numeric(nrow(subject_df))

    row_counter <- 1L
    for (t in seq_len(Ti)) {
      for (r in seq_len(R)) {
        mu_itr <- mean_mat[t, r] + b_i + u_i[r] # # Mean + subject + spatial deviation
        if (i %in% abnormal_ids) {
          mu_itr <- mu_itr + delta_vec[r]
        }

        y_itr <- rnorm(1, mean = mu_itr, sd = sigma)

        mu_vec[row_counter] <- mu_itr
        y_vec[row_counter] <- y_itr
        row_counter <- row_counter + 1L
      }
    }

    subject_df <- subject_df %>%
      mutate(
        y = y_vec,
        mu_true = mu_vec,
        abnormal_region = as.integer(region %in% abnormal_object$abnormal_regions),
        u_true = u_i[region],
        b_true = b_i,
        delta_true = ifelse(i %in% abnormal_ids, delta_vec[region], 0)
      )

    if (missing_followup && Ti > 1) {
      keep_flag <- rep(TRUE, nrow(subject_df))

      visit_age_df <- subject_df %>%
        distinct(visit, age) %>%
        arrange(visit)

      age_mean <- mean(visit_age_df$age, na.rm = TRUE)
      age_sd <- sd(visit_age_df$age, na.rm = TRUE)
      if (is.na(age_sd) || age_sd == 0) age_sd <- 1

      later_visits <- visit_age_df$visit[visit_age_df$visit > 1]

      for (vv in later_visits) {
        idx <- which(subject_df$visit == vv)

        age_vv <- visit_age_df$age[visit_age_df$visit == vv][1]
        age_z <- (age_vv - age_mean) / age_sd

        prior_idx <- which(subject_df$visit < vv)
        prior_burden <- if (length(prior_idx) == 0) {
          0
        } else {
          mean(abs(subject_df$u_true[prior_idx]), na.rm = TRUE)
        }

        if (is.na(prior_burden)) prior_burden <- 0

        logit_p_drop <- -2 + missing_strength * age_z + 0.8 * prior_burden
        p_drop <- plogis(logit_p_drop)

        if (is.na(p_drop) || !is.finite(p_drop)) p_drop <- 0
        p_drop <- min(max(p_drop, 0), 1)

        if (runif(1) < p_drop) {
          keep_flag[idx] <- FALSE
        }
      }

      subject_df <- subject_df[keep_flag, , drop = FALSE]
    }

    long_data[[i]] <- subject_df
  }

  long_data <- bind_rows(long_data) %>%
    arrange(subject, visit, region) %>%
    mutate(
      subject = as.integer(subject),
      visit = as.integer(visit),
      region = as.integer(region)
    )

  truth <- list(
    beta_mat = beta_mat,
    b_true = truth_b,
    u_true = truth_u,
    abnormal_subjects = abnormal_ids,
    abnormal_regions = abnormal_object$abnormal_regions,
    delta = delta_vec,
    W = W,
    Q = Q,
    Sigma_u = Sigma_u,
    params = list(
      nonlinear = nonlinear,
      rho = rho,
      sigma_b = sigma_b,
      tau_u = tau_u,
      sigma = sigma
    )
  )

  list(data = long_data, truth = truth)
}


# ##Benchmark model A: independent cross-sectional model
fit_model_A <- function(dat, nonlinear = FALSE) {
  regions <- sort(unique(dat$region))
  fits <- vector("list", length(regions))
  names(fits) <- as.character(regions)

  for (r in regions) {
    d_r <- dat %>% filter(region == r)

    if (!nonlinear) {
      form <- as.formula("y ~ age + sex")
    } else {
      form <- as.formula("y ~ age + I(age^2) + sex")
    }

    fit <- lm(form, data = d_r)

    sigma_r <- summary(fit)$sigma
    if (!is.finite(sigma_r) || sigma_r <= 0) sigma_r <- sd(d_r$y)
    if (!is.finite(sigma_r) || sigma_r <= 0) sigma_r <- 1

    pred <- predict(fit, newdata = d_r)
    residuals_r <- d_r$y - pred

    u_hat_sr <- d_r %>%
      mutate(resid = residuals_r) %>%
      group_by(subject, region) %>%
      summarise(u_hat = mean(resid), .groups = "drop")

    fits[[as.character(r)]] <- list(
      fit = fit,
      sigma = sigma_r,
      pred_df = d_r %>% mutate(mu_hat = pred, z_hat = (y - pred) / sigma_r),
      u_hat_sr = u_hat_sr
    )
  }

  pred_df <- bind_rows(lapply(fits, function(x) x$pred_df)) %>%
    arrange(subject, visit, region)

  u_hat_df <- bind_rows(lapply(fits, function(x) x$u_hat_sr)) %>%
    arrange(subject, region)

  list(
    model = fits,
    pred_df = pred_df,
    u_hat_df = u_hat_df,
    name = "Model_A_Independent_CS"
  )
}


# ##Benchmark model B: longitudinal non-spatial model
fit_model_B <- function(dat, nonlinear = FALSE) {
  regions <- sort(unique(dat$region))
  fits <- vector("list", length(regions))
  names(fits) <- as.character(regions)

  ctrl <- lmerControl(
    check.rankX = "silent.drop.cols",
    check.conv.singular = "ignore",
    check.conv.grad = "ignore",
    check.conv.hess = "ignore"
  )

  for (r in regions) {
    d_r <- dat %>% filter(region == r)

    if (!nonlinear) {
      form <- as.formula("y ~ age + sex + (1 | subject)")
    } else {
      form <- as.formula("y ~ age + I(age^2) + sex + (1 | subject)")
    }

    fit <- suppressWarnings(
      lmer(form, data = d_r, REML = TRUE, control = ctrl)
    )

    pred <- predict(fit, newdata = d_r, allow.new.levels = TRUE)
    sigma_r <- sigma(fit)
    if (!is.finite(sigma_r) || sigma_r <= 0) sigma_r <- sd(d_r$y)
    if (!is.finite(sigma_r) || sigma_r <= 0) sigma_r <- 1

    ranef_df <- ranef(fit)$subject %>%
      tibble::rownames_to_column("subject") %>%
      rename(b_hat = `(Intercept)`) %>%
      mutate(subject = as.integer(subject))

    u_hat_sr <- d_r %>%
      left_join(ranef_df, by = "subject") %>%
      mutate(resid = y - pred) %>%
      group_by(subject, region) %>%
      summarise(u_hat = mean(resid), .groups = "drop")

    fits[[as.character(r)]] <- list(
      fit = fit,
      sigma = sigma_r,
      pred_df = d_r %>% mutate(mu_hat = pred, z_hat = (y - pred) / sigma_r),
      u_hat_sr = u_hat_sr
    )
  }

  pred_df <- bind_rows(lapply(fits, function(x) x$pred_df)) %>%
    arrange(subject, visit, region)

  u_hat_df <- bind_rows(lapply(fits, function(x) x$u_hat_sr)) %>%
    arrange(subject, region)

  list(
    model = fits,
    pred_df = pred_df,
    u_hat_df = u_hat_df,
    name = "Model_B_Longitudinal_NoSpatial"
  )
}



# Proposed Bayesian longitudinal spatial model in Stan
brain_spatial_stan_code <- '
data {
  int<lower=1> N;
  int<lower=1> N_subj;
  int<lower=1> R;
  int<lower=1> p;
  matrix[N, p] X;
  vector[N] y_std;
  array[N] int<lower=1, upper=N_subj> subj;
  array[N] int<lower=1, upper=R> region;
  matrix[R, R] U_base;
}
parameters {
  matrix[p, R] beta;
  vector[N_subj] b_raw;
  matrix[N_subj, R] u_raw;
  real<lower=1e-4> sigma;
  real<lower=1e-4> sigma_b;
  real<lower=1e-4> tau_u;
}
transformed parameters {
  vector[N_subj] b;
  matrix[N_subj, R] u;

  b = sigma_b * b_raw;
  u = (u_raw * U_base) * tau_u;
}
model {
  to_vector(beta) ~ normal(0, 1);
  b_raw ~ normal(0, 1);
  to_vector(u_raw) ~ normal(0, 1);

  sigma ~ normal(0, 0.5);
  sigma_b ~ normal(0, 0.5);
  tau_u ~ normal(0, 0.5);

  for (n in 1:N) {
    real mu_n;
    mu_n = dot_product(X[n], beta[, region[n]]) + b[subj[n]] + u[subj[n], region[n]];
    y_std[n] ~ normal(mu_n, sigma);
  }
}
generated quantities {
  vector[N] mu_std;
  vector[N] z_hat;

  for (n in 1:N) {
    mu_std[n] = dot_product(X[n], beta[, region[n]]) + b[subj[n]] + u[subj[n], region[n]];
    z_hat[n] = (y_std[n] - mu_std[n]) / sigma;
  }
}
'

write_brain_spatial_stan <- function(stan_file = "brain_spatial_normative_model.stan") {
  writeLines(brain_spatial_stan_code, con = stan_file)
  invisible(stan_file)
}



## # Stan data preparation and posterior extraction
prepare_stan_data <- function(dat, Q, nonlinear = FALSE) {
  X_raw <- make_design_matrix(age = dat$age, sex = dat$sex, nonlinear = nonlinear)
  X_scaled <- standardize_design_matrix(X_raw)

  y_mean <- mean(dat$y, na.rm = TRUE)
  y_sd <- sd(dat$y, na.rm = TRUE)
  if (is.na(y_sd) || y_sd == 0) y_sd <- 1

  Sigma_base <- solve(Q + diag(1e-5, nrow(Q)))
  U_base <- chol(Sigma_base)

  list(
    stan_data = list(
      N = nrow(dat),
      N_subj = length(unique(dat$subject)),
      R = length(unique(dat$region)),
      p = ncol(X_scaled$X),
      X = unname(X_scaled$X),
      y_std = as.numeric((dat$y - y_mean) / y_sd),
      subj = dat$subject,
      region = dat$region,
      U_base = U_base
    ),
    transform = list(
      y_mean = y_mean,
      y_sd = y_sd,
      x_means = X_scaled$means,
      x_sds = X_scaled$sds
    )
  )
}

extract_posterior_means <- function(fit, variable) {
  dr <- fit$draws(variable = variable, format = "draws_matrix")
  colMeans(dr)
}

fit_model_C <- function(dat,
                        Q,
                        nonlinear = FALSE,
                        stan_file = "brain_spatial_normative_model.stan",
                        mod = NULL,
                        chains = 2,
                        parallel_chains = 2,
                        iter_warmup = 500,
                        iter_sampling = 500,
                        refresh = 100,
                        seed = 89,
                        adapt_delta = 0.95,
                        max_treedepth = 12) {
  if (is.null(mod)) {
    if (!file.exists(stan_file)) write_brain_spatial_stan(stan_file)
    mod <- cmdstan_model(stan_file)
  }

  prep <- prepare_stan_data(dat, Q = Q, nonlinear = nonlinear)

  fit <- tryCatch(
    mod$sample(
      data = prep$stan_data,
      seed = seed,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = refresh,
      show_messages = FALSE,
      show_exceptions = TRUE,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    pred_df <- dat %>%
      mutate(mu_hat = NA_real_, z_hat = NA_real_)

    u_hat_df <- expand.grid(
      subject = sort(unique(dat$subject)),
      region = sort(unique(dat$region)),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      mutate(u_hat = NA_real_) %>%
      arrange(subject, region)

    return(list(
      model = NULL,
      pred_df = pred_df,
      u_hat_df = u_hat_df,
      name = "Model_C_Bayesian_Longitudinal_Spatial",
      diagnostics = list(failed = TRUE)
    ))
  }

  mu_std_hat <- extract_posterior_means(fit, "mu_std")
  z_hat <- extract_posterior_means(fit, "z_hat")
  u_mean <- extract_posterior_means(fit, "u")

  mu_hat <- prep$transform$y_mean + prep$transform$y_sd * mu_std_hat

  pred_df <- dat %>%
    mutate(
      mu_hat = mu_hat,
      z_hat = z_hat
    )

  u_index <- tibble(
    param = names(u_mean),
    u_hat = as.numeric(u_mean)
  ) %>%
    mutate(
      inside = stringr::str_match(param, "u\\[(\\d+),(\\d+)\\]"),
      subject = as.integer(inside[, 2]),
      region = as.integer(inside[, 3])
    ) %>%
    dplyr::select(subject, region, u_hat) %>%
    arrange(subject, region)

  diag_summary <- fit$summary()  # Convergence and ESS diagnostics

  rhat_vals <- diag_summary$rhat
  ess_bulk_vals <- diag_summary$ess_bulk
  ess_tail_vals <- diag_summary$ess_tail

  rhat_vals <- rhat_vals[is.finite(rhat_vals)]
  ess_bulk_vals <- ess_bulk_vals[is.finite(ess_bulk_vals)]
  ess_tail_vals <- ess_tail_vals[is.finite(ess_tail_vals)]

  sampler_diagnostics <- list(
    failed = FALSE,
    max_rhat = if (length(rhat_vals) == 0) NA_real_ else max(rhat_vals),
    min_ess_bulk = if (length(ess_bulk_vals) == 0) NA_real_ else min(ess_bulk_vals),
    min_ess_tail = if (length(ess_tail_vals) == 0) NA_real_ else min(ess_tail_vals)
  )

  list(
    model = fit,
    pred_df = pred_df,
    u_hat_df = u_index,
    name = "Model_C_Bayesian_Longitudinal_Spatial",
    diagnostics = sampler_diagnostics
  )
}


# #Simulation evaluation metrics
compute_normative_metrics <- function(pred_df, reference_only = TRUE) {
  eval_df <- pred_df

  if (reference_only) {
    eval_df <- eval_df %>% filter(abnormal_subject == 0)
  }

  list(
    mean_bias = mean(eval_df$mu_hat - eval_df$mu_true, na.rm = TRUE),
    mean_mse = mean((eval_df$mu_hat - eval_df$mu_true)^2, na.rm = TRUE),
    z_mean = mean(eval_df$z_hat, na.rm = TRUE),
    z_var = var(eval_df$z_hat, na.rm = TRUE),
    z_tail_prop = mean(abs(eval_df$z_hat) > 1.96, na.rm = TRUE)
  )
}

compute_map_mse <- function(u_hat_df, truth_u_df) {
  merged <- truth_u_df %>%
    left_join(u_hat_df, by = c("subject", "region"))

  if (all(is.na(merged$u_hat))) {
    return(NA_real_)
  }

  mean((merged$u_hat - merged$u_true)^2, na.rm = TRUE)
}

compute_detection_metrics <- function(pred_df, threshold = 1.96) {
  abnormal_eval <- pred_df %>%
    mutate(flag_hat = as.integer(abs(z_hat) > threshold))

  tp <- sum(abnormal_eval$flag_hat == 1 &
              abnormal_eval$abnormal_region == 1 &
              abnormal_eval$abnormal_subject == 1,
            na.rm = TRUE)

  fp <- sum(abnormal_eval$flag_hat == 1 &
              !(abnormal_eval$abnormal_region == 1 &
                  abnormal_eval$abnormal_subject == 1),
            na.rm = TRUE)

  fn <- sum(abnormal_eval$flag_hat == 0 &
              abnormal_eval$abnormal_region == 1 &
              abnormal_eval$abnormal_subject == 1,
            na.rm = TRUE)

  tn <- sum(abnormal_eval$flag_hat == 0 &
              !(abnormal_eval$abnormal_region == 1 &
                  abnormal_eval$abnormal_subject == 1),
            na.rm = TRUE)

  list(
    sensitivity = if ((tp + fn) == 0) NA_real_ else tp / (tp + fn),
    specificity = if ((tn + fp) == 0) NA_real_ else tn / (tn + fp),
    ppv = if ((tp + fp) == 0) NA_real_ else tp / (tp + fp),
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn
  )
}

assemble_truth_u_df <- function(sim_obj) {
  expand.grid(
    subject = seq_len(nrow(sim_obj$truth$u_true)),
    region = seq_len(ncol(sim_obj$truth$u_true)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(u_true = as.vector(t(sim_obj$truth$u_true))) %>%
    arrange(subject, region)
}

collect_model_metrics <- function(fit_obj, truth_u_df) {
  metrics <- compute_normative_metrics(fit_obj$pred_df)
  map_mse <- compute_map_mse(fit_obj$u_hat_df, truth_u_df)
  det_metrics <- compute_detection_metrics(fit_obj$pred_df)

  tibble(
    model = fit_obj$name,
    mean_bias = metrics$mean_bias,
    mean_mse = metrics$mean_mse,
    z_mean = metrics$z_mean,
    z_var = metrics$z_var,
    z_tail_prop = metrics$z_tail_prop,
    map_mse = map_mse,
    sensitivity = det_metrics$sensitivity,
    specificity = det_metrics$specificity,
    ppv = det_metrics$ppv,
    stan_failed = isTRUE(fit_obj$diagnostics$failed %||% FALSE),
    stan_max_rhat = fit_obj$diagnostics$max_rhat %||% NA_real_,
    stan_min_ess_bulk = fit_obj$diagnostics$min_ess_bulk %||% NA_real_,
    stan_min_ess_tail = fit_obj$diagnostics$min_ess_tail %||% NA_real_
  )
}


# ##Monte Carlo replicate and simulation-study wrappers
run_one_replicate <- function(sim_args,
                              stan_file = "brain_spatial_normative_model.stan",
                              mod = NULL,
                              chains = 2,
                              parallel_chains = 2,
                              iter_warmup = 500,
                              iter_sampling = 500,
                              refresh = 100,
                              seed = 1234,
                              adapt_delta = 0.95,
                              max_treedepth = 12) {
  sim_obj <- do.call(simulate_brain_longitudinal, sim_args)
  dat <- sim_obj$data
  truth_u_df <- assemble_truth_u_df(sim_obj)
  nonlinear <- isTRUE(sim_args$nonlinear %||% FALSE)

  fit_A <- fit_model_A(dat, nonlinear = nonlinear)
  fit_B <- fit_model_B(dat, nonlinear = nonlinear)
  fit_C <- fit_model_C(
    dat = dat,
    Q = sim_obj$truth$Q,
    nonlinear = nonlinear,
    stan_file = stan_file,
    mod = mod,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = refresh,
    seed = seed,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )

  summary_tbl <- bind_rows(
    collect_model_metrics(fit_A, truth_u_df),
    collect_model_metrics(fit_B, truth_u_df),
    collect_model_metrics(fit_C, truth_u_df)
  )

  list(
    sim_obj = sim_obj,
    fits = list(fit_A = fit_A, fit_B = fit_B, fit_C = fit_C),
    summary = summary_tbl
  )
}

run_simulation_study <- function(
    M = 10,
    scenario_name = "Scenario_1_NoSpatial",
    sim_args = list(),
    stan_file = "brain_spatial_normative_model.stan",
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 500,
    iter_sampling = 500,
    refresh = 100,
    seed = 2025,
    adapt_delta = 0.95,
    max_treedepth = 12
) {
  if (is.null(sim_args)) {
    stop("sim_args is NULL. Use a valid scenario name.", call. = FALSE)
  }

  if (!file.exists(stan_file)) {
    write_brain_spatial_stan(stan_file)
  }

  mod <- cmdstan_model(stan_file)
  out <- vector("list", M)

  for (m in seq_len(M)) {
    cat("Running", scenario_name, "replicate", m, "of", M, "\n")

    sim_args_m <- sim_args
    sim_args_m$seed <- seed + m

    out[[m]] <- run_one_replicate(
      sim_args = sim_args_m,
      stan_file = stan_file,
      mod = mod,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      refresh = refresh,
      seed = seed + m,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    )$summary %>%
      mutate(
        replicate = m,
        scenario = scenario_name
      )
  }

  results <- bind_rows(out)  ### Replicate-level results

  summary_table <- results %>%
    group_by(scenario, model) %>%
    summarise(
      across(
        c(
          mean_bias, mean_mse, z_mean, z_var, z_tail_prop,
          map_mse, sensitivity, specificity, ppv
        ),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          mcse = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
        ),
        .names = "{.col}_{.fn}"
      ),
      stan_failure_rate = mean(stan_failed, na.rm = TRUE),
      stan_max_rhat_mean = mean(stan_max_rhat, na.rm = TRUE),
      stan_min_ess_bulk_mean = mean(stan_min_ess_bulk, na.rm = TRUE),
      stan_min_ess_tail_mean = mean(stan_min_ess_tail, na.rm = TRUE),
      .groups = "drop"
    )

  list(raw = results, summary = summary_table)
}


# Simulation scenario definitions
get_simulation_scenarios <- function() {
  scenarios <- list(
    Scenario_1_NoSpatial = list(
      n_subj = 120,
      R = 20,
      visits = 3,
      variable_visits = FALSE,
      nonlinear = FALSE,
      rho = 0.0,
      sigma_b = 0.5,
      tau_u = 0.8,
      sigma = 0.6,
      proportion_abnormal = 0.25,
      abnormal_pattern = "localized",
      abnormal_effect_size = 0.8,
      missing_followup = FALSE
    ),
    Scenario_2_ModerateSpatial = list(
      n_subj = 120,
      R = 20,
      visits = 3,
      variable_visits = FALSE,
      nonlinear = FALSE,
      rho = 0.4,
      sigma_b = 0.5,
      tau_u = 0.8,
      sigma = 0.6,
      proportion_abnormal = 0.25,
      abnormal_pattern = "localized",
      abnormal_effect_size = 0.8,
      missing_followup = FALSE
    ),
    Scenario_3_StrongSpatial = list(
      n_subj = 120,
      R = 20,
      visits = 3,
      variable_visits = FALSE,
      nonlinear = FALSE,
      rho = 0.8,
      sigma_b = 0.5,
      tau_u = 0.8,
      sigma = 0.6,
      proportion_abnormal = 0.25,
      abnormal_pattern = "localized",
      abnormal_effect_size = 0.8,
      missing_followup = FALSE
    ),
    Scenario_4_VariableVisits = list(
      n_subj = 120,
      R = 20,
      variable_visits = TRUE,
      visits_range = c(2L, 5L),
      nonlinear = FALSE,
      rho = 0.5,
      sigma_b = 0.5,
      tau_u = 0.8,
      sigma = 0.6,
      proportion_abnormal = 0.25,
      abnormal_pattern = "localized",
      abnormal_effect_size = 0.8,
      missing_followup = FALSE
    ),
    Scenario_5_MissingFollowup = list(
      n_subj = 120,
      R = 20,
      variable_visits = TRUE,
      visits_range = c(3L, 5L),
      nonlinear = FALSE,
      rho = 0.5,
      sigma_b = 0.5,
      tau_u = 0.8,
      sigma = 0.6,
      proportion_abnormal = 0.25,
      abnormal_pattern = "localized",
      abnormal_effect_size = 0.8,
      missing_followup = TRUE,
      missing_strength = 0.25
    ),
    Scenario_6_NonlinearMean = list(
      n_subj = 120,
      R = 20,
      visits = 3,
      variable_visits = FALSE,
      nonlinear = TRUE,
      rho = 0.5,
      sigma_b = 0.5,
      tau_u = 0.8,
      sigma = 0.6,
      proportion_abnormal = 0.25,
      abnormal_pattern = "diffuse",
      abnormal_effect_size = 0.8,
      missing_followup = FALSE
    )
  )

  scenarios$Scenario_4_UnequalVisits <- scenarios$Scenario_4_VariableVisits
  scenarios$Scenario_6_Nonlinear <- scenarios$Scenario_6_NonlinearMean

  scenarios
}


