data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=1> R;
  int<lower=1> K_version;

  array[N] int<lower=1, upper=I> subject_id;
  array[N] int<lower=1, upper=R> region_id;
  array[N] int<lower=1, upper=K_version> fs_version_id;

  vector[N] age_z;
  vector[N] y;

  matrix[R, R] Q;
}

transformed data {
  matrix[R, R] Q_jitter;
  matrix[R, R] L_cov;

  Q_jitter = Q + diag_matrix(rep_vector(1e-3, R));
  L_cov = cholesky_decompose(inverse_spd(Q_jitter));
}

parameters {
  vector[R] alpha;
  vector[R] beta_age;
  vector[K_version] gamma_version_raw;

  vector[I] b_raw;
  matrix[R, I] z_u;

  real<lower=1e-4, upper=3> sigma;
  real<lower=1e-4, upper=3> sigma_b;
  real<lower=1e-4, upper=3> tau_u;
  real<lower=1e-4, upper=3> sigma_version;
}

transformed parameters {
  vector[I] b;
  vector[K_version] gamma_version;
  matrix[R, I] u;

  b = sigma_b * b_raw;

  gamma_version = sigma_version * gamma_version_raw;
  gamma_version = gamma_version - mean(gamma_version);

  for (i in 1:I) {
    vector[R] temp_u;
    temp_u = tau_u * L_cov * z_u[, i];
    u[, i] = temp_u - mean(temp_u);
  }
}

model {
  vector[N] mu;

  alpha ~ normal(0, 1);
  beta_age ~ normal(0, 0.5);
  gamma_version_raw ~ normal(0, 1);

  b_raw ~ normal(0, 1);
  to_vector(z_u) ~ normal(0, 1);

  sigma ~ normal(0, 0.5);
  sigma_b ~ normal(0, 0.5);
  tau_u ~ normal(0, 0.5);
  sigma_version ~ normal(0, 0.5);

  for (n in 1:N) {
    mu[n] =
      alpha[region_id[n]] +
      beta_age[region_id[n]] * age_z[n] +
      gamma_version[fs_version_id[n]] +
      b[subject_id[n]] +
      u[region_id[n], subject_id[n]];
  }

  y ~ normal(mu, sigma);
}
