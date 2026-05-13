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
