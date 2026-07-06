functions {
  // no custom functions in v1
}
data {
  int<lower=0> N;
  int<lower=1> D;
  matrix[N, D] Z;
  vector[D] mu0;
  corr_matrix[D] R0;
}
parameters {
  vector[D] mu;
  cholesky_factor_corr[D] L;
  vector<lower=0>[D] sigma;
}
transformed parameters {
  matrix[D, D] L_sigma = diag_pre_multiply(sigma, L);
}
model {
  // Priors shrunk toward Phase 1 estimates
  mu ~ normal(mu0, 0.1);
  L ~ lkj_corr_cholesky(2);
  sigma ~ cauchy(0, 2.5);
  for (n in 1:N) {
    Z[n] ~ multi_normal_cholesky(mu, L_sigma);
  }
}
generated quantities {
  vector[D] mu_post_mean;
  mu_post_mean = mu;
}
