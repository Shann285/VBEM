data {
int<lower=1> N; // sample size
int<lower=1> J; // number of items
array[N,J] int<lower=0, upper=1> Y; // data matrix of order [N,J]
int<lower=1> P; // number of covariates 
matrix[N,P] X;
real Nmupr;
real<lower=0> Nsigpr;
real alpha_a;
real alpha_b;
}

transformed data {
matrix[P, N] X_t = X';
array[N*J] int<lower=0, upper=1> Y_1d = to_array_1d(Y);
}

parameters {
vector[J] inu; // baseline intercept
vector<lower=0>[J] lam; // baseline slope
matrix[J,P] kap;  // intercept DIF coefficients
vector[N] fac_dist_helper; // helper for non-centered sampling
vector[P] mu_imp; // factor mean impact
vector[P] phi_imp; // factor sd impact
matrix<lower=0>[J,P] alpha_i; // Laplace variance on intercept DIF
}

model {
vector[N] fac_scor = X*mu_imp + exp(0.5*(X*phi_imp)) .* fac_dist_helper;
matrix[J,N] mu = rep_matrix(inu, N) + (kap*X_t) + (lam*fac_scor');

// the priors
to_vector(alpha_i) ~ gamma(alpha_a, alpha_b);
inu ~ normal(Nmupr, Nsigpr);
lam ~ normal(1, 0.5);
fac_dist_helper ~ std_normal();
to_vector(kap) ~ double_exponential(0, inv_sqrt(to_vector(alpha_i)));

mu_imp ~ normal(Nmupr, Nsigpr);
phi_imp ~ normal(Nmupr, Nsigpr);

// The likelihood
Y_1d ~ bernoulli_logit(to_vector(mu));
}
