functions {
  real dirichlet_multinomial_lpmf(array[] int y, vector alpha) {
    real alpha_0 = sum(alpha);
    real n = sum(y);
    return lgamma(alpha_0) - lgamma(n + alpha_0) + 
      + lgamma(n+1) - sum(lgamma(to_vector(y)+1)) +
      sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha));
  }
}
data {
  int<lower=0> N; // # of baseline variants
  int<lower=0> K; // # of bins
  array[N, K] int<lower=0> y; // FACS bin count
}
parameters { 
  real<lower=0> phi; // dispersion in dirichlet
  simplex[K] q; // baseline 
}
model { 
  phi ~ exponential(0.1);
  for (n in 1:N) {
    y[n] ~ dirichlet_multinomial(phi*q);
  }
}
