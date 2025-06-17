functions {
//   real dirichlet_multinomial_lpmf(array[] int y, vector alpha) {
//     real alpha_0 = sum(alpha);
//     real n = sum(y);
//     return lgamma(alpha_0) - lgamma(n + alpha_0) + 
//       + lgamma(n+1) - sum(lgamma(to_vector(y)+1)) +
//       sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha));
//   }
  real dirichlet_multinomial_lpmf(array[,] int y_arr, matrix alpha) {
    int N = dims(y_arr)[1];
    int K = dims(y_arr)[2];
    matrix[N, K] y = to_matrix(y_arr);
    vector[K] v = rep_vector(1.0, K);
    vector[N] alpha_0 = alpha * v;
    vector[N] n = y * v;
    vector[N] y_rowsum_1  = lgamma(y+1) * v;
    vector[N] y_rowsum_alpha  = lgamma(y+alpha) * v;
    vector[N] alpha_rowsum = lgamma(alpha) * v;
    return sum(lgamma(alpha_0) - lgamma(n + alpha_0) + lgamma(n + 1) -
        y_rowsum_1 + y_rowsum_alpha - alpha_rowsum);
  }
//   real dirichlet_multinomial_lpdf(matrix y, matrix alpha) {
//     int N = dims(y)[1];
//     int K = dims(y)[2];
//     vector[K] v = rep_vector(1.0, K);
//     vector[N] alpha_0 = alpha * v;
//     vector[N] n = y * v;
//     vector[N] y_rowsum_1  = lgamma(y+1) * v;
//     vector[N] y_rowsum_alpha  = lgamma(y+alpha) * v;
//     vector[N] alpha_rowsum = lgamma(alpha) * v;
//     return sum(lgamma(alpha_0) - lgamma(n + alpha_0) + lgamma(n + 1) -
//         y_rowsum_1 + y_rowsum_alpha - alpha_rowsum);
//   }
}
data {
  int<lower=0> V; // # of baseline variants
  int<lower=0> S; // # of synonymous variants
  int<lower=0> N; // # of observation
  int<lower=0> N_syn; // # of syn observation
  int<lower=0> P; // # of positions
  int<lower=0> R; // # of reps
  int<lower=0> G; // # of groups
  array[N] int nMAPv;
  array[N] int nMAPr;
  array[N] int nMAPg;
  array[V] int vMAPp;
  array[V] int vMAPs;
//   array[N] int<lower=0> n_counts; // variant counts vector
  vector[N] n_counts; // variant counts vector
  int<lower=0> K; // # of bins
  array[N, K] int<lower=0> y; // FACS bin count

  array[N_syn] int sMAPr;
  array[N_syn] int sMAPg;
  vector[N_syn] n_syn_counts;
  array[N_syn, K] int<lower=0> y_syn; // FACS bin count
//   matrix<lower=0> [N, K]  y; // FACS bin count
//   simplex[K] q; // baseline 
}
// transformed data {
// //   vector[K - 1] t0; // latent cutoff
// //   t0 = inv_Phi(head(cumulative_sum(q), K - 1));
  
// }
parameters {
  simplex[K] q; 
  vector<lower=0> [P-1] sigma;
  real<lower=0> sigma_syn;
  vector[P-1] theta;
  vector[V] z; // shifted mu for each variant
//   vector<lower=0> [S] tau;
  vector [S] theta_syn;
  vector<lower=0> [G] phi_groups; // dispersion in dirichlet
//   vector<lower=0> [R] a;
//   vector<lower=0> [R] b;
}
transformed parameters {
  vector[V] mu; // shifted mu for each variant
  for (v in 1:V) {
    if (vMAPp[v] == 1) { // if synonymous, use variant specific tau as variance (~ cauchy (0, sigma_syn))
    //   mu[v] = theta[vMAPp[v]] + sigma[vMAPs[v]]*z[v];
    //   mu[v] = theta_syn[vMAPs[v]] + sigma_syn*z[v];
        mu[v] = theta_syn[vMAPs[v]] + sigma_syn*z[v];
    } else {
      mu[v] = theta[vMAPp[v]-1] + sigma[vMAPp[v]-1]*z[v];
    }
  }
}
model { 
  // priors
  theta ~ std_normal();
  theta_syn ~ std_normal();
  sigma_syn ~ inv_gamma(1, 1);
  sigma ~ inv_gamma(1, 1);
  phi_groups ~ exponential(0.01);
  z ~ std_normal();
  q ~ uniform(0, 1);

  // overdispersion
  vector[N] phi;
  phi = phi_groups[nMAPg];
//   phi = a[nMAPr] .* n_counts + b[nMAPr];
  matrix [N, K] phi_mat = rep_matrix(phi, K);

  vector[N_syn] phi_syn;
  phi_syn = phi_groups[sMAPg];
//   phi_syn = a[sMAPr] .* n_syn_counts + b[sMAPr];
  matrix [N_syn, K] phi_syn_mat = rep_matrix(phi_syn, K);
  // baseline 
  y_syn ~ dirichlet_multinomial(phi_syn_mat .* rep_matrix(to_row_vector(q), N_syn));
  vector[K - 1] t0; // latent cutoff
  t0 = inv_Phi(head(cumulative_sum(q), K - 1));

  // latent normal transform
  matrix[N, K-1] t2cdf;
  matrix[N, K] p0;
  t2cdf = Phi(rep_matrix(t0, N)' - rep_matrix(mu[nMAPv], K-1));
  p0[,1] = t2cdf[,1];
  for (k in 2:(K - 1)) {
    p0[,k] = t2cdf[,k] - t2cdf[,k - 1];
  }
  p0[,K] = 1 - t2cdf[,K - 1];

  // sample
  y ~ dirichlet_multinomial(phi_mat .* p0);
}
