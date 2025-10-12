source('../src/build_sim.R')
source('../src/fit_func.R')
source('../src/plot_func.R')
library(mclust)

args <- commandArgs(trailingOnly = TRUE)

data <- readRDS(args[1]) %>% ungroup()
desc <- yaml::read_yaml(args[2])
output_file <- args[3]
yaml_file <- args[4]


K <- sum(startsWith(colnames(data), "c_")) # num bins in data
n_bins <- K
# TODO: generalize to not 4 bins

breaks <- desc$breaks
vars <- unique(data$hgvs_exp)
data$ML_mean <- NA
data$ML_sd <- NA
# fit normal ML on each 3 rep variant 
for (i in 1:length(vars)) {
  var <- vars[i]
  var_df <- data[data$hgvs_exp==var,]
  counts <- var_df %>% dplyr::select(starts_with("c_"))
  params <- c(NA, NA)
  opt_out <- optim(c(1000, 6.5), function(theta) cdf_lik_log_multi(theta, counts, breaks, cdf=plnorm))
  if (opt_out$convergence == 0) {
    params <- opt_out$par
  } 
  data[data$hgvs_exp==var,]$ML_mean <- params[1]
  data[data$hgvs_exp==var,]$ML_sd <- exp(params[2])
}

data <- data[!is.na(data$ML_mean),]
data_control <- data[data$type=="synonymous",]
control_counts <- data_control %>% ungroup() %>% dplyr::select(starts_with("c_")) 
mean_control_counts <- colMeans(control_counts)
control_ML_mean <- mean(data_control$ML_mean)
control_ML_mean_var <- var(data_control$ML_mean)

data$ML_effects <- NA
for (i in 1:length(vars)) {
  var <- vars[i]
  var_df <- data[data$hgvs_exp==var,]
  counts <- var_df[c("c_0", "c_1", "c_2", "c_3")]
  ML_effects <- var_df$ML_mean - control_ML_mean
  ML_effect_mean <- mean(ML_effects)
  ML_effect_se <- sd(ML_effects)
  ML_effect_se_analytic <- mean(sqrt(var_df$ML_sd^2 / var_df$n_counts + mean(data_control$ML_sd^2 / data_control$n_counts))) # biased estimator of sd (underestimate), also not taking into account variance in control var estimate
  data[data$hgvs_exp==var,"ML_effects"] <- ML_effects
  data[data$hgvs_exp==var,c("ML_effect_mean","ML_effect_se", "ML_effect_se_analytic")] <- matrix(c(ML_effect_mean, ML_effect_se, ML_effect_se_analytic), nrow=nrow(var_df), ncol=3, byrow=T)

  # get data weight effects
  coeffs <- seq(1/K, 1, length.out=K)
  control_weights <- as.numeric(as.matrix(control_counts) %*% coeffs / data_control$n_counts)
  # control_weights <- as.numeric(unlist((0.25*control_counts[,1] + 0.5*control_counts[,2] +
  #                      0.75*control_counts[,3] + control_counts[,4]) / data_control$n_counts))
  control_weight_mean <- mean(control_weights, na.rm=T)
  weights <- as.numeric(as.matrix(counts) %*% coeffs / var_df$n_counts)
  # weights <- as.numeric(unlist((0.25*counts[,1] + 0.5*counts[,2] + 0.75*counts[,3] + counts[,4]) / var_df$n_counts))
  weight_effects <- weights - control_weight_mean
  weight_effect_mean <- mean(weight_effects, na.rm=T)
  weight_effect_se <- sd(weight_effects, na.rm=T)
  data[data$hgvs_exp==var,"weight_effects"] <- weight_effects
  data[data$hgvs_exp==var,c("weight_effect_mean","weight_effect_se")] <- matrix(c(weight_effect_mean, weight_effect_se), nrow=nrow(var_df), ncol=2, byrow=T)
}


syn_counts <- data[data$type=="synonymous",] %>% ungroup() %>% select(starts_with("c_"))
print(head(syn_counts))
opt_out <- optim(c(1000, 6.5), function(theta) cdf_lik_log_multi(theta, as.matrix(syn_counts), breaks, cdf=plnorm))
# estimate LOF effect mean + sd, GOF effect mean + sd
effects <- data$ML_effect_mean

# sample effect size distribution by taking sampled effect size and rejecting as effect 
# with probability according to likelihood of coming from synonymous fitted distr
N <- 10000
dx <- density(effects)
i <- 1
sampled_effect_sizes <- rep(NA, N)
x <- sample(effects, N, replace=T)
# compute p-val of observing under synonymous
fit_syn <- MASS::fitdistr(effects[data$type=="synonymous"], densfun="normal")
p_syn <- 2*pnorm(abs((x-fit_syn$estimate["mean"])/fit_syn$estimate["sd"]), lower.tail=F)
# reject with probability of observance
sampled_effect_sizes <- ifelse(runif(N) > p_syn, x, 0)

prop_LOF <- sum(sampled_effect_sizes < 0) / N
prop_GOF <- sum(sampled_effect_sizes > 0) / N

# variant level overdispersion
data[c("q_0", "q_1", "q_2", "q_3", "phi", "rep_avg_n", "rep_total_n")] <- NA
vars <- unique(data$hgvs_exp)

for (i in 1:length(vars)) {
  if (i %% 1000 == 0) {
    print(i)
  }
  var <- vars[i]
  var_df <- data[data$hgvs_exp==var,]
  counts <- var_df[c("c_0", "c_1", "c_2", "c_3")]
  avg_n <- mean(var_df$n_counts)
  total_n <- sum(var_df$n_counts)
  opt_out <- optim(c(1,1,1,1), function(alpha) dirmult_NLL_constrained(alpha, counts))
  phi <- sum(exp(opt_out$par))
  q <- exp(opt_out$par) / phi
  data[data$hgvs_exp==var,c("q_0", "q_1", "q_2", "q_3", "phi", "rep_avg_n", "rep_total_n")] <- matrix(c(q, phi, avg_n, total_n), nrow=nrow(var_df), ncol=7, byrow=T)
}

# estimate synonymous phi more consistently by fixing q
data_phi <- data[data$n_counts != 0, ] %>% group_by(hgvs_exp) %>% add_count(name="n_rep")
vars <- unique(data_phi[data_phi$type=="synonymous" & data_phi$n_rep == max(data_phi$n_rep),]$hgvs_exp)
syn_q <- colMeans((data_phi[data_phi$type=="synonymous",] %>% ungroup %>% select(starts_with("c_"))) / (data_phi[data_phi$type=="synonymous",]$n_counts))

for (i in 1:length(vars)) {
  if (i %% 100 == 0) {
    print(i)
  }
  var <- vars[i]
  var_df <- data_phi[data_phi$hgvs_exp==var,]
  counts <- var_df[c("c_0", "c_1", "c_2", "c_3")]
  avg_n <- mean(var_df$n_counts)
  total_n <- sum(var_df$n_counts)
  opt_out <- optim(c(1), function(phi) dirmult_NLL_fixed_alpha(phi, syn_q, counts), method="Brent", lower=0, upper=1000)
  phi <- opt_out$par
  data_phi[data_phi$hgvs_exp==var,c("q_0", "q_1", "q_2", "q_3", "phi", "rep_avg_n", "rep_total_n")] <- matrix(c(q, phi, avg_n, total_n), nrow=nrow(var_df), ncol=7, byrow=T)
}


# estimate position effect variance
data <- data %>% group_by(position, rep) %>% mutate(pos_sd=sd(ML_effect_mean)) %>% ungroup()


# saved vals
# data w/ estimated params (effects, variant level overdispersion estimates)
saveRDS(data, output_file)
# prop_LOF and GOF
yaml::write_yaml(list(prop_LOF=prop_LOF, prop_GOF=prop_GOF, breaks=desc$breaks), yaml_file)