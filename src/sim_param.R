source('../src/build_sim.R')
source('../src/fit_func.R')
source('../src/methods/method_utils.R')

start <- Sys.time()

args <- commandArgs(trailingOnly = TRUE)

# get sim parameters from real dataset
ref_data_file <- args[1]
yaml_file <- yaml::read_yaml(args[2])
# read in specific param to toggle
n_pos <- as.numeric(args[3])
param_string <- args[4]
param_value <- as.numeric(args[5]) # could also add prop LOF/GOF here
sim_data_outfile <- args[6]
sim_data_rescaled_outfile <- args[7]
sim_obj_outfile <- args[8]
iter <- args[9]
set.seed(iter)
# read from config
cell_count_factor <- 200
cell_count_phi <- -1 # -1 means equal cell samples / o.w. draw from dirmult
variant_phi <- -2 # -1 means use real samples o.w use this is replicate level overdispersion, if -2 estimate linear fxn for phi
FACS_gating <- "log_equal"
effect_prop <- -1 # -1 means empirical, otherwise draw from distribution
# baseline_effect_prop <- 0.3 # lower bound on # of effects (at least keeps this proportion of small effects)
reads_phi <- -1 # -1 means use real samples o.w. this is bin read overdispersion
position_effect <- NA # set
position_effect_mode <- "variance" # if sampling: 1 means uniform, 0 means point mass
breaks <- yaml_file$breaks
gate_offset <- 0 # shift of gate quantiles from equal proportions
prop_small <- NA
small_effect_sd <- NA
small_effect_cutoff <- NA
large_effect_sd <- NA
large_effect_cutoff <- NA
reads_per_obs <- NA
bin_phi <- -1
rep_effect <- 0

a <- NA
b <- NA

# grab parameters from fitted parameter dataset
print("getting parameters")
data <- readRDS(ref_data_file)
K <- sum(startsWith(colnames(data), "c_")) # num bins in data
n_bins <- K
total_pos <- length(unique(data$position))
n_variants <- length(unique(data$mutation))
n_replicates <- length(unique(data$rep))
n_reads <- NA


# variant bin variance (variance goes in all directions instead of just left-right shift from the effect variance)
data <- data[data$n_counts !=0, ] %>% group_by(hgvs_exp) %>% add_count(name="n_rep")
variant_phi_quantiles <- quantile(data[data$n_rep==3,]$phi, prob=c(0.1,0.9)) # take middle 80% of values to avoid outliers
variant_phi_samples <- data[data$n_rep==3 & data$phi > variant_phi_quantiles[1] & data$phi < variant_phi_quantiles[2],]$phi

# get ML effect distiributions
library(mclust)
data_effect_estim <- data[!is.na(data$ML_effect_mean),]
gmm <- Mclust(data_effect_estim$ML_effect_mean, G=2)
data_effect_estim$ML_class <- gmm$classification

syn_labels <- data_effect_estim[data_effect_estim$type=="synonymous",]$ML_class
syn_class <- sort(table(syn_labels), decreasing=T)[1]
syn_class_label <- names(syn_class)
print(syn_class_label)

if (is.na(small_effect_cutoff)) {
  # p small cutoff
  small_effects <- data_effect_estim[data_effect_estim$ML_class == syn_class_label,]$ML_effect_mean
  pars <- MASS::fitdistr(small_effects, "normal")
  small_effect_mean <- pars$estimate[1]
  small_effect_sd <- pars$estimate[2]
  small_effect_cutoff <- min(abs(c(small_effect_mean + 2*small_effect_sd, small_effect_mean - 2*small_effect_sd)))
}
# proportion of small effects
if (is.na(prop_small)) {
  prop_small <- (table(data_effect_estim$ML_class) / nrow(data_effect_estim))[syn_class_label]
}

if (is.na(large_effect_cutoff)) {
  # p large cutoff
  large_effects <- data_effect_estim[data_effect_estim$ML_class != syn_class_label,]$ML_effect_mean
  pars <- MASS::fitdistr(large_effects, "normal")
  large_effect_mean <- pars$estimate[1]
  large_effect_sd <- pars$estimate[2]
  large_effect_cutoff <- min(abs(c(large_effect_mean + 2*large_effect_sd, large_effect_mean - 2*large_effect_sd)))
  p_large <- 1 - prop_small
}

# variance in effect at a position
if (is.na(position_effect)) {
  position_effect <- mean((data_effect_estim %>% group_by(position) %>% mutate(pos_sd=sd(ML_effect_mean)))$pos_sd)
}

# get position effect pattern from weighted bin mean from data using 2*sd cutoff
print(head(data))
pos_mut_df <- unique(data[,c("position", "mutation", "weight_effect_mean", "type")])
pos_mut_df$position <- as.numeric(pos_mut_df$position)
# sample random n_pos positions and renumber as 1:n_pos
sampled_positions <- sample(unique(pos_mut_df$position), n_pos)
pos_mut_df <- pos_mut_df[pos_mut_df$position %in% sampled_positions,]

# set positions to numeric
pos_mut_df <- pos_mut_df %>% 
  group_by(position) %>% 
  mutate(position = cur_group_id()) 
# set mutations to numeric
pos_mut_df <- pos_mut_df %>% 
  group_by(mutation) %>% 
  mutate(mutation = cur_group_id()) 
weight_syn_mean <- mean(pos_mut_df[pos_mut_df$type=="synonymous",]$weight_effect_mean)
weight_syn_sd <- sd(pos_mut_df[pos_mut_df$type=="synonymous",]$weight_effect_mean)
pos_mut_df$is_effect <- (pos_mut_df$weight_effect_mean > weight_syn_mean + 2 * weight_syn_sd) | (pos_mut_df$weight_effect_mean < weight_syn_mean - 2 * weight_syn_sd)
print(pos_mut_df)


# get syn fluorescence distribution from syn counts
syn_counts <- data[data$type=="synonymous",paste0("c_",0:(K-1))]
fit_syn_direct <- optim(c(1000, 6.5), function(theta) cdf_lik_log_multi(theta, as.matrix(syn_counts), breaks, cdf=plnorm))
syn_mean <- fit_syn_direct$par[1]
measurement_variance <- exp(fit_syn_direct$par[2])

# variance in reads 
bin_reads_phi_means <- c()
n_counts_phi_means <- c()
for (r in 1:n_replicates) {
  rep <- unique(data$rep)[r]
  bin_phis <- array(NA, c(100,K))
  n_counts_phis <- rep(NA, 100)
  for (i in 1:100) {
    # positions <- sample(data[data$rep==rep,]$position, n_pos, replace=T) # randomly select number of pos being simmed and get variant overdsps
    positions <- sample(data[data$rep==rep,]$position, total_pos, replace=T) # randomly select number of pos being simmed and get variant overdsps
    for (k in 1:K) {
      n_counts <- data[data$position %in% positions & data$rep==rep,][[paste0("c_",k-1)]]
      num_variants <- length(n_counts)
      bin_phis[i,k] <- optim(10, function(phi) {dirmult_NLL_onlyphi(phi, num_variants, n_counts)}, method="Brent", lower=0, upper=1e10)$par
    }
    # sample n count phi
    n_counts <- data[data$position %in% positions & data$rep==rep,]$n_counts
    num_variants <- length(n_counts)
    n_counts_phis[i] <- optim(10, function(phi) {dirmult_NLL_onlyphi(phi, num_variants, n_counts)}, method="Brent", lower=0, upper=1e10)$par
  }
  for (k in 1:K) {
    bin_reads_phi_means <- rbind(bin_reads_phi_means, c(r, k, mean(bin_phis[,k]), sd(bin_phis[,k])))
  }
  n_counts_phi_means <- rbind(n_counts_phi_means, c(r, mean(n_counts_phis), sd(n_counts_phis)))
}
bin_reads_phi_means <- as.data.frame(bin_reads_phi_means)
colnames(bin_reads_phi_means) <- c("rep", "bin", "mean_eta", "sd_eta")
print(bin_reads_phi_means)
n_counts_phi_means <- as.data.frame(n_counts_phi_means)
colnames(n_counts_phi_means) <- c("rep", "mean_eta", "sd_eta")
print(n_counts_phi_means)

# estimate read phi corr
data_3rep <- data
data_3rep <- data_3rep %>% group_by(hgvs) %>% mutate(n_rep = n())
data_3rep <- data_3rep[data_3rep$n_rep==n_replicates,]
print(head(data_3rep))
print(table(data_3rep$rep))
print(table(data_3rep$n_rep))
corrs <- c()
rep_string <- sort(unique(data_3rep$rep))[1:n_replicates]
rep_combo <- combn(rep_string, 2)
for (i in 1:ncol(rep_combo)) {
  rep1 <- rep_combo[,i][1]
  rep2 <- rep_combo[,i][2]
  print(rep1)
  print(rep2)
  for (bin in paste0("c_", 1:K - 1)) {
    print(bin)
    print(head(data_3rep[data_3rep$rep==rep1,][[bin]]))
    print(length(data_3rep[data_3rep$rep==rep1,][[bin]]))
    print(head(data_3rep[data_3rep$rep==rep2,][[bin]]))
    print(length(data_3rep[data_3rep$rep==rep2,][[bin]]))
    cor <- cor(data_3rep[data_3rep$rep==rep1,][[bin]], data_3rep[data_3rep$rep==rep2,][[bin]])
    corrs <- rbind(corrs, c(rep1, rep2, bin, cor))
  }
}
corr_df <- data.frame(corrs)
colnames(corr_df) <- c("rep1", "rep2", "bin", "crossrep_cor")
mean_crossrep_cor <- mean(as.numeric(corr_df$crossrep_cor))

# read count
if (is.na(n_reads)) {
  n_reads <- sum(data$n_counts) / n_replicates * n_pos / total_pos # scale # reads per replicate to number of positions
} else {
  n_reads <- n_reads * n_pos * n_variants # read coverage * variants to get reads per replicate
}

# if param string not empty - set parameter to param
print(paste("assigning", param_string, param_value))
assign(param_string, param_value)

# cell count
n_cells <- n_pos * n_variants * cell_count_factor
cell_counts <- rep(cell_count_factor, n_pos*n_variants)


cat(paste("prop_effect", effect_prop, "\n", "measurement_variance", measurement_variance, "\n", 
    "position_effect", position_effect, "\n", "cell_count_factor", cell_count_factor, "\n",
    "variant_phi", variant_phi, "\n", "reads_phi", reads_phi, "\n", "bin_phi", bin_phi, "\n", "K", K))

print("Sampling effect sizes...")

n_data_variants <- max(length(unique(pos_mut_df$mutation)), n_variants)
is_effect_mat <- matrix(FALSE, n_data_variants, n_pos)
print(dim(is_effect_mat))
print(head(pos_mut_df))
print(sort(unique(pos_mut_df$position)))
print(sort(unique(pos_mut_df$mutation)))
is_effect_mat[as.matrix(pos_mut_df[2:1])] <- pos_mut_df$is_effect
is_effect_mat <- is_effect_mat[1:n_variants,] # if diff number of variants set (hopefully less than in the datset lol)
print(is_effect_mat)
n_effects <- sum(is_effect_mat)
print(n_effects)
# build effect distribution based on low and high prop
n_small <- round(prop_small * n_effects)
print(n_small)
n_large <- n_effects - round(prop_small * n_effects)
small_effects <- rnorm(10000, small_effect_mean, small_effect_sd)
small_effects <- sample(small_effects[abs(small_effects) < small_effect_cutoff], n_small, replace=T)
large_effects <- rnorm(10000, large_effect_mean, large_effect_sd)
large_effects <- sample(large_effects[abs(large_effects) > large_effect_cutoff], n_large, replace=T)
effects <- c(small_effects, large_effects)

# sample effect sizes 
# sample effect size distribution by taking sampled effect size and rejecting as effect 
# with probability according to likelihood of coming from synonymous fitted distr
N <- n_pos * n_variants
i <- 1
position_effect_indices <- sample.int(length(effects), n_pos, replace=T)
sampled_effect_sizes <- rep(NA, N)
for (p in 1:length(position_effect_indices)) {
    i <- position_effect_indices[p]
    if (position_effect_mode == "sampling") {
        # position effect = 1 means uniform sampling -> 0 for point mass same effect
        effects_sorted <- sort(effects)
        index_probs <- sapply(1:length(effects), function(j) {position_effect^(abs(j-i)*0.001)})
        index_probs <- index_probs / sum(index_probs) # normalize to probability
        pos_effects <- sample(effects_sorted, n_variants, replace=T, prob=index_probs)
    } else if (position_effect_mode == "variance") {
        pos_effects <- rnorm(n_variants, effects[i], position_effect)
    } else {
        stop("position effect mode not recognized")
    }
    pos_effects[1] <- 0 # synonymous mutation zero effect
    sampled_effect_sizes[((p-1)*n_variants+1):(p*n_variants)] <- pos_effects
}
effect_mat <- matrix(sampled_effect_sizes, n_variants, n_pos)
if (effect_prop > 0) {
    # take is effect mat and randomly increase or decrease effects to match effect prop
    n_effects <- sum(is_effect_mat)
    n_target <- round(effect_prop * length(effect_mat))
    if (n_effects < n_target) {
        set_vec <- rep(F, sum(!is_effect_mat))
        set_vec[1:(n_target - n_effects)] <- T
        is_effect_mat[!is_effect_mat] <- sample(set_vec)
    } else if (n_effects > n_target) {
        set_vec <- rep(T, sum(is_effect_mat))
        set_vec[1:(n_effects - n_target)] <- F
        is_effect_mat[is_effect_mat] <- sample(set_vec)
    }
}
# set non effects to zero
effect_mat[!is_effect_mat] <- 0

print(do.call(paste, as.list(c("Effect mat dim:", dim(effect_mat)))))
print(paste("Prop small effects:", sum(abs(effect_mat[effect_mat != 0]) < small_effect_cutoff)/sum(effect_mat != 0)))
print(paste("Prop large effects:", sum(abs(effect_mat[effect_mat != 0]) > large_effect_cutoff)/sum(effect_mat != 0)))
print(paste("Prop zero:", sum(effect_mat == 0)/ length(effect_mat)))


# if set to -2, estimate linear fxn of count dsp from data
if (variant_phi == -2) {
    syn_data <- data[data$type=="synonymous" & data$n_rep==max(data$n_rep),]
    data_phi_quant <- quantile(syn_data$phi, probs=0.9) # trim outliers
    phi_fit <- lm(phi ~ rep_avg_n, data=syn_data[syn_data$phi < data_phi_quant,])
    a <- phi_fit$coefficients[2]
    b <- phi_fit$coefficients[1]
}


# rep effects
rep_effects <- seq(-1*rep_effect, rep_effect, length.out=n_replicates)
print(rep_effects)
print("Running simulation...")
# create sim object
sim_obj <- create_sim(n_pos, n_variants, n_replicates, cell_counts, cell_count_factor,
                        effect_mat, syn_mean, measurement_variance, variant_phi, variant_phi_samples, a, b,
                        n_reads, reads_phi, bin_phi, bin_reads_phi_means, n_counts_phi_means, reads_per_obs, mean_crossrep_cor, FACS_gating, gate_offset, n_bins, rep_effects)
# sim data
sim_obj_simmed <- run_sim(sim_obj)
FACS_df <- sim_obj_simmed$FACS_df
obs_df <- sim_obj_simmed$obs_df
print(colnames(obs_df))
obs_df$n_counts <- apply(obs_df[paste0("c_",0:(n_bins-1))], 1, sum)
obs_df$exp <- "exp1"
# obs_df$hgvs_exp <- obs_df$hgvs
obs_df <- obs_df %>% unite("hgvs_exp", hgvs, exp, remove = FALSE)
saveRDS(obs_df, sim_data_outfile)
yaml::write_yaml(list(breaks=sim_obj_simmed$cutoffs), paste0(tools::file_path_sans_ext(sim_data_outfile), '.yaml'))
saveRDS(sim_obj_simmed, sim_obj_outfile)


# save rescaled RDS
count_df <- sim_obj_simmed$FACS_df %>% dplyr::select(starts_with("c_")) 
FACS_trace_percent <- colSums(count_df) / sum(count_df) # cell count proportions
print("Rescaling...")
print(FACS_trace_percent)
for (rep in unique(obs_df$rep)) {
    rep_total <- sum(obs_df[obs_df$rep==rep,paste0("c_",0:(n_bins-1))])
    for (i in 1:n_bins) {
        bin <- paste0("c_", 0:(n_bins-1))[i]
        rescaled_bin_total <- FACS_trace_percent * rep_total
        obs_df[obs_df$rep==rep,][[bin]] <- ceiling(obs_df[obs_df$rep==rep,][[bin]] / sum(obs_df[obs_df$rep==rep,][[bin]]) * rescaled_bin_total[i])
    }
}
  
obs_df$n_counts <- apply(obs_df[paste0("c_",0:(n_bins-1))], 1, sum)
# obs_df <- obs_df %>% group_by(hgvs) %>% mutate(total_counts=sum(n_counts))
saveRDS(obs_df, sim_data_rescaled_outfile)
yaml::write_yaml(list(breaks=sim_obj_simmed$cutoffs), paste0(tools::file_path_sans_ext(sim_data_rescaled_outfile), '.yaml'))

end <- Sys.time()
print(end - start)