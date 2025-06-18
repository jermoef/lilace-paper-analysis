library(ggplot2)
library(dplyr)
source('../src/build_sim.R')
source('../src/fit_func.R')
source('../src/methods/method_utils.R')
source('../src/methods/diagnostics.R')
source('../src/plot_func.R')


# cmd line args
args <- commandArgs(trailingOnly = TRUE)
modelfile <- args[1]
modelname <- tools::file_path_sans_ext(basename(modelfile))
baseline_model_file <- args[2]
data <- readRDS(args[3])
dataname <- tools::file_path_sans_ext(basename(args[3]))
baseline_data <- readRDS(args[4])
desc <- yaml::read_yaml(args[5])
use_selfphi <- as.logical(args[6])
output_fit_file <- args[7]
output_df_file <- args[8]
input_file <- args[9]
plot_dir <- args[10]

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

print(paste("using selfphi:", use_selfphi))

K <- sum(startsWith(colnames(data), "c_")) # num bins
G_size <- as.numeric(desc$G_size)
if (is.null(desc$G_size)) {
  G_size <- 100 # not used anyway unless running grouped model
}
print(paste("Using G size:", G_size))

# filter by > 15 total counts for variant
print(table(data$type))
data <- data %>% group_by(hgvs_exp) %>% mutate(total_counts=sum(n_counts))
data <- data[data$total_counts>=15,]
print(table(data$type))

# filter zero rows
print(paste("Removing", nrow(data) - nrow(data[data$n_counts > 0,]), "zero rows"))
data <- data[data$n_counts > 0,]
print(table(data$type))
# run baseline model
baseline_file <- paste0(tools::file_path_sans_ext(output_df_file), '_baseline.tsv')
base_stanfit <- fit_baseline(baseline_data %>% ungroup(), baseline_model_file, baseline_file)
saveRDS(base_stanfit, paste0(tools::file_path_sans_ext(output_df_file), '_baseline_stanfit.RData'))
param <- read_tsv(baseline_file)
q <- param[startsWith(param$param, "q"),2][[1]]
q <- q/sum(q)
print("baseline estimate:")
print(q)

# stop if baseline divergences
sampler_params <- get_sampler_params(base_stanfit, inc_warmup = FALSE)
mean_accept_stat <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
mean_step <- sapply(sampler_params, function(x) mean(x[, "stepsize__"]))
max_tree <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
mean_leapfrog <- sapply(sampler_params, function(x) sum(x[, "n_leapfrog__"]))
chain_divergences <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
mean_energy <- sapply(sampler_params, function(x) mean(x[, "energy__"]))
print("Num divergences in each baseline chain: ")
print(chain_divergences)
stopifnot(sum(chain_divergences)==0)

# plot n_eff and Rhat for baseline
p <- ggplot() + geom_histogram(aes(x=summary(base_stanfit)$summary[,"n_eff"]))
ggsave(p, file=paste0(plot_dir, "/baseline_n_eff.png"))
p <- ggplot() + geom_histogram(aes(x=summary(base_stanfit)$summary[,"Rhat"]))
ggsave(p, file=paste0(plot_dir, "/baseline_Rhat.png"))

# write model stat file
model_stats <- list(mean_accept_stat=mean_accept_stat, mean_step=mean_step, max_tree=max_tree, 
mean_leapfrog=mean_leapfrog, chain_divergences=chain_divergences, mean_energy=mean_energy)
yaml::write_yaml(model_stats, paste0(plot_dir, "/baseline_model_stats.yaml"))

# run baseline on same condition to get selfphi
selfphi <- NULL
if (use_selfphi) {
  if (grepl("fixedphi_byreplicate", modelfile)) { # get separate baseline phi on each replicate
    rep_names <- unique(data$rep)
    n_rep <- length(rep_names)
    selfphi <- rep(NA, n_rep)
    for (i in 1:n_rep) {
      repphi_file <- paste0(tools::file_path_sans_ext(output_df_file), '_selfphi', i, '.tsv')
      repphi_stanfit <- fit_baseline(data[data$rep==rep_names[i],] %>% ungroup(), baseline_model_file, repphi_file)
      param <- read_tsv(repphi_file)
      selfphi[i] <- param[startsWith(param$param, "phi"),2][[1]]
    }
  } else { # single phi across replicates in baseline
    selfphi_file <- paste0(tools::file_path_sans_ext(output_df_file), '_selfphi.tsv')
    selfphi_stanfit <- fit_baseline(data %>% ungroup(), baseline_model_file, selfphi_file)
    param <- read_tsv(selfphi_file)
    selfphi <- param[startsWith(param$param, "phi"),2][[1]]
  }
}
print(paste("selfphi value:", selfphi))

# fit full model
start <- Sys.time()
model_diagnostics <- fit_overall(data %>% ungroup(), baseline_file, modelfile, output_fit_file, output_df_file, input_file, G_size, selfphi)
end <- Sys.time()
difftime <- end - start
runtime <- paste(difftime, attr(difftime, "units"))
model_diagnostics$runtime <- runtime
chain_divergences <- model_diagnostics$num_divergent
print("Num divergences in each full model chain: ")
print(chain_divergences)
yaml::write_yaml(model_diagnostics, paste0(plot_dir, "/model_stats.yaml"))