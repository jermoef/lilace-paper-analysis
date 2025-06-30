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

K <- sum(startsWith(colnames(data), "c_")) # num bins
G_size <- as.numeric(desc$G_size)
if (is.null(desc$G_size)) {
  G_size <- 100 # not used anyway unless running grouped model
}

# filter by > 15 total counts for variant
print(table(data$type))
data <- data %>% group_by(hgvs_exp) %>% mutate(total_counts=sum(n_counts))
data <- data[data$total_counts>=15,]
print(table(data$type))

# filter zero rows
print(paste("Removing", nrow(data) - nrow(data[data$n_counts > 0,]), "zero rows"))
data <- data[data$n_counts > 0,]
print(table(data$type))

# fit full model
start <- Sys.time()
model_diagnostics <- fit_overall(data %>% ungroup(), modelfile, output_fit_file, output_df_file, input_file, G_size, selfphi)
end <- Sys.time()
difftime <- end - start
runtime <- paste(difftime, attr(difftime, "units"))
model_diagnostics$runtime <- runtime
chain_divergences <- model_diagnostics$num_divergent
print("Num divergences in each full model chain: ")
print(chain_divergences)
yaml::write_yaml(model_diagnostics, paste0(plot_dir, "/model_stats.yaml"))