library(ggplot2)
library(dplyr)
source('../src/build_sim2.R')
source('../src/run_sim.R')
source('../src/methods/method_utils.R')
source('../src/methods/diagnostics.R')
source('../src/plot_func.R')

args <- commandArgs(trailingOnly = TRUE)
stan_model_df_file <- args[1]
filename <- tools::file_path_sans_ext(basename(stan_model_df_file))
modelname <- substr(filename, 1, nchar(filename) - 3)
stan_model_df <- readRDS(stan_model_df_file)
output_df_VI_file <- args[2]
plot_dir <- args[3]
conda_path <- args[4]

## write table and run VI version
input_file <- paste0(plot_dir, "/pre_VI.csv")
write_csv(stan_model_df, input_file)


guide_type <- "normal"
output_file <- paste0(plot_dir, paste0("/VI_variant_scores_", guide_type, ".tsv"))

stdout_file <- tempfile()
stderr_file <- tempfile()

VI_diagnostics <- list()
VI_start <- Sys.time()
if (modelname == "FACS_double_sample_repq_nopos") {
  status <- system2(conda_path, args=c("run", "-n", "DMS", "python", "../src/run_model_VI.py",
                        "--input_file ", input_file,
                        "--output_file ", output_file,
                        "--guide_type ", guide_type,
                        "--plot_dir ", plot_dir,
                        "--nopos"),
        stdout = stdout_file,
        stderr = stderr_file)
} else {
  status <- system2(conda_path, args=c("run", "-n", "DMS", "python", "../src/run_model_VI.py",
                        "--input_file ", input_file,
                        "--output_file ", output_file,
                        "--guide_type ", guide_type,
                        "--plot_dir ", plot_dir),
        stdout = stdout_file,
        stderr = stderr_file)
}
VI_end <- Sys.time()
VI_difftime <- VI_end - VI_start
VI_runtime <- paste(VI_difftime, attr(VI_difftime, "units"))
VI_diagnostics$normal_runtime <- VI_runtime

stdout_output <- readLines(stdout_file)
stderr_output <- readLines(stderr_file)
unlink(c(stdout_file, stderr_file))

print(stdout_output)
print(stderr_output)

## write table and run multinormal VI version
input_file <- paste0(plot_dir, "/fitted_df_VI.csv")

guide_type <- "multivariate_normal"
output_file <- paste0(plot_dir, paste0("/multinormal_VI_variant_scores_", guide_type, ".tsv"))

stdout_file <- tempfile()
stderr_file <- tempfile()

VI_start <- Sys.time()
if (modelname == "FACS_double_sample_repq_nopos") {
  status <- system2(conda_path, args=c("run", "-n", "DMS", "python", "../src/run_model_VI.py",
                        "--input_file ", input_file,
                        "--output_file ", output_file,
                        "--guide_type ", guide_type,
                        "--plot_dir ", plot_dir,
                        "--nopos"),
        stdout = stdout_file,
        stderr = stderr_file)
} else {
  status <- system2(conda_path, args=c("run", "-n", "DMS", "python", "../src/run_model_VI.py",
                        "--input_file ", input_file,
                        "--output_file ", output_file,
                        "--guide_type ", guide_type,
                        "--plot_dir ", plot_dir),
        stdout = stdout_file,
        stderr = stderr_file)
}
VI_end <- Sys.time()
VI_difftime <- VI_end - VI_start
VI_runtime <- paste(VI_difftime, attr(VI_difftime, "units"))
VI_diagnostics$multinormal_runtime <- VI_runtime

stdout_output <- readLines(stdout_file)
stderr_output <- readLines(stderr_file)
unlink(c(stdout_file, stderr_file))

print(stdout_output)
print(stderr_output)

results_df_with_VI <- read_csv(paste0(plot_dir, "/fitted_df_VI.csv"))


saveRDS(results_df_with_VI, output_df_VI_file)
yaml::write_yaml(VI_diagnostics, paste0(plot_dir, "/VI_model_stats.yaml"))