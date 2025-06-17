library(dplyr)
library(yaml)


args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
input_yaml <- args[2]
output_file <- args[3]
output_yaml <- args[4]

prop_masked <- as.numeric(args[5])
iter <- as.numeric(args[6])

set.seed(iter)

data_exp <- readRDS(input_file)
desc <- yaml::read_yaml(input_yaml)
data_exp_masked <- data_exp %>% group_by(hgvs_exp) %>% mutate(var_id = cur_group_id())
syn_ids <- unique(data_exp_masked[data_exp_masked$type == "synonymous",]$var_id)
masked_syn_ids <- sample(syn_ids, ceiling(length(syn_ids) * prop_masked))
data_exp_masked[data_exp_masked$var_id %in% masked_syn_ids,]$type <- "masked_synonymous"
data_exp_masked <- data_exp_masked %>% ungroup()
saveRDS(data_exp_masked, output_file)
yaml::write_yaml(desc, output_yaml)
