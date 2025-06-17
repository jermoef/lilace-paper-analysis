
write_enrich2_files <- function(df, path, ref_type, condition_name, K=4) {
    n_rep <- length(unique(df$rep))
    # for each replicate, for each bin make tsv
    for (r in 1:n_rep) {
        rep_name <- unique(df$rep)[r]
        for (k in 1:K) {
            count <- as.data.frame(df[df$rep==rep_name, paste0("c_", 0:(K-1))][,k] + 1)
            if (ref_type == "synonymous") {
                count[nrow(count) + 1,] <- sum(count[df[df$rep==rep_name,"type"]=="synonymous",], na.rm = TRUE)
            } else {
                # supply wt
                count[nrow(count) + 1,] <- wt[k]
            }
            colnames(count) <- "count"
            rownames(count) <- c(df[df$rep==rep_name,]$hgvs_exp, "_wt")
            write.table(count, quote = FALSE, sep="\t", file = file.path(path, paste0("count_rep", rep_name, "_c", k-1, ".tsv")))
        }
    }
    # json format
    experiment <- list('name' = 'simulation', 'output directory' = path, 'conditions' = list())
    condition <- list('name' = condition_name, selections = list())
    for (r in 1:n_rep) {
        rep_name <- unique(df$rep)[r]
        FACS <- list('name' = paste("rep", rep_name, sep = ""), 'libraries' = list())
        for (k in 1:K) {
            seqlib <- list('counts file' = file.path(path, paste0("count_rep", rep_name, "_c", k-1, ".tsv")),
                        'identifiers' = list(),
                        'name' = paste0("rep", rep_name, "_c", k-1),
                        'report filtered read' = FALSE,
                        'timepoint' = k - 1)
            FACS$libraries[[k]] <- seqlib
        }
        condition$selections[[r]] <- FACS
    }
    experiment$conditions[[1]] <- condition

    ### save json file
    jsonData <- rjson::toJSON(experiment, indent = 2)
    write(jsonData, file = file.path(path, "config.json"))
}

run_enrich2 <- function(shell, rcfile, enrich_config) {
    system(
        paste0(shell, " -c
        . $HOME/", rcfile, "
        conda activate enrich2
        enrich_cmd ", enrich_config, " WLS wt --no-plots
        conda activate DMS
        ")
    )
}

run_enrich2_ratio <- function(shell, rcfile, enrich_config) {
    system(
        paste0(shell, " -c
        . $HOME/", rcfile, "
        conda activate enrich2
        enrich_cmd ", enrich_config, " ratios wt --no-plots
        conda activate DMS
        ")
    )
}

# read enrich2 data
extract_enrich_exp <- function(dir) {
  # scores
  file <- file.path(dir, "tsv/simulation_exp/main_identifiers_scores.tsv")
  enrich_score <- read_tsv(file, skip = 2, col_names = FALSE)

  # change column name
  title <- read_tsv(file, n_max = 2, col_names = FALSE)
  print(title)
  print(as.character(apply(as.matrix(title), 2, str_c, collapse = ".")))
  colnames(enrich_score) <- as.character(apply(as.matrix(title), 2, str_c, collapse = "."))
  colnames(enrich_score)[1] <- "hgvs_exp"
  rm(title)

  # clean the data frame
  print(enrich_score)
  print(enrich_score %>% pivot_longer(!hgvs_exp))
  enrich_score <- enrich_score %>%
    pivot_longer(!hgvs_exp) %>%
    separate(name, into = c("exp", "stats"), sep="[.]") %>%
    pivot_wider(names_from = stats, values_from = value)
  colnames(enrich_score)[3:5] <- str_c("enrich", colnames(enrich_score)[3:5], sep = "_")

  # pvalues
  file <- file.path(dir, "tsv/simulation_exp/main_identifiers_scores_pvalues_wt.tsv")
  enrich_pval <- read_tsv(file, skip = 2, col_names = FALSE)

  # change column name
  title <- read_tsv(file, n_max = 2, col_names = FALSE)
  colnames(enrich_pval) <- as.character(apply(as.matrix(title), 2, str_c, collapse = "."))
  colnames(enrich_pval)[1] <- "hgvs_exp"
  rm(title)

  # clean the data frame
  enrich_pval <- enrich_pval %>%
    pivot_longer(!hgvs_exp) %>%
    separate(name, into = c("exp", "stats"), sep = "[.]") %>%
    pivot_wider(names_from = stats, values_from = value)

  colnames(enrich_pval)[3:4] <- str_c("enrich", c("pval", "z"), sep = "_")

  # combine score and pval
  enrich_res <- full_join(enrich_pval, enrich_score)

  return(enrich_res)
}


gsub_models <- function(labels) {
  labels <- gsub("FACS_phi_widepos_syndist_nfunc_rep_syn_recalibrated_mu_mean", "FACS (fixed syn)", labels)
  labels <- gsub("FACS_phi_widepos_syngrouped_nfunc_rep_syn_recalibrated_mu_mean", "FACS (grouped syn)", labels)
  labels <- gsub("FACS_phi_widepos_nfunc_rep_syn_recalibrated_mu_mean", "FACS (syn reg)", labels)
  labels <- gsub("FACS_only_phi_nfunc_rep_syn_recalibrated_mu_mean", "FACS (no pos)", labels)
  labels <- gsub("enrich_score_self", "Enrich2", labels)
  return(labels)
}