

#' Computes posterior predictive variant obs outside CI counts for a given dataset, posterior draws, and model
posterior_samples <- function(n_sims, data, sim, modelname, K=4) {
    n_draws <- length(sim$lp__)
    if (modelname == "FACS_phi_pos_grouped") {
        return(posterior_CI_FACS_phi_pos_grouped(n_sims, n_draws, data, sim, K=K))
    }
    print("Posterior Predictive Model name not recognized")
    return(c())
}


posterior_CI_FACS_phi_pos_grouped <- function(n_sims, n_draws, data, sim, K) {
    scores <- c()
    y_rep <- array(NA, c(nrow(data), n_sims, K))
    for (i in 1:nrow(data)) {
        if (i %% 100 == 0) {
            print(i)
        }
        score <- 0
        s_indices <- sample.int(n_draws, n_sims, replace=T)
        for (j in 1:n_sims) {
            s <- s_indices[j]
            phi <- sim$phi[s]
            p <- brms::rdirichlet(1, phi*sim$p0[s,i,])
            y_rep[i,j,] <- rmultinom(1, data$n_counts[i], p)
        }
        for (k in 1:K) {
            bounds <- quantile(y_rep[i,,k], probs=c(0.05, 0.95))
            if (data[i,paste0("c_",k-1)] < bounds[1] || data[i,paste0("c_",k-1)] > bounds[2]) {
                score <- score+1
            }
        }
        scores <- c(scores, score)
    }
    return(scores)
}

# calling significance fxn
label_significant <- function(df, effect_cols, sd_cols, c_thresh) {
  for (i in 1:length(effect_cols)) {
    method <- effect_cols[i]
    if (method %in% colnames(df)) {
      if (grepl("enrich_score", method, fixed=T)) {
        # adjust enrich p-pvalues
        # print(hist(df$enrich_pval))
        df$enrich_pval_adjusted <- p.adjust(df$enrich_pval, method="BH")
        # print(hist(df$enrich_pval_adjusted))
        df[paste0(method, "_disc")] <- df$enrich_pval_adjusted <= 1 - c_thresh
        df[paste0(method, "_isnull")] <- df$enrich_pval_adjusted > 1 - c_thresh
        # df[paste0(method, "_rank")] <- rank(df$enrich_pval_adjusted, ties.method="first")
        rank_name <- paste0(method, "_rank")
        rank_fit_df <- df %>% group_by(hgvs) %>% filter(row_number() == 1) %>% ungroup() %>% 
          mutate(!! rank_name := dense_rank(enrich_pval_adjusted)) %>% select(hgvs, one_of(rank_name))
        df <- merge(df, rank_fit_df)
        # df[paste0(method, "_rank")] <- rank(-abs(df$effect))
        # add syn sd cut
        method_sd_cut <- paste0(method, "_syn_sd")
        syn_muts <- df[df$type=="synonymous", method]
        dist <- MASS::fitdistr(syn_muts[!is.na(syn_muts)], "normal")
        syn_mean <- dist$estimate[1]
        syn_sd <- dist$estimate[2]
        df[paste0(method_sd_cut, "_disc")] <- 2*pnorm(-abs((df[[method]] - syn_mean)/syn_sd)) < 1 - c_thresh
        df[paste0(method_sd_cut, "_isnull")] <- !df[paste0(method_sd_cut, "_disc")]
        rank_name <- paste0(method_sd_cut, "_rank")
        df[paste0(method_sd_cut, "_pval")] <- 2*pnorm(-abs((df[[method]] - syn_mean)/syn_sd))
        rank_fit_df <- df %>% group_by(hgvs) %>% filter(row_number() == 1) %>% ungroup() %>% 
          mutate(!! rank_name := dense_rank(get(paste0(method_sd_cut, "_pval")))) %>% select(hgvs, one_of(rank_name))
        df <- merge(df, rank_fit_df)
        # df[paste0(method_sd_cut, "_rank")] <- rank(2*pnorm(-abs((df[[method]] - syn_mean)/syn_sd)))
      } else if (grepl("mu_mean", method, fixed=T)) {
        model <- strsplit(method, "_mu_mean")[[1]]
        lfsr <- df[[paste0(model, "_mu_lfsr")]]
        df[paste0(method, "_disc")] <- lfsr <= 1 - c_thresh
        df[paste0(method, "_isnull")] <- lfsr > 1 - c_thresh

        print(model)
        df[[paste0(model, "_mu_lfsr_approx")]] <-  sapply(1:nrow(df), function(i) {
          p <- pnorm(0, mean = df[[i,paste0(model, "_mu_mean")]], sd = df[[i,paste0(model, "_mu_sd")]])
          lfsr <- min(p, 1-p)
          return(lfsr)
        })
        rank_name <- paste0(method, "_rank")
        rank_fit_df <- df %>% group_by(hgvs) %>% filter(row_number() == 1) %>% ungroup() %>% 
          mutate(!! rank_name := dense_rank(get(paste0(model, "_mu_lfsr_approx")))) %>% select(hgvs, one_of(rank_name))
        df <- merge(df, rank_fit_df)
        # df[paste0(method, "_rank")] <- rank(df[[paste0(model, "_mu_lfsr_approx")]], ties.method="first")

    #   } else if (method == "ML_effect_mean") {
    #     print("Using ML pval")
    #     df[paste0(method, "_disc")] <- df$ML_pval < 1 - c_thresh
    #     df[paste0(method, "_isnull")] <- df$ML_pval > 1 - c_thresh
    #     df[paste0(method, "_rank")] <- rank(df$ML_pval)
    #   } else if (method == "shep_effect_mean") {
    #     print("Using shep pval")
    #     df[paste0(method, "_disc")] <- df$shep_pval < 1 - c_thresh
    #     df[paste0(method, "_isnull")] <- df$shep_pval > 1 - c_thresh
    #     df[paste0(method, "_rank")] <- rank(df$shep_pval)
    #   } 
      } else if (is.na(sd_cols[i])) {
        # fit normal to synonymous to get syn mean and syn sd
        syn_muts <- df[df$type=="synonymous", method]
        dist <- MASS::fitdistr(syn_muts[!is.na(syn_muts)], "normal")
        syn_mean <- dist$estimate[1]
        syn_sd <- dist$estimate[2]
        # df[paste0(method_sd_cut, "_disc")] <- (mean < syn_mean - 2*syn_sd) | (mean > syn_mean + 2*syn_sd)
        df[paste0(method_sd_cut, "_disc")] <- 2*pnorm(-abs((mean - syn_mean)/syn_sd)) < 1 - c_thresh
        df[paste0(method_sd_cut, "_isnull")] <- !df[paste0(method_sd_cut, "_disc")]
        df[paste0(method_sd_cut, "_pval")] <- 2*pnorm(-abs((df[[method]] - syn_mean)/syn_sd))
        rank_name <- paste0(method_sd_cut, "_rank")
        rank_fit_df <- df %>% group_by(hgvs) %>% filter(row_number() == 1) %>% ungroup() %>% 
          mutate(!! rank_name := dense_rank(get(paste0(method_sd_cut, "_pval")))) %>% select(hgvs, one_of(rank_name))
        df <- merge(df, rank_fit_df)
        # df[paste0(method_sd_cut, "_rank")] <- rank(2*pnorm(-abs((mean - syn_mean)/syn_sd)))
      } else {
        # lfsr analgoues
        mean <- df[[effect_cols[i]]]
        sd <- df[[sd_cols[i]]]
        lfsr <- apply(cbind(pnorm(0, mean, sd), 1-pnorm(0, mean, sd)), 1, min)
        df[paste0(method, "_disc")] <- lfsr <= 1 - c_thresh
        df[paste0(method, "_isnull")] <- lfsr > 1 - c_thresh
        df[paste0(method, "_lfsr")] <- lfsr
        rank_name <- paste0(method, "_rank")
        rank_fit_df <- df %>% group_by(hgvs) %>% filter(row_number() == 1) %>% ungroup() %>% 
          mutate(!! rank_name := dense_rank(get(paste0(method, "_lfsr")))) %>% select(hgvs, one_of(rank_name))
        df <- merge(df, rank_fit_df)
        # df[paste0(method, "_rank")] <- rank(lfsr, ties.method="first")
        method_sd_cut <- paste0(method, "_syn_sd")
        # syn_mean <- mean(df[df$type=="synonymous", method])
        # syn_sd <- sd(df[df$type=="synonymous", method])

        # fit normal to synonymous to get syn mean and syn sd
        syn_muts <- df[df$type=="synonymous", method]
        dist <- MASS::fitdistr(syn_muts[!is.na(syn_muts)], "normal")
        syn_mean <- dist$estimate[1]
        syn_sd <- dist$estimate[2]
        # df[paste0(method_sd_cut, "_disc")] <- (mean < syn_mean - 2*syn_sd) | (mean > syn_mean + 2*syn_sd)
        df[paste0(method_sd_cut, "_disc")] <- 2*pnorm(-abs((mean - syn_mean)/syn_sd)) < 1 - c_thresh
        df[paste0(method_sd_cut, "_isnull")] <- !df[paste0(method_sd_cut, "_disc")]
        df[paste0(method_sd_cut, "_pval")] <- 2*pnorm(-abs((mean - syn_mean)/syn_sd))
        rank_name <- paste0(method_sd_cut, "_rank")
        rank_fit_df <- df %>% group_by(hgvs) %>% filter(row_number() == 1) %>% ungroup() %>% 
          mutate(!! rank_name := dense_rank(get(paste0(method_sd_cut, "_pval")))) %>% select(hgvs, one_of(rank_name))
        df <- merge(df, rank_fit_df)
        # df[paste0(method_sd_cut, "_rank")] <- rank(2*pnorm(-abs((mean - syn_mean)/syn_sd)))
      }
    }
  }
  return(df)
}

# get top n discoveries
label_top <- function(df, effect_cols, sd_cols, n_disc, rank_by_effect_size=TRUE) {
  for (i in 1:length(effect_cols)) {
    method <- effect_cols[i]
    if (method %in% colnames(df)) {
      # get sig cutoff by ranking by method and taking top n_discovery ranks
      if (grepl("enrich_score", method, fixed=T)) {
        # create ranks by significance then effect size
        uniq_df <- unique(df %>% group_by(hgvs) %>% select(enrich_pval, matches(method)))
        max_effect <- max(abs(uniq_df[[method]]))
        # uniq_df[paste0(method, "_rank")] <- data.table::frank(list(uniq_df$enrich_pval, abs(uniq_df[[method]]) - max_effect), ties.method = "min")
        uniq_df <- uniq_df %>% arrange(uniq_df$enrich_pval, max_effect - abs(uniq_df[[method]]))
        uniq_df[paste0(method, "_rank")] <- 1:nrow(uniq_df)
      } else if (grepl("mu_mean", method, fixed=T)) {
        model <- strsplit(method, "_mu_mean")[[1]]
        uniq_df <- unique(df %>% group_by(hgvs) %>% select(paste0(model, "_mu_lfsr"), matches(method)))
        max_effect <- max(abs(uniq_df[[method]]))
        uniq_df <- uniq_df %>% arrange(paste0(model, "_mu_lfsr"), max_effect - abs(uniq_df[[method]]))
        # uniq_df[paste0(method, "_rank")] <- data.table::frank(list(uniq_df[[paste0(model, "_mu_lfsr")]], abs(uniq_df[[method]]) - max_effect), ties.method = "min")
        uniq_df[paste0(method, "_rank")] <- 1:nrow(uniq_df)
      } else if (method == "ML_effect_mean") {
        mean <- df[[method]]
        sd <- df[["ML_effect_se"]]
        lfsr <- apply(cbind(pnorm(0, mean, sd), 1-pnorm(0, mean, sd)), 1, min)
        hgvs <- df$hgvs
        uniq_df <- unique(data.frame(hgvs, lfsr, mean))
        max_effect <- max(abs(uniq_df$mean))
        # uniq_df[paste0(method, "_rank")] <- data.table::frank(list(uniq_df$lfsr, abs(uniq_df$mean) - max_effect), ties.method = "min")
        uniq_df <- uniq_df %>% arrange(uniq_df$lfsr, max_effect - abs(uniq_df$mean))
        uniq_df[paste0(method, "_rank")] <- 1:nrow(uniq_df)
        # using null sd equiv- rank by effect size
        # uniq_df[paste0(method, "_effect_rank")] <- data.table::frank(list(abs(uniq_df$mean) - max_effect), ties.method = "min")
        uniq_df <- uniq_df %>% arrange(max_effect - abs(uniq_df$mean))
        uniq_df[paste0(method, "_syn_sd_rank")] <- 1:nrow(uniq_df)
        uniq_df <- uniq_df[, !names(uniq_df) %in% c("mean", "lfsr")]
      } else if (method == "weight_effect_mean") {
        mean <- df[[method]]
        sd <- df[["weight_effect_se"]]
        lfsr <- apply(cbind(pnorm(0, mean, sd), 1-pnorm(0, mean, sd)), 1, min)
        hgvs <- df$hgvs
        uniq_df <- unique(data.frame(hgvs, lfsr, mean))
        max_effect <- max(abs(uniq_df$mean))
        uniq_df <- uniq_df %>% arrange(uniq_df$lfsr, max_effect - abs(uniq_df$mean))
        uniq_df[paste0(method, "_rank")] <- 1:nrow(uniq_df)
        uniq_df <- uniq_df %>% arrange(max_effect - abs(uniq_df$mean))
        uniq_df[paste0(method, "_syn_sd_rank")] <- 1:nrow(uniq_df)
        # uniq_df[paste0(method, "_rank")] <- data.table::frank(list(uniq_df$lfsr, abs(uniq_df$mean) - max_effect), ties.method = "min")
        # uniq_df[paste0(method, "_effect_rank")] <- data.table::frank(list(abs(uniq_df$mean) - max_effect), ties.method = "min")
        uniq_df <- uniq_df[, !names(uniq_df) %in% c("mean", "lfsr")]
      } else if (method == "shep_effect_mean") {
        mean <- df[[method]]
        sd <- df[["shep_effect_se"]]
        lfsr <- apply(cbind(pnorm(0, mean, sd), 1-pnorm(0, mean, sd)), 1, min)
        hgvs <- df$hgvs
        uniq_df <- unique(data.frame(hgvs, lfsr, mean))
        max_effect <- max(abs(uniq_df$mean))
        uniq_df <- uniq_df %>% arrange(uniq_df$lfsr, max_effect - abs(uniq_df$mean))
        uniq_df[paste0(method, "_rank")] <- 1:nrow(uniq_df)
        uniq_df <- uniq_df %>% arrange(max_effect - abs(uniq_df$mean))
        uniq_df[paste0(method, "_syn_sd_rank")] <- 1:nrow(uniq_df)
        uniq_df <- uniq_df[, !names(uniq_df) %in% c("mean", "lfsr")]
        # uniq_df[paste0(method, "_rank")] <- ...
        # uniq_df[paste0(method, "_effect_rank")] <- ...

      } else {
        print(paste("method", method, "not found"))
      }
       # label discoveries by rank
      df <- merge(df, uniq_df)
      df[[paste0(method, "_disc")]] <- df[[paste0(method, "_rank")]] <= n_disc
      df[[paste0(method, "_isnull")]] <- df[[paste0(method, "_rank")]] > n_disc
      if (method == "ML_effect_mean" || method == "weight_effect_mean" || method == "shep_effect_mean") {
        df[[paste0(method, "_syn_sd_disc")]] <- df[[paste0(method, "_syn_sd_rank")]] <= n_disc
        df[[paste0(method, "_syn_sd_isnull")]] <- df[[paste0(method, "_syn_sd_rank")]] > n_disc
      }
    }
  }
  return(df)
}


# calculate FDR fxn
calc_FDR <- function(df, effect_cols) {
  method_FDRs <- c()
  for (method in effect_cols) {
    # colnames
    
    isnull_col <- paste0(method, "_isnull")
    print(isnull_col)
    
    syn_df <- df[df$type=="synonymous",]
    effect_df <- df[df$ML_called_effect,]
    disco_df <- df[!df[[isnull_col]],]
    
    # syn_df <- rank_df[rank_df$type=="synonymous",]
    # effect_df <- rank_df[rank_df$ML_called_effect,]
    # disco_df <- rank_df[!rank_df[[isnull_col]],]
    
    # FDR raw classification + threshold classification
    FDR <- sum(disco_df$type=="synonymous", na.rm=T)/nrow(disco_df) # false positives / total discoveries
    FNR <- sum(effect_df[[paste0(method, "_isnull")]], na.rm=T) / nrow(effect_df) # false negatives (classified zero or uncertain) / total true effects
    FPR <- sum(syn_df[[paste0(method, "_disc")]]) / nrow(syn_df) # false positives / total null variants
    
    FDR_vec <- c(FDR, FNR, FPR)
    names(FDR_vec) <- paste0(method, '_', c('FDR', 'FNR', 'FPR'))
    method_FDRs <- c(method_FDRs, FDR_vec)
  }
  return(method_FDRs)
}



# get sim model stats
get_sim_model_stats <- function(yaml, res_dir) {
  sim_model_stats <- c()
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(yaml$params)) {
        param <- yaml$params[i]
        for (param_val in yaml$param_vals[[i]]) {
          for (iter in 0:(yaml$n_iter-1)) {
            for (rescal in c("", "_rescaled")) {
              try({
                sim_name <- paste0("sim_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter, rescal)
                sim_res_dir <- paste0(res_dir, "/sim_results/", sim_name)
                # get fit model stats
                for (model in yaml$models) {
                  stats_file <- paste0(sim_res_dir, "/baseline_", sim_name, "/", model, "_plots/model_stats.yaml")
                  print(stats_file)
                  model_stats <- data.frame(yaml::read_yaml(stats_file))
                  stat_cols <- colnames(model_stats)
                  model_stats <- cbind(dataset, n_pos, param, param_val, iter, model_stats)
                  model_stats$rescaled <- rescal == "_rescaled"
                  sim_model_stats <- rbind(sim_model_stats, model_stats)
                }
              })
            }
          }
        }
      }
    }
  }
  return(sim_model_stats)
}

get_cor_df <- function(yaml, res_dir, effect_quantile=0) {
  cor_df <- c()
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(yaml$params)) {
        param <- yaml$params[i]
        for (param_val in yaml$param_vals[[i]]) {
          for (iter in 0:(yaml$n_iter-1)) {
            for (rescal in c("_rescaled")) {
              try({
                sim_name <- paste0("sim_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter, rescal)
                sim_res_dir <- paste0(res_dir, "/sim_results/", sim_name)
                # get fit model stats
                for (model in yaml$models) {
                  sim_obj_name <- paste0("sim_obj_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter)
                  sim_obj <- readRDS(paste0(res_dir, "/param_data/", sim_obj_name, ".RData"))
                  effect_cutoff <- quantile(abs(sim_obj$effect_mat[sim_obj$effect_mat != 0]), probs=effect_quantile)
                  # get fit df
                  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
                  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
                  fit_df <- readRDS(paste0(sim_res_dir, "/baseline_", sim_name, "/fitted_df.RData"))
                  # fit_df_called <- label_significant(fit_df, 
                  #                                   effect_cols,
                  #                                   sd_cols, 
                  #                                   c_thresh)
                  # make pos - variant - effect - beta - lfsr - disc mat
                  effect_df <- cbind(expand.grid(position = 1:ncol(sim_obj$effect_mat), mutation = 1:nrow(sim_obj$effect_mat)), effect=unlist(c(t(sim_obj$effect_mat))))
                  effect_df_merged <- merge(fit_df, effect_df)
                  # filter to one obs per variant to get FDR on variant scale
                  effect_df_merged <- effect_df_merged %>% group_by(hgvs) %>% filter(row_number() == 1)
                  for (method in effect_cols) {
                    is_rescaled <- rescal == "_rescaled"
                    corr <- cor(effect_df_merged$effect, effect_df_merged[[method]], use="complete.obs")
                    model_corr <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, corr)
                    cor_df <- rbind(cor_df, model_corr)
                  }
                }
              })
            }
          }
        }
      }
    }
  }
  return(cor_df)
}


# get FDR df from yaml
get_FDR_df <- function(yaml, res_dir, c_thresh=0.95, effect_quantile=0) {
  # extract fits
  FDR_df <- c()
  rank_FDR_df <- c()
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(yaml$params)) {
        param <- yaml$params[i]
        for (param_val in yaml$param_vals[[i]]) {
          for (iter in 0:(yaml$n_iter-1)) {
            for (rescal in c("_rescaled")) {
              try({
              sim_name <- paste0("sim_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter, rescal)
              sim_res_dir <- paste0(res_dir, "/sim_results/", sim_name)
              # get fit model stats
              for (model in yaml$models) {
                stats_file <- paste0(sim_res_dir, "/baseline_", sim_name, "/", model, "_plots/model_stats.yaml")
                print(stats_file)
                model_stats <- data.frame(yaml::read_yaml(stats_file))
                stat_cols <- colnames(model_stats)
                model_stats <- cbind(dataset, n_pos, param, param_val, iter, model_stats)
                model_stats$rescaled <- rescal == "_rescaled"

                # get sim obj
                sim_obj_name <- paste0("sim_obj_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter)
                sim_obj <- readRDS(paste0(res_dir, "/param_data/", sim_obj_name, ".RData"))
                effect_cutoff <- quantile(abs(sim_obj$effect_mat[sim_obj$effect_mat != 0]), probs=effect_quantile)
                # get fit df
                effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
                sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
                fit_df <- readRDS(paste0(sim_res_dir, "/baseline_", sim_name, "/fitted_df.RData"))
                fit_df_called <- label_significant(fit_df, 
                                                  effect_cols,
                                                  sd_cols, 
                                                  c_thresh)
                # make pos - variant - effect - beta - lfsr - disc mat
                effect_df <- cbind(expand.grid(position = 1:ncol(sim_obj$effect_mat), mutation = 1:nrow(sim_obj$effect_mat)), effect=unlist(c(t(sim_obj$effect_mat))))
                effect_df_merged <- merge(fit_df_called, effect_df)
                # filter to one obs per variant to get FDR on variant scale
                effect_df_merged <- effect_df_merged %>% group_by(hgvs) %>% filter(row_number() == 1)

                # FDR df -- dataset n pos param iter model FDR FNR div 
                # method <- paste0(model, "_syn_recalibrated_mu_mean")
                for (method in effect_cols) {
                  # FDR--# of zero effect discoveries / # discoveries
                  disco_df <- effect_df_merged[!is.na(effect_df_merged[[paste0(method, "_disc")]]) & effect_df_merged[[paste0(method, "_disc")]],]
                  FDR <- sum(disco_df$effect==0) / nrow(disco_df)
                  if (nrow(disco_df) == 0) {
                    FDR <- 0
                  }
                  # FNR--# of nonzero effect nondiscoveries / # of nonzero effect
                  nonzero_effect_df <- effect_df_merged[abs(effect_df_merged$effect) > effect_cutoff,]
                  FNR <- sum(!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]) & nonzero_effect_df[[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
                  # FPR--# of zero effect discoveries / # of zero effects
                  FPR <- sum(disco_df$effect==0) / sum(effect_df_merged$effect==0)
                  # compute sim metrics
                  is_rescaled <- rescal == "_rescaled"
                  model_FDR <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, FDR, FNR, FPR)
                  FDR_df <- rbind(FDR_df, model_FDR)

                  # add sd effects
                  if (paste0(method, "_syn_sd_disc") %in% colnames(effect_df_merged)) {
                    method <- paste0(method, "_syn_sd")
                     # FDR--# of zero effect discoveries / # discoveries
                    disco_df <- effect_df_merged[!is.na(effect_df_merged[[paste0(method, "_disc")]]) & effect_df_merged[[paste0(method, "_disc")]],]
                    FDR <- sum(disco_df$effect==0) / nrow(disco_df)
                    if (nrow(disco_df) == 0) {
                      FDR <- 0
                    }
                    # FNR--# of nonzero effect nondiscoveries / # of nonzero effect
                    nonzero_effect_df <- effect_df_merged[abs(effect_df_merged$effect) > effect_cutoff,]
                    FNR <- sum(!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]) & nonzero_effect_df[[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
                    # FPR--# of zero effect discoveries / # of zero effects
                    FPR <- sum(disco_df$effect==0) / sum(effect_df_merged$effect==0)
                    # compute sim metrics
                    is_rescaled <- rescal == "_rescaled"
                    model_FDR <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, FDR, FNR, FPR)
                    FDR_df <- rbind(FDR_df, model_FDR)
                  }
                }
              }
              })
            }
          }
        }
      }
    }
  }
  return(FDR_df)
}


# get rank df from yaml + param setting (default: none zero)
get_rank_FDR_df <- function(yaml, res_dir, rank_param, rank_param_val, c_thresh=0.95, effect_cutoff=0) {
  rank_FDR_df <- c()
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(yaml$params)) {
        param <- yaml$params[i]
        for (param_val in yaml$param_vals[[i]]) {
          if (param == rank_param && param_val == rank_param_val) {
            for (iter in 0:(yaml$n_iter-1)) {
              for (rescal in c("", "_rescaled")) {
                try({
                  sim_name <- paste0("sim_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter, rescal)
                  sim_res_dir <- paste0(res_dir, "/sim_results/", sim_name)
                  # get fit model stats
                  for (model in yaml$models) {
                    stats_file <- paste0(sim_res_dir, "/baseline_", sim_name, "/", model, "_plots/model_stats.yaml")
                    print(stats_file)
                    model_stats <- data.frame(yaml::read_yaml(stats_file))
                    stat_cols <- colnames(model_stats)
                    model_stats <- cbind(dataset, n_pos, param, param_val, iter, model_stats)
                    model_stats$rescaled <- rescal == "_rescaled"

                    # get sim obj
                    sim_obj_name <- paste0("sim_obj_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter)
                    sim_obj <- readRDS(paste0(res_dir, "/param_data/", sim_obj_name, ".RData"))
                    # get fit df
                    effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
                    sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
                    fit_df <- readRDS(paste0(sim_res_dir, "/baseline_", sim_name, "/fitted_df.RData"))
                    fit_df_called <- label_significant(fit_df,
                                                      effect_cols,
                                                      sd_cols,
                                                      c_thresh)
                    # make pos - variant - effect - beta - lfsr - disc mat
                    effect_df <- cbind(expand.grid(position = 1:ncol(sim_obj$effect_mat), mutation = 1:nrow(sim_obj$effect_mat)), effect=unlist(c(t(sim_obj$effect_mat))))
                    effect_df_merged <- merge(fit_df_called, effect_df)
                    # filter to one obs per variant to get FDR on variant scale
                    effect_df_merged <- effect_df_merged %>% group_by(hgvs) %>% filter(row_number() == 1)

                    print("ranking")
                    print(iter)
                    for (rank_cutoff in seq(100, max(100, nrow(effect_df_merged)), 100)) {
                      for (method in effect_cols) {
                        # ranks <- rank(effect_df_merged[[method]])
                        # rank_df <- effect_df_merged[ranks < rank_cutoff,]
                        rank_df <- effect_df_merged[effect_df_merged[[paste0(method, "_rank")]] < rank_cutoff,]

                        # disco_df <- rank_df[rank_df[[paste0(method, "_disc")]],]
                        # discoveries is everything less than rank
                        disco_df <- rank_df
                        FDR <- sum(disco_df$effect==0) / nrow(disco_df)
                        # FNR--# of nonzero effect nondiscoveries / # of nonzero effect
                        FNR <- 1 - sum(abs(rank_df$effect) > effect_cutoff) / sum(abs(effect_df_merged$effect) > effect_cutoff)
                        # FNR <- sum(nonzero_effect_df[[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
                        # FPR--# of zero effect discoveries / # of zero effects
                        disco_df <- rank_df[rank_df[[paste0(method, "_disc")]],]
                        FPR <- sum(disco_df$effect==0) / sum(rank_df$effect==0)
                        sens <- sum(disco_df$effect!=0) / sum(rank_df$effect!=0)
                        # compute sim metrics
                        is_rescaled <- rescal == "_rescaled"
                        model_FDR <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, FDR, FNR, FPR, sens, rank_cutoff)
                        rank_FDR_df <- rbind(rank_FDR_df, model_FDR)
                        # add sd effects
                        if (paste0(method, "_syn_sd_disc") %in% colnames(effect_df_merged)) {
                          method <- paste0(method, "_syn_sd")
                          rank_df <- effect_df_merged[effect_df_merged[[paste0(method, "_rank")]] < rank_cutoff,]
                          # FDR--# of zero effect discoveries / # discoveries
                          # disco_df <- rank_df[rank_df[[paste0(method, "_disc")]],]
                          disco_df <- rank_df
                          FDR <- sum(disco_df$effect==0) / nrow(disco_df)
                          # FNR--# of nonzero effect nondiscoveries / # of nonzero effect
                          FNR <- 1 - sum(abs(rank_df$effect) > effect_cutoff) / sum(abs(effect_df_merged$effect) > effect_cutoff)
                          # nonzero_effect_df <- rank_df[abs(rank_df$effect) > effect_cutoff,]
                          # FNR <- sum(nonzero_effect_df[[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
                          # FPR--# of zero effect discoveries / # of zero effects
                          disco_df <- rank_df[rank_df[[paste0(method, "_disc")]],]
                          FPR <- sum(disco_df$effect==0) / sum(rank_df$effect==0)
                          sens <- sum(disco_df$effect!=0) / sum(rank_df$effect!=0)
                          # compute sim metrics
                          is_rescaled <- rescal == "_rescaled"
                          model_FDR <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, FDR, FNR, FPR, sens, rank_cutoff)
                          rank_FDR_df <- rbind(rank_FDR_df, model_FDR)
                        }
                      }
                    }
                  }
                })
              }
            }
          }
        }
      }
    }
  }
  return(rank_FDR_df)
}


# get sig cutoff df from yaml + param setting
get_sig_FDR_df <- function(yaml, res_dir, sig_param, sig_param_val, effect_cutoff=0) {
  FDR_FNR_df <- c()
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(yaml$params)) {
        param <- yaml$params[i]
        for (param_val in yaml$param_vals[[i]]) {
          if (param == sig_param & param_val == sig_param_val) {
            for (iter in 0:(yaml$n_iter-1)) {
              for (rescal in c("", "_rescaled")) {
                try({
                  sim_name <- paste0("sim_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter, rescal)
                  sim_res_dir <- paste0(res_dir, "/sim_results/", sim_name)
                  # get fit model stats
                  for (model in yaml$models) {
                    stats_file <- paste0(sim_res_dir, "/baseline_", sim_name, "/", model, "_plots/model_stats.yaml")
                    print(stats_file)
                    model_stats <- data.frame(yaml::read_yaml(stats_file))
                    stat_cols <- colnames(model_stats)
                    model_stats <- cbind(dataset, n_pos, param, param_val, iter, model_stats)
                    model_stats$rescaled <- rescal == "_rescaled"

                    # get sim obj
                    sim_obj_name <- paste0("sim_obj_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter)
                    sim_obj <- readRDS(paste0(res_dir, "/param_data/", sim_obj_name, ".RData"))
                    # get fit df
                    effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
                    sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
                    fit_df <- readRDS(paste0(sim_res_dir, "/baseline_", sim_name, "/fitted_df.RData"))
                    for (sig_cutoff in seq(0.999, 0.5, length=100)) {
                        fit_df_called <- label_significant(fit_df, 
                                                      effect_cols,
                                                      sd_cols, 
                                                      sig_cutoff)
                      # make pos - variant - effect - beta - lfsr - disc mat
                      effect_df <- cbind(expand.grid(position = 1:ncol(sim_obj$effect_mat), mutation = 1:nrow(sim_obj$effect_mat)), effect=unlist(c(t(sim_obj$effect_mat))))
                      effect_df_merged <- merge(fit_df_called, effect_df)
                      # filter to one obs per variant to get FDR on variant scale
                      effect_df_merged <- effect_df_merged %>% group_by(hgvs) %>% filter(row_number() == 1)

                      for (method in effect_cols) {
                        rank_df <- effect_df_merged
                        
                        disco_df <- rank_df[!is.na(rank_df[[paste0(method, "_disc")]]) & rank_df[[paste0(method, "_disc")]],]
                        FDR <- sum(disco_df$effect==0) / nrow(disco_df)
                        # FNR--# of nonzero effect nondiscoveries / # of nonzero effect
                        nonzero_effect_df <- rank_df[abs(rank_df$effect) > effect_cutoff,]
                        FNR <- sum(!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]) & nonzero_effect_df[[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
                        # FPR--# of zero effect discoveries / # of zero effects
                        FPR <- sum(disco_df$effect==0) / sum(rank_df$effect==0)
                        # sens
                        sens <- sum(abs(disco_df$effect) > effect_cutoff) / sum(abs(rank_df$effect) > effect_cutoff)
                        # compute sim metrics
                        is_rescaled <- rescal == "_rescaled"
                        model_FDR <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, FDR, FNR, FPR, sens, sig_cutoff)
                        FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
                        # add sd effects
                        if (paste0(method, "_syn_sd_disc") %in% colnames(effect_df_merged)) {
                          method <- paste0(method, "_syn_sd")
                          disco_df <- rank_df[!is.na(rank_df[[paste0(method, "_disc")]]) & rank_df[[paste0(method, "_disc")]],]
                          FDR <- sum(disco_df$effect==0) / nrow(disco_df)
                          # FNR--# of nonzero effect nondiscoveries / # of nonzero effect
                          nonzero_effect_df <- rank_df[abs(rank_df$effect) > effect_cutoff,]
                          FNR <- sum(!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]) & nonzero_effect_df[[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
                          # FPR--# of zero effect discoveries / # of zero effects
                          FPR <- sum(disco_df$effect==0) / sum(rank_df$effect==0)
                          # sens
                          sens <- sum(abs(disco_df$effect) > effect_cutoff) / sum(abs(rank_df$effect) > effect_cutoff)
                          # compute sim metrics
                          is_rescaled <- rescal == "_rescaled"
                          model_FDR <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, FDR, FNR, FPR, sens, sig_cutoff)
                          FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
                        }
                      }
                    }
                  }
                })
              }
            }
          }
        }
      }
    }
  }
  return(FDR_FNR_df)
}

# get n discovery df from yaml + param setting
get_n_discovery_FDR_df <- function(yaml, res_dir, disc_param, disc_param_val, n_discoveries, effect_cutoff=0) {
  n_disc_FDR_FNR_df <- c()
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(yaml$params)) {
        param <- yaml$params[i]
        for (param_val in yaml$param_vals[[i]]) {
          if (param == disc_param && param_val == disc_param_val) {
            for (iter in 0:(yaml$n_iter-1)) {
              for (rescal in c("", "_rescaled")) {
                try({
                sim_name <- paste0("sim_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter, rescal)
                sim_res_dir <- paste0(res_dir, "/sim_results/", sim_name)
                # get fit model stats
                for (model in yaml$models) {
                  stats_file <- paste0(sim_res_dir, "/baseline_", sim_name, "/", model, "_plots/model_stats.yaml")
                  print(stats_file)
                  model_stats <- data.frame(yaml::read_yaml(stats_file))
                  stat_cols <- colnames(model_stats)
                  model_stats <- cbind(dataset, n_pos, param, param_val, iter, model_stats)
                  model_stats$rescaled <- rescal == "_rescaled"

                  # get sim obj
                  sim_obj_name <- paste0("sim_obj_", dataset, "_", n_pos, "pos", "_", param, "_", param_val, "_", iter)
                  sim_obj <- readRDS(paste0(res_dir, "/param_data/", sim_obj_name, ".RData"))
                  # get fit df
                  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
                  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
                  fit_df <- readRDS(paste0(sim_res_dir, "/baseline_", sim_name, "/fitted_df.RData"))
                  
                    
                      fit_df_called <- label_top(fit_df, 
                                              effect_cols,
                                              sd_cols, 
                                              n_discoveries)
                      effect_df <- cbind(expand.grid(position = 1:ncol(sim_obj$effect_mat), mutation = 1:nrow(sim_obj$effect_mat)), effect=unlist(c(t(sim_obj$effect_mat))))
                      effect_df_merged <- merge(fit_df_called, effect_df)
                      # filter to one obs per variant to get FDR on variant scale
                      effect_df_merged <- effect_df_merged %>% group_by(hgvs) %>% filter(row_number() == 1)
                      
                      for (method in effect_cols) {
                          
                          # make pos - variant - effect - beta - lfsr - disc mat
                          rank_df <- effect_df_merged
                          
                          disco_df <- rank_df[rank_df[[paste0(method, "_disc")]],]
                          FDR <- sum(disco_df$effect==0) / nrow(disco_df)
                          # FNR--# of nonzero effect nondiscoveries / # of nonzero effect
                          nonzero_effect_df <- rank_df[abs(rank_df$effect) > effect_cutoff,]
                          FNR <- sum(nonzero_effect_df[[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
                          # FPR--# of zero effect discoveries / # of zero effects
                          FPR <- sum(disco_df$effect==0) / sum(rank_df$effect==0)
                          # sens
                          sens <- sum(disco_df$effect!=0) / sum(rank_df$effect!=0)
                          # compute sim metrics
                          is_rescaled <- rescal == "_rescaled"
                          model_FDR <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, FDR, FNR, FPR, sens)
                          n_disc_FDR_FNR_df <- rbind(n_disc_FDR_FNR_df, model_FDR)
                          
                          if (method == "ML_effect_mean" || method == "weight_effect_mean" || method == "shep_effect_mean") {
                            method <- paste0(method, "_syn_sd")
                            disco_df <- rank_df[rank_df[[paste0(method, "_disc")]],]
                            FDR <- sum(disco_df$effect==0) / nrow(disco_df)
                            # FNR--# of nonzero effect nondiscoveries / # of nonzero effect
                            nonzero_effect_df <- rank_df[abs(rank_df$effect) > effect_cutoff,]
                            FNR <- sum(nonzero_effect_df[[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
                            # FPR--# of zero effect discoveries / # of zero effects
                            FPR <- sum(disco_df$effect==0) / sum(rank_df$effect==0)
                            # sens
                            sens <- sum(disco_df$effect!=0) / sum(rank_df$effect!=0)
                            # compute sim metrics
                            is_rescaled <- rescal == "_rescaled"
                            model_FDR <- data.frame(dataset, n_pos, param, param_val, iter, method, is_rescaled, FDR, FNR, FPR, sens)
                            n_disc_FDR_FNR_df <- rbind(n_disc_FDR_FNR_df, model_FDR)
                          }
                      }
                }
                })
              }
            }
          }
        }
      }
    }
  }
  return(n_disc_FDR_FNR_df)
}

get_masked_FDR_df <- function(fit_df, model="FACS_double_sample") {
  FDR_FNR_df <- c()
  # get sim obj
  # get fit df
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  for (sig_cutoff in seq(0.999, 0.8, length=100)) {
    fit_df_called <- label_significant(fit_df, 
                                  effect_cols,
                                  sd_cols, 
                                  sig_cutoff)
    # filter to one obs per variant to get FDR on variant scale
    fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1)
    for (method in effect_cols) {
      rank_df <- fit_df_called
      
      disco_df <- rank_df[!is.na(rank_df[[paste0(method, "_disc")]]) & rank_df[[paste0(method, "_disc")]],]
      # FDR <- sum(disco_df$type=="masked_synonymous") / nrow(disco_df)
      N_syn_called <- sum(disco_df$type=="masked_synonymous") / sum(rank_df$type=="masked_synonymous")
      # compute sim metrics
      model_FDR <- data.frame(method, N_syn_called, sig_cutoff)
      FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      # add sd effects
      if (paste0(method, "_syn_sd_disc") %in% colnames(rank_df)) {
        method <- paste0(method, "_syn_sd")
        disco_df <- rank_df[!is.na(rank_df[[paste0(method, "_disc")]]) & rank_df[[paste0(method, "_disc")]],]
        # FDR <- sum(disco_df$type=="masked_synonymous")  / nrow(disco_df)
        N_syn_called <- sum(disco_df$type=="masked_synonymous") / sum(rank_df$type=="masked_synonymous")
        model_FDR <- data.frame(method, N_syn_called, sig_cutoff)
        FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      }
    }
  }
  return(FDR_FNR_df)
}

get_ranked_masked_FDR_df <- function(fit_df, model="FACS_double_sample") {
  FDR_FNR_df <- c()
  # get sim obj
  # get fit df
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  fit_df_called <- label_significant(fit_df, 
                                  effect_cols,
                                  sd_cols,
                                  c_thresh=0.95)
  # filter to one obs per variant to get FDR on variant scale
  fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1)
  for (rank_cutoff in seq(100, max(100, nrow(fit_df_called)), 100)) {
    for (method in effect_cols) {
      rank_df <- fit_df_called
      disco_df <- rank_df[!is.na(rank_df[[paste0(method, "_rank")]]) & rank_df[[paste0(method, "_rank")]] < rank_cutoff,]
      N_syn_called <- sum(disco_df$type=="masked_synonymous") / sum(rank_df$type=="masked_synonymous")
      # compute sim metrics
      model_FDR <- data.frame(method, N_syn_called, rank_cutoff)
      FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      # add sd effects
      if (paste0(method, "_syn_sd_disc") %in% colnames(rank_df)) {
        method <- paste0(method, "_syn_sd")
        disco_df <- rank_df[!is.na(rank_df[[paste0(method, "_rank")]]) & rank_df[[paste0(method, "_rank")]] < rank_cutoff,]
        N_syn_called <- sum(disco_df$type=="masked_synonymous") / sum(rank_df$type=="masked_synonymous")
        model_FDR <- data.frame(method, N_syn_called, rank_cutoff)
        FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      }
    }
  }
  return(FDR_FNR_df)
}

get_am_FNR_df <- function(fit_df, model="FACS_double_sample", threshold="label") {
  FDR_FNR_df <- c()
  # get sim obj
  # get fit df
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  for (sig_cutoff in seq(0.999, 0.8, length=100)) {
    fit_df_called <- label_significant(fit_df, 
                                  effect_cols,
                                  sd_cols, 
                                  sig_cutoff)
    # filter to one obs per variant to get FDR on variant scale
    fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1)
    for (method in effect_cols) {
      rank_df <- fit_df_called
      
      if (threshold == "label") {
        nonzero_effect_df <- rank_df[!is.na(rank_df$class) & rank_df$class == "pathogenic",]
      } else if (threshold == "10percentile") {
        nonzero_effect_df <- rank_df[!is.na(rank_df$class) & rank_df$val > quantile(rank_df$val, probs=0.9, na.rm=T),]
      } else {
        stop("threshold category not found")
      }
      FNR <- sum(nonzero_effect_df[!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]),][[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
      # compute sim metrics
      model_FDR <- data.frame(method, FNR, sig_cutoff)
      FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      # add sd effects
      if (paste0(method, "_syn_sd_disc") %in% colnames(rank_df)) {
        method <- paste0(method, "_syn_sd")
        if (threshold == "label") {
          nonzero_effect_df <- rank_df[!is.na(rank_df$class) & rank_df$class == "pathogenic",]
        } else if (threshold == "10percentile") {
          nonzero_effect_df <- rank_df[!is.na(rank_df$class) & rank_df$val > quantile(rank_df$val, probs=0.9, na.rm=T),]
        } else {
          stop("threshold category not found")
        }
        FNR <- sum(nonzero_effect_df[!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]),][[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
        # compute sim metrics
        model_FDR <- data.frame(method, FNR, sig_cutoff)
        FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      }
    }
  }
  return(FDR_FNR_df)
}

make_am_box_df <- function(fit_df, am_df) {
  FDR_FNR_df <- c()
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  am_df <- am_df[am_df$am_class == "likely_benign" | am_df$am_class == "likely_pathogenic",]
  sig_cutoff <- 0.95

  fit_df_called <- label_significant(fit_df, 
                                  effect_cols,
                                  sd_cols, 
                                  sig_cutoff)
  # filter to one obs per variant to get FDR on variant scale
  fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1) %>% arrange(position)
  if (protein == "Kir21") {
    # correct for FLAG tag
    am_df$kcnj2_hgvs <- am_df$hgvs
    fit_df_called$kcnj2_pos <- ifelse(fit_df_called$position < 115, fit_df_called$position, -1)
    fit_df_called$kcnj2_pos <- ifelse(fit_df_called$position > 123, fit_df_called$position - 8, fit_df_called$kcnj2_pos)
    fit_df_called$kcnj2_hgvs <- paste0("p.(", fit_df_called$wildtype, fit_df_called$kcnj2_pos, fit_df_called$mutation, ")")
    fit_df_called_am <- merge(fit_df_called, am_df, by=c("kcnj2_hgvs", "wildtype", "mutation"))
    fit_df_called_am$position <- fit_df_called_am$position.x
  } else {
    fit_df_called_am <- merge(fit_df_called, am_df)
  }
  methods <- effect_cols
  for (method in effect_cols) {
    if (paste0(method, "_syn_sd_disc") %in% colnames(fit_df_called_am)) {
        methods <- c(methods, paste0(method, "_syn_sd"))
    }
  }
  for (method in methods) {
    disco_df <- fit_df_called_am[!is.na(fit_df_called_am[[paste0(method, "_disc")]]) & fit_df_called_am[[paste0(method, "_disc")]],]
    nonzero_effect_df <- fit_df_called_am[fit_df_called_am$am_class=="likely_pathogenic",]

    FDR <- sum(disco_df$am_class=="likely_benign") / nrow(disco_df)
    FNR <- sum(nonzero_effect_df[!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]),][[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
    sens <- 1-FNR
    model_FDR <- data.frame(method, FDR, FNR, sens)
    FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
  }
  return(FDR_FNR_df)
}

make_clinvar_box_df <- function(fit_df, clinvar_df) {
  FDR_FNR_df <- c()
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  sig_cutoff <- 0.95

  fit_df_called <- label_significant(fit_df, 
                                  effect_cols,
                                  sd_cols, 
                                  sig_cutoff)
  # filter to one obs per variant to get FDR on variant scale
  fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1) %>% arrange(position)
  if (protein == "Kir21") {
    # correct for FLAG tag
    clinvar_df$kcnj2_hgvs <- clinvar_df$hgvs
    fit_df_called$kcnj2_pos <- ifelse(fit_df_called$position < 115, fit_df_called$position, -1)
    fit_df_called$kcnj2_pos <- ifelse(fit_df_called$position > 123, fit_df_called$position - 8, fit_df_called$kcnj2_pos)
    fit_df_called$kcnj2_hgvs <- paste0("p.(", fit_df_called$wildtype, fit_df_called$kcnj2_pos, fit_df_called$mutation, ")")
    fit_df_called_clinvar <- merge(fit_df_called, clinvar_df, by=c("kcnj2_hgvs", "wildtype", "mutation"))
    fit_df_called_clinvar$position <- fit_df_called_clinvar$position.x
  } else {
    fit_df_called_clinvar <- merge(fit_df_called, clinvar_df)
  }
  methods <- effect_cols
  for (method in effect_cols) {
    if (paste0(method, "_syn_sd_disc") %in% colnames(fit_df_called_clinvar)) {
        methods <- c(methods, paste0(method, "_syn_sd"))
    }
  }
  for (method in methods) {
    disco_df <- fit_df_called_clinvar[!is.na(fit_df_called_clinvar[[paste0(method, "_disc")]]) & fit_df_called_clinvar[[paste0(method, "_disc")]],]
    nonzero_effect_df <- fit_df_called_clinvar[fit_df_called_clinvar$class=="pathogenic",]

    FDR <- sum(disco_df$class=="benign") / nrow(disco_df)
    FNR <- sum(nonzero_effect_df[!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]),][[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
    sens <- 1-FNR
    model_FDR <- data.frame(method, FDR, FNR, sens)
    FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
  }
  return(FDR_FNR_df)

  # include num benign and num pathogenic
}

make_am_curve_df <- function(fit_df, am_df, model) {
  FDR_FNR_df <- c()
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  am_df <- am_df[am_df$am_class == "likely_benign" | am_df$am_class == "likely_pathogenic",]
  sig_cutoff <- 0.95

  fit_df_called <- label_significant(fit_df, 
                                  effect_cols,
                                  sd_cols, 
                                  sig_cutoff)
  # filter to one obs per variant to get FDR on variant scale
  fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1) %>% arrange(position)
  if (protein == "Kir21") {
    # correct for FLAG tag
    am_df$kcnj2_hgvs <- am_df$hgvs
    fit_df_called$kcnj2_pos <- ifelse(fit_df_called$position < 115, fit_df_called$position, -1)
    fit_df_called$kcnj2_pos <- ifelse(fit_df_called$position > 123, fit_df_called$position - 8, fit_df_called$kcnj2_pos)
    fit_df_called$kcnj2_hgvs <- paste0("p.(", fit_df_called$wildtype, fit_df_called$kcnj2_pos, fit_df_called$mutation, ")")
    fit_df_called_am <- merge(fit_df_called, am_df, by=c("kcnj2_hgvs", "wildtype", "mutation"))
    fit_df_called_am$position <- fit_df_called_am$position.x
  } else {
    fit_df_called_am <- merge(fit_df_called, am_df)
  }
  methods <- effect_cols
  for (method in effect_cols) {
    if (paste0(method, "_syn_sd_disc") %in% colnames(fit_df_called_am)) {
        methods <- c(methods, paste0(method, "_syn_sd"))
    }
  }
  for (method in methods) {
    TP <- 0
    FP <- 0
    TN <- 0
    FN <- 0
    rank_df <- fit_df_called_am %>% arrange(get(paste0(method, "_rank")))
    for (i in 1:nrow(rank_df)) {
      obs <- rank_df[i,]
      if (obs$am_class == "likely_benign") {
        FP <- FP + 1
      } else if (obs$am_class == "likely_pathogenic") {
        TP <- TP + 1
      } else {
        print("not benign or pathogenic")
      }
      # false positives / total discoveries
      FDR <- FP / (FP + TP)
      # false negatives / total negatives
      FNR <- FN / sum(fit_df_called_am$am_class == "likely_benign")
      # true positive / total positives
      TPR <- TP / sum(fit_df_called_am$am_class == "likely_pathogenic")
      # false positive / total negatives
      FPR <- FP / sum(fit_df_called_am$am_class == "likely_benign")
      model_FDR <- data.frame(method, FDR, FNR, TPR, FPR, i)
      FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      
    }
  }
  return(FDR_FNR_df)

} 

get_nonsense_FNR_df <- function(fit_df, model="FACS_double_sample_repq") {
  FDR_FNR_df <- c()
  # get sim obj
  # get fit df
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  sig_cutoff <- 0.95
  fit_df_called <- label_significant(fit_df, 
                                effect_cols,
                                sd_cols, 
                                sig_cutoff)
  # filter to one obs per variant to get FDR on variant scale
  fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1)
  for (method in effect_cols) {
    rank_df <- fit_df_called
    nonzero_effect_df <- rank_df[rank_df$type == "nonsense",]
    FNR <- sum(nonzero_effect_df[!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]),][[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
    sens <- 1-FNR
    # compute sim metrics
    model_FDR <- data.frame(method, FNR, sens, sig_cutoff)
    FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
    # add sd effects
    if (paste0(method, "_syn_sd_disc") %in% colnames(rank_df)) {
      method <- paste0(method, "_syn_sd")
      FNR <- sum(nonzero_effect_df[!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]),][[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
      sens <- 1-FNR
      # compute sim metrics
      model_FDR <- data.frame(method, FNR, sens, sig_cutoff)
      FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
    }
  }
  # }
  return(FDR_FNR_df)
}

get_ranked_nonsense_FNR_df <- function(fit_df, model="FACS_double_sample") {
  FDR_FNR_df <- c()
  # get sim obj
  # get fit df
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  fit_df_called <- label_significant(fit_df, 
                                  effect_cols,
                                  sd_cols, 
                                  sig_cutoff)
  # filter to one obs per variant to get FDR on variant scale
  fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1)
  for (rank_cutoff in seq(100, max(100, nrow(fit_df_called)), 100)) {
    for (method in effect_cols) {
      rank_df <- fit_df_called
      disco_df <- rank_df[!is.na(rank_df[[paste0(method, "_rank")]]) & rank_df[[paste0(method, "_rank")]] < rank_cutoff,]
      N_nonsense_called <- sum(disco_df$type=="nonsense") / sum(rank_df$type=="nonsense")
      # compute sim metrics
      model_FDR <- data.frame(method, N_nonsense_called, rank_cutoff)
      FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      # add sd effects
      if (paste0(method, "_syn_sd_disc") %in% colnames(rank_df)) {
        method <- paste0(method, "_syn_sd")
        disco_df <- rank_df[!is.na(rank_df[[paste0(method, "_rank")]]) & rank_df[[paste0(method, "_rank")]] < rank_cutoff,]
        N_nonsense_called <- sum(disco_df$type=="nonsense") / sum(rank_df$type=="nonsense")
        model_FDR <- data.frame(method, N_nonsense_called, rank_cutoff)
        FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      }
    }
  }
  return(FDR_FNR_df)
}

get_clinvar_FNR_df <- function(fit_df, model="FACS_double_sample") {
  FDR_FNR_df <- c()
  # get sim obj
  # get fit df
  effect_cols <- c(paste0(model, "_mu_mean"), paste0(model, "_syn_recalibrated_mu_mean"), "enrich_score", "weight_effect_mean", "shep_effect_mean", "ML_effect_mean")
  sd_cols <- c(NA, NA, "enrich_SE", "weight_effect_se", "shep_effect_se", "ML_effect_se")
  for (sig_cutoff in seq(0.999, 0.8, length=100)) {
    fit_df_called <- label_significant(fit_df, 
                                  effect_cols,
                                  sd_cols, 
                                  sig_cutoff)
    # filter to one obs per variant to get FDR on variant scale
    fit_df_called <- fit_df_called %>% group_by(hgvs) %>% filter(row_number() == 1)
    for (method in effect_cols) {
      rank_df <- fit_df_called
      nonzero_effect_df <- rank_df[!is.na(rank_df$`Germline classification`) & 
        (rank_df$`Germline classification` == "Pathogenic" | rank_df$`Germline classification` == "Pathogenic/Likely pathogenic" | rank_df$`Germline classification` == "Likely pathogenic"),]
      FNR <- sum(nonzero_effect_df[!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]),][[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
      # compute sim metrics
      model_FDR <- data.frame(method, FNR, sig_cutoff)
      FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      # add sd effects
      if (paste0(method, "_syn_sd_disc") %in% colnames(rank_df)) {
        method <- paste0(method, "_syn_sd")
        FNR <- sum(nonzero_effect_df[!is.na(nonzero_effect_df[[paste0(method, "_isnull")]]),][[paste0(method, "_isnull")]]) / nrow(nonzero_effect_df)
        # compute sim metrics
        model_FDR <- data.frame(method, FNR, sig_cutoff)
        FDR_FNR_df <- rbind(FDR_FNR_df, model_FDR)
      }
    }
  }
  return(FDR_FNR_df)
}


# rename methods
rename_methods <- function(FDR_FNR_df) {
  FDR_FNR_df <- FDR_FNR_df[FDR_FNR_df$method != "shep_effect_mean" &
    FDR_FNR_df$method != "shep_effect_mean_syn_sd" &  
    FDR_FNR_df$method != "FACS_phi_widepos_syndist_nfunc_rep_optim_mu_mean" & 
    FDR_FNR_df$method != "FACS_phi_widepos_syndist_nfunc_rep_optim_syn_recalibrated_mu_mean",]
  FDR_FNR_df[FDR_FNR_df$method=="FACS_double_sample_repq_mu_mean",]$method <- "Lilace\n(unrecalibrated)"
  FDR_FNR_df[FDR_FNR_df$method=="FACS_double_sample_repq_syn_recalibrated_mu_mean",]$method <- "Lilace"
  FDR_FNR_df[FDR_FNR_df$method=="enrich_score",]$method <- "enrich2\np-value"
  FDR_FNR_df[FDR_FNR_df$method=="enrich_score_syn_sd",]$method <- "enrich2\n2sd"
  FDR_FNR_df[FDR_FNR_df$method=="weight_effect_mean",]$method <- "mean bin\nnorm approx"
  FDR_FNR_df[FDR_FNR_df$method=="weight_effect_mean_syn_sd",]$method <- "mean bin\n2sd"
  FDR_FNR_df[FDR_FNR_df$method=="ML_effect_mean_syn_sd",]$method <- "ML\n2sd"
  FDR_FNR_df[FDR_FNR_df$method=="ML_effect_mean",]$method <- "ML\nnorm approx"
  return(FDR_FNR_df)
}

rename_params <- function(FDR_FNR_df) {
  FDR_FNR_df[FDR_FNR_df$param=="none",]$param <- "default"
  # FDR_FNR_df[FDR_FNR_df$param=="n_variants",]$param <- "Number of\nVariants"
  FDR_FNR_df[FDR_FNR_df$param=="effect_prop",]$param <- "Proportion of\nEffects"
  FDR_FNR_df[FDR_FNR_df$param=="measurement_variance",]$param <- "Fluorescence\nVariance"
  # FDR_FNR_df[FDR_FNR_df$param=="position_effect",]$param <- "Position SD"
  FDR_FNR_df[FDR_FNR_df$param=="cell_count_factor",]$param <- "Cell Count\nper\nObservation"
  FDR_FNR_df[FDR_FNR_df$param=="variant_phi",]$param <- "Variant\nOverdispersion"
  # FDR_FNR_df[FDR_FNR_df$param=="n_reads",]$param <- "Read Count"
  # FDR_FNR_df[FDR_FNR_df$param=="reads_phi",]$param <- "Read\nOverdispersion"
  FDR_FNR_df[FDR_FNR_df$param=="reads_per_obs",]$param <- "Reads per\nObservation"
  FDR_FNR_df[FDR_FNR_df$param=="bin_phi",]$param <- "Bin\nOverdispersion"
  FDR_FNR_df[FDR_FNR_df$param=="n_bins",]$param <- "Number of\nBins"
  FDR_FNR_df[FDR_FNR_df$param=="gate_offset",]$param <- "Gate Offset"
  FDR_FNR_df[FDR_FNR_df$param=="prop_small",]$param <- "Proportion\nSmall Effects"
  FDR_FNR_df[FDR_FNR_df$param=="n_replicates",]$param <- "Number of\nReplicates"
  FDR_FNR_df[FDR_FNR_df$param=="rep_effect",]$param <- "Replicate\nEffect"
  return(FDR_FNR_df)
  }