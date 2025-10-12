library(dplyr)
source('../src/build_sim.R')
source('../src/fit_func.R')
source('../src/methods/method_utils.R')
source('../src/plot_func.R')


# cmd line args
args <- commandArgs(trailingOnly = TRUE)
df <- readRDS(args[1]) %>% ungroup()
desc <- yaml::read_yaml(args[2])
methods <- strsplit(args[3], ",")[[1]]
res_dir <- args[4]
outfile <- args[5]
shell <- args[6]
rcfile <- args[7]
# df <- data

K <- sum(startsWith(colnames(df), "c_"))

# filter zero rows
print(paste("Removing", nrow(df) - nrow(df[df$n_counts > 0,]), "zero rows"))
df <- df[df$n_counts > 0,]


if ("enrich2" %in% methods) {
    print("Running enrich2...")
    start <- Sys.time()

    # dataset <- tools::file_path_sans_ext(basename(args[1]))
    # if (dataset == "tpmt_dropNA" || dataset == "pten_dropNA") {
    #     df$rep <- paste0("R", df$rep)
    # }

    df <- df %>% group_by(hgvs) %>% mutate(rep_list=toString(sort(rep)))
    df$exp <- as.character(df$exp)
    df[,c("enrich_pval", "enrich_z", "enrich_SE", "enrich_epsilon", "enrich_score")] <- as.numeric(NA)
    for (rep_group in unique(df$rep_list)) { # or n rep group as R1, R2, R3 combos
        print(rep_group)
        enrich_dir <- paste0(res_dir, "/enrich2_", rep_group)
        enrich_dir <- gsub(" ", "", enrich_dir)
        if (!dir.exists(enrich_dir)) {
            dir.create(enrich_dir)
        }
        enrich_config <- paste0(enrich_dir, '/config.json')
        condition_name <- paste(unlist(unique(df$exp)), collapse='_')
        run_df <- df[df$rep_list==rep_group,]
        write_enrich2_files(run_df, enrich_dir, "synonymous", condition_name, K=K)
        if (K == 2) {
            run_enrich2_ratio(shell, rcfile, enrich_config)
        } else {
            run_enrich2(shell, rcfile, enrich_config)
        }
        enrich2_scores <- extract_enrich_exp(enrich_dir)
        print(rep_group)
        print(head(enrich2_scores))
        # enrich2_scores$exp <- as.numeric(enrich2_scores$exp)
        df <- rows_patch(df, enrich2_scores, by = "hgvs_exp", unmatched="ignore")
    }
    df <- df %>% ungroup()
    # enrich_dir <- paste0(res_dir, "/enrich2")
    # if (!dir.exists(enrich_dir)) {
    #     dir.create(enrich_dir)
    # }
    # enrich_config <- paste0(enrich_dir, '/config.json')
    # condition_name <- paste(unlist(unique(data$exp)), collapse='_')
    # write_enrich2_files(data, enrich_dir, "synonymous", condition_name, K=K)
    # if (K == 2) {
    #     run_enrich2_ratio(shell, rcfile, enrich_config)
    # } else {
    #     run_enrich2(shell, rcfile, enrich_config)
    # }

    # enrich2_scores <- extract_enrich_exp(enrich_dir)
    # df <- merge(data, enrich2_scores, all.x = TRUE)
    end <- Sys.time()
    print(end-start)
}

if ("ML" %in% methods) {
    print("Running ML...")
    start <- Sys.time()
    breaks <- desc$breaks
    vars <- unique(df$hgvs_exp)
    df$ML_mean <- NA
    df$ML_sd <- NA

    df$syn_ML_mean <- NA
    df$syn_ML_sd <- NA
    for (rep in unique(df$rep)) {
        syn_counts <- df[df$type=="synonymous" & df$rep==rep,] %>% ungroup() %>% select(starts_with("c_"))
        print(head(syn_counts))
        opt_out <- optim(c(1000, 6.5), function(theta) cdf_lik_log_multi(theta, as.matrix(syn_counts), breaks, cdf=plnorm))
        if (opt_out$convergence == 0) {
            params <- opt_out$par
        } else {
            params <- c(NA, NA)
        }
        df[df$type=="synonymous" & df$rep==rep,]$syn_ML_mean <- params[1]
        df[df$type=="synonymous" & df$rep==rep,]$syn_ML_sd <- exp(params[2])
    }

    for (i in 1:length(vars)) {
        var <- vars[i]
        var_df <- df[df$hgvs_exp==var,]
        counts <- as.matrix(var_df %>% dplyr::select(starts_with("c_")))
        params <- c(NA, NA)
        opt_out <- optim(c(1000, 6.5), function(theta) cdf_lik_log_multi(theta, counts, breaks, cdf=plnorm))
        if (opt_out$convergence == 0) {
            params <- opt_out$par
        } 
        df[df$hgvs_exp==var,]$ML_mean <- params[1]
        df[df$hgvs_exp==var,]$ML_sd <- exp(params[2])
    }
    df$ML_mean_indiv <- NA
    df$ML_sd_indiv <- NA
    for (i in 1:nrow(df)) {
        counts <- as.numeric(df[i,] %>% ungroup() %>% select(starts_with("c_")))
        params <- c(NA, NA)
        opt_out <- optim(c(1000, 6.5), function(theta) cdf_lik_log(theta, counts, breaks, cdf=plnorm))
        if (opt_out$convergence == 0) {
            params <- opt_out$par
        } 
        df[i,]$ML_mean_indiv <- params[1]
        df[i,]$ML_sd_indiv <- exp(params[2])
    }
    # get likelihoods
    df$ML_syn_lik <- NA
    df$ML_param_lik <- NA
    for (rep in unique(df$rep)) {
        # get syn ML params for each replicate
        syn_ML_mean <- unique(df[df$type=="synonymous" & df$rep == rep,]$syn_ML_mean)
        syn_ML_sd <- unique(df[df$type=="synonymous" & df$rep == rep,]$syn_ML_sd)
        for (i in 1:nrow(df[df$rep==rep,])) {
            counts <- as.numeric(df[df$rep==rep,][i,] %>% ungroup() %>% select(starts_with("c_")))
            param_ML_mean <- df[df$rep==rep,][i,]$ML_mean_indiv
            param_ML_sd <- df[df$rep==rep,][i,]$ML_sd_indiv
            lik_syn <- cdf_lik_log(c(syn_ML_mean, log(syn_ML_sd)), counts, breaks, cdf=plnorm)
            lik_param <- cdf_lik_log(c(param_ML_mean, log(param_ML_sd)), counts, breaks, cdf=plnorm)
            df[df$rep==rep,][i,"ML_syn_lik"] <- lik_syn
            df[df$rep==rep,][i,"ML_param_lik"] <- lik_param
        }
    }
    end <- Sys.time()
    print(end-start)
}

if ("MoM" %in% methods) {
    print("Running MoM...")
    start <- Sys.time()
    df[c("shep_mean", "shep_var")] <- t(apply(df %>% dplyr::select(starts_with("c_")), 1, function(x) {shep_moments(x, desc$breaks)}))
    end <- Sys.time()
    print(end-start)
}

# get effects from ML,MoM,weight distribution estimation
print("Calculating effects...")
start <- Sys.time()
vars <- unique(df$hgvs_exp)
for (i in 1:length(vars)) {
    var <- vars[i]
    var_df <- df[df$hgvs_exp==var,]
    counts <- var_df %>% dplyr::select(starts_with("c_"))

    data_control <- df %>% filter(type == "synonymous")
    control_counts <- data_control %>% dplyr::select(starts_with("c_"))
    mean_control_counts <- colMeans(control_counts)

    if ("ML" %in% methods) {
        control_ML_mean <- mean(data_control$ML_mean, na.rm=T)
        # control_ML_mean_var <- var(data_control$ML_mean, na.rm=T)
        control_ML_var <- mean(data_control$ML_sd, na.rm=T)^2
        ML_effects <- var_df$ML_mean - control_ML_mean
        ML_effect_mean <- mean(ML_effects, na.rm=T)
        ML_effect_se <- sd(ML_effects, na.rm=T)
        ML_effect_se_analytic <- mean(sqrt(var_df$ML_sd^2 / var_df$n_counts + mean(data_control$ML_sd^2 / data_control$n_counts))) # biased estimator of sd (underestimate), also not taking into account variance in control var estimate
        df[df$hgvs_exp==var,"ML_effects"] <- ML_effects
        df[df$hgvs_exp==var,c("ML_effect_mean","ML_effect_se", "ML_effect_se_analytic")] <- matrix(c(ML_effect_mean, ML_effect_se, ML_effect_se_analytic), nrow=nrow(var_df), ncol=3, byrow=T)
        # get p value by combining 3 z-tests w/ fisher's method
        # m1 <- exp(control_ML_mean + control_ML_var / 2)
        # m2 <-  exp(var_df$ML_mean + var_df$ML_sd^2 / 2)
        # v1 <- (exp(control_ML_var)-1) * exp(2*control_ML_mean + control_ML_var)
        # v2 <- (exp(var_df$ML_sd^2)-1) * exp(2*var_df$ML_mean + var_df$ML_sd^2)
        #  z <- (m2 - m1) / sqrt(v2 + v1)
        z <- (var_df$ML_mean - control_ML_mean) / sqrt(var_df$ML_sd^2 + control_ML_var)
        p <- 2*pnorm(-abs(z))
        df[df$hgvs_exp==var,"ML_obs_pval"] <- p
        df[df$hgvs_exp==var,"ML_pval"] <- pchisq(-2 * sum(log(p)), 2*length(p), lower.tail=F)
    }

    if ("MoM" %in% methods) {
        control_shep_mean <- mean(data_control$shep_mean, na.rm=T)
        # control_shep_mean_var <- var(data_control$shep_mean, na.rm=T)
        control_shep_var <- mean(data_control$shep_var, na.rm=T)
        shep_effects <- var_df$shep_mean - control_shep_mean
        shep_effect_mean <- mean(shep_effects, na.rm=T)
        shep_effect_se <- sd(shep_effects, na.rm=T)
        df[df$hgvs_exp==var,"shep_effects"] <- shep_effects
        df[df$hgvs_exp==var,c("shep_effect_mean","shep_effect_se")] <- matrix(c(shep_effect_mean, shep_effect_se), nrow=nrow(var_df), ncol=2, byrow=T)
        # get p value by combining 3 z-tests w/ fisher's method
        z <- (var_df$shep_mean - control_shep_mean) / sqrt(var_df$shep_var + control_shep_var)
        p <- 2*pnorm(-abs(z))
        df[df$hgvs_exp==var,"shep_obs_pval"] <- p
        df[df$hgvs_exp==var,"shep_pval"] <- pchisq(-2 * sum(log(p)), 2*length(p), lower.tail=F)
    }

    if ("weight_mean" %in% methods) {
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
        df[df$hgvs_exp==var,"weight_effects"] <- weight_effects
        df[df$hgvs_exp==var,c("weight_effect_mean","weight_effect_se")] <- matrix(c(weight_effect_mean, weight_effect_se), nrow=nrow(var_df), ncol=2, byrow=T)
    }
}

saveRDS(df, outfile)
end <- Sys.time()
print(end-start)