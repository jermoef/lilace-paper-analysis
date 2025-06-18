# Functions to construct and generate FACS simulated data
# create_sim -- read data parameters into object
# sample_variant -- given index in variant-pos-condition mat, 
#                  draw variant fluorescence distribution
# sort_samples -- given fluorescence distribution, sort into bins
# add_FACS_df -- simulate whole sorting + sequencing process (calls sort_samples)
# sim_FACS -- set FACS cutoffs based on given strategy and (calls add_FACS_df)
# run_sim -- draws samples (calls sample_variant + sim_FACS)


library(dplyr)


#' Specify simulation parameters
#' @param n_pos number of positions in protein
#' @param n_variants number of variants at each position
#' @param n_conditions number of conditions across experiments
#' @param n_replicates number of replicates for each variant
#' @param ref_loc location parameter of reference (wildtype) phenotype
#' @param ref_scale scale parameter of reference phenotype
#' @param pheno_rfunc_mat phenotype distribution sampling function for each variant
#' @param effect_mat effect size at each variant, position, condition
#' @param n_cells total number of cells for experiment
#' @param measurement_bias_avg avg bias in fluorescence distr
#' @param measurement_bias_sd sd of bias in fluorescence distr
#' @param measurement_noise variance in fluorescence tagging
#' @param sample_phi overdispersion parameter for DirMult sampling of cell counts 
#' @param variant_phi variant level overdspn param
#' @param reads_phi overdispersion parameter for DirMult sampling of reads
#' @param PCR_factor PCR amplification factor
#' @param FACS_strat FACS gating strategy: 'overall', 'log_overall', 'modes'
#' @param error_rate FACS error rate
#' @param n_bins number of FACS bins
#' @param n_reads number of sequencing reads
#' @param equal_samples all variants have same number of cells (o.w. draw sample sizes from multinomial)
#' @param variant_phi_samples a vector of overdispersion (phi) values
create_sim <- function(n_pos, n_variants, n_replicates, cell_counts, cell_count_factor, 
                        effect_mat, syn_mean, measurement_variance, variant_phi, variant_phi_samples, a, b,
                        n_reads, reads_phi, bin_phi, bin_reads_phi_means, n_counts_phi_means, reads_per_obs, 
                        mean_crossrep_cor, FACS_gating, gate_offset, n_bins, rep_effects) {
    # list w/ sim parameters
    sim_obj <- as.list(environment())
    return(sim_obj)
}

#' Given sim parameters and a specific variant, sample full fluorescence distribution
sample_variant <- function(sim_obj, i, j, r) {
    list2env(sim_obj, envir = environment())
    if (sample_counts[i, j] == 0) {
        return(NA)
    }
    return(rlnorm(sample_counts[i, j], syn_mean + effect_mat[i, j] + rep_effects[r], measurement_variance))
}

#' Given fluorescence distribution for variant, sort + bin
sort_samples <- function(sampled_sim_obj, i, j, r) {
    list2env(sampled_sim_obj, envir = environment())
    variant_samples <- samples[!is.na(samples[,1]) & samples[,2]==i & samples[,3]==j & samples[,4]==r, 1]
    levels <- cut(variant_samples, cutoffs)
    print(table(levels, useNA="ifany"))
    bin_counts <- as.numeric(table(levels))
    # assert bin counts matches number of bins
    stopifnot(length(bin_counts)==n_bins)
    return(bin_counts)
}

#' sim FACS helper
add_FACS_df <- function(sampled_sim_obj) {
    list2env(sampled_sim_obj, envir = environment())
    FACS_df <- data.frame()

    # create true bin counts (before sequencing)
    print(paste("Sorting samples into", n_bins, "bins..."))
    for (i in 1:n_variants) {
        for (j in 1:n_pos) {
            for (r in 1:n_replicates) {
                # sort samples
                bin_counts <- sort_samples(sampled_sim_obj, i, j, r)
                type <- ifelse(i==1, "synonymous",
                            ifelse(i <= 20, "missense", "indel"))
                hgvs <- paste0("p.(", 1, 'p', j, 'm', i, ")")
                FACS_df <- rbind(FACS_df, c(hgvs, j, 1, i, type, bin_counts, r))
            }
        }
    }
    colnames(FACS_df) <- c("hgvs", "position", "wildtype", "mutation", "type", 
                        paste0("c_", 0:(n_bins-1)), "rep")
    FACS_df <- FACS_df %>% mutate_at(vars(starts_with("c_")), as.numeric)
    sampled_sim_obj$FACS_df <- FACS_df

    obs_df <- FACS_df

    # sample read count vector
    n_obs_per_replicate <- n_variants * n_pos
    read_counts <- matrix(NA, n_obs_per_replicate, n_replicates)
    for (r in 1:n_replicates) {
        if (reads_phi == -1) {
            est_reads_phi <- n_counts_phi_means[n_counts_phi_means$rep==r,]$mean_eta
        } else {
            est_reads_phi <- reads_phi
        }
        print(paste("est read phi:", est_reads_phi))
        sampled_sim_obj$est_reads_phi <- est_reads_phi
        p_read_count <- brms::rdirichlet(1, est_reads_phi*rep(1/n_obs_per_replicate, n_obs_per_replicate))
        read_counts[,r] <- rmultinom(1, n_reads, p_read_count)
    }
    print(head(read_counts))

    # sample counts from multinomial based on true proportions
    count_df <- obs_df %>% dplyr::select(starts_with("c_"))
    print("Sampling variants...")
    count_props <- sweep(count_df, 1, rowSums(count_df), '/')
    data_count_indices <- sample.int(nrow(data), nrow(count_df), replace=T)
    observed_bin_counts <- t(sapply(1:nrow(count_df), function(i) {
        props <- as.numeric(count_props[i,])
        # cell_count <- sum(count_df[i,])
        rep <- obs_df[i,]$rep
        read_count_index <- sum(obs_df[1:i,]$rep==rep) # number of times this rep covered before and including current index
        read_count <- read_counts[read_count_index, as.numeric(rep)]
        if (reads_phi == -1) {
            read_count <- data$n_counts[data_count_indices[i]]
        }
        if (!is.na(reads_per_obs)) { # if set reads per obs use that instead of sampling
            read_count <- reads_per_obs
        }
        if (bin_phi != -1) {
            read_count <- cell_count_factor
        }
        if (read_count == 0) {
            return(rep(0, n_bins))
        }
        props[props==0] <- props[props==0] + 1e-16 # avoid zero alpha
        if (variant_phi == -2) {
            variant_phi_sample <- a*read_count + b
            p <- brms::rdirichlet(1, variant_phi_sample*props)
            return(rmultinom(1, read_count, p))
        } else if (variant_phi != -1) { # phi value supplied -> use it
            p <- brms::rdirichlet(1, variant_phi*props)
            return(rmultinom(1, read_count, p))
        } else { # randomly draw a phi from real data variant phi vector
            variant_phi_sample <- variant_phi_samples[sample.int(length(variant_phi_samples), 1)]
            if (variant_phi_sample < 0.01) {
                variant_phi_sample <- 0.01 # if too small, rdirichlet results in NAs
            }
            p <- brms::rdirichlet(1, variant_phi_sample*props)
            return(rmultinom(1, read_count, p))
        }
    }))
    obs_df[, paste0("c_", 0:(n_bins-1))] <- as.data.frame(observed_bin_counts)
    sampled_sim_obj$variant_df <- obs_df
    print(head(obs_df))
    if (bin_phi != -1) {
        print("Sampling bin wise reads")
        for (r in 1:n_replicates) {
            count_props <- sweep(count_df[obs_df$rep == r,], 2, colSums(count_df[obs_df$rep == r,]), '/') # note: using count_df here
            count_props[count_props==0] <- count_props[count_props==0] + 1e-16 # avoid zero alpha
            observed_bin_counts <- sapply(1:n_bins, function(k) {
                p <- brms::rdirichlet(1, bin_phi*count_props[,k])
                return(round(rmultinom(1, n_reads/n_bins, p)))
            })
            obs_df[obs_df$rep==r, paste0("c_", 0:(n_bins-1))] <- as.data.frame(observed_bin_counts)
        }
    }
    print(head(obs_df))
    obs_df <- obs_df %>% mutate_at(vars(starts_with("c_")), as.numeric)
    sampled_sim_obj$obs_df <- obs_df

    return(sampled_sim_obj)
}


#' Given sim_obj that already has phenotype samples -> sim FACS experiment to get observed df
sim_FACS <- function(sampled_sim_obj) {
    list2env(sampled_sim_obj, envir = environment())
    if (FACS_gating == "log_equal") {
        qtiles <- quantile(log(samples[,1]), na.rm=T, probs=seq(0,1,1/n_bins))
        qtiles <- qtiles - gate_offset
        cutoffs <- c(0,exp(qtiles[2:(length(qtiles)-1)]), max(samples[,1]))
    } else if (FACS_gating == "equal") {
        qtiles <- quantile(samples[,1], na.rm=T)
        qtiles <- qtiles - gate_offset
        cutoffs <- c(-Inf,qtiles[2:(length(qtiles)-1)], max(samples[,1]))
    } else if (FACS_gating == "modes") {
        # TODO: set gates based on modes of bimodal overall distribution
        print("modes strat not implemented")
        stop()
    } else {
        print("FACS gating strategy not recognized")
        stop()
    }
    sampled_sim_obj$cutoffs <- cutoffs
    sampled_sim_obj <- add_FACS_df(sampled_sim_obj)
    return(sampled_sim_obj)
}

run_sim <- function(sim_obj) {
    list2env(sim_obj, envir = environment())

    sample_counts <- matrix(cell_counts, n_variants, n_pos)
    sim_obj$sample_counts <- sample_counts
    # draw true cell fluorescence values including synonymous mutations
    print("Drawing fluorescence distributions...")
    samples <- c()
    for (i in 1:n_variants) {
        for (j in 1:n_pos) {
            for (r in 1:n_replicates) {
                variant_samples <- cbind(sample_variant(sim_obj, i, j, r), i, j, r) # i, j, r is repeated for each sampeld fluorescence value
                samples <- rbind(samples, variant_samples)
            }
        }
    }
    sim_obj$samples <- samples

    # experimental process
    sim_obj <- sim_FACS(sim_obj)
    return(sim_obj)
}

