# Functions to fit methods on simulation

library(tidyverse)
library(rstan)
library(cmdstanr)
library(dplyr)
library(ggplot2)


# fit stan model baseline
fit_baseline <- function(obs_df, model_file, baseline_file, G_size=100, breaks=NULL) {
    # fit stan model baseline
    data_control <- obs_df %>% filter(type == "synonymous")
    data_control  <- data_control %>% mutate(across(starts_with("c_"), ~ .x + 1))
    N <- nrow(data_control)
    sMAPr <- as.numeric(factor(data_control$rep))
    R <- length(unique(data_control$rep))
    y <- as.matrix(data_control %>% dplyr::select(starts_with("c_")))
    K <- ncol(y)
    if (basename(model_file) == "syn_groups_phi.stan" ||
        basename(model_file) == "syn_groups_sigma.stan" ||
        basename(model_file) == "baseline_grouped.stan") {
      data_control <- data_control %>% group_by(hgvs_exp) %>% mutate(total_counts=sum(n_counts))
      data_control <- data_control %>% group_by(hgvs_exp) %>% mutate(mean_counts=mean(n_counts))
      # nMAPg <- as.numeric(cut_number(data_control$total_counts, as.integer(length(data_control$hgvs_exp) / G_size)))
      nMAPg <- as.numeric(cut_number(data_control$mean_counts, as.integer(length(data_control$hgvs_exp) / G_size)))
      G <- max(nMAPg)
      print(G)
      input <- list(y = y, N = N, G=G, nMAPg=nMAPg, K = K)
    } else if (basename(model_file) == "baseline_nfunc.stan") {
      print("using baseline nfunc")
      n_counts <- data_control$n_counts
      input <- list(y = y, N = N, K = K, R=R, n_counts=n_counts, sMAPr=sMAPr)
    } else if (is.null(breaks)) {
        input <- list(y = y, N = N, K = K)
    } else {
        input <- list(y = y, N = N, K = K, breaks=breaks)
    }
    
    rm(N, K, y)

    mod <- cmdstan_model(model_file)
    fit <- mod$sample(
    data = input,
    chains = 4,
    parallel_chains = 4,
    refresh = 500
    )
    stanfit <- rstan::read_stan_csv(fit$output_files())

    param <- summary(stanfit, pars = c("q"), probs = c(0.025, 0.25, 0.5, 0.75, 0.975))$summary
    param <- data.frame(param)
    param$param <- rownames(param)
    param <- param %>% dplyr::select(11, all_of(1:10))
    write_tsv(param, file = baseline_file)
    return(stanfit)
}

fit_overall <- function(obs_df, model_file, 
                        stanfit_file, res_file, input_file, G_size=100, selfphi=NULL,
                        seed=NULL, noinit=F) {
    pseudocount <- T # use pseudocount
    data_exp <- obs_df %>% unite("hgvs_exp", hgvs, exp, remove = FALSE)
    if (pseudocount) {
      data_exp  <- data_exp %>% mutate(across(starts_with("c_"), ~ .x + 1))
    }
    
    V <- length(unique(data_exp$hgvs_exp))
    N <- nrow(data_exp)
    nMAPv <- as.numeric(factor(data_exp$hgvs_exp))
    n_counts <- as.numeric(data_exp$n_counts)
    n_syn_counts <- as.numeric(data_exp[data_exp$type=="synonymous",]$n_counts)
    distinct_df <- (distinct(data.frame(hgvs_exp = data_exp$hgvs_exp, nMAPv = nMAPv, 
                                total_counts = data_exp$total_counts, type=data_exp$type)) %>%
              arrange(nMAPv))
    n_variant_counts <- distinct_df[["total_counts"]]
    n_syn_variant_counts <- distinct_df[distinct_df$type=="synonymous",][["total_counts"]]

    C <- length(unique(data_exp$exp))
    vMAPc <- distinct(data.frame(hgvs_exp = data_exp$hgvs_exp, nMAPv = nMAPv, 
                                cond = data_exp$exp)) %>%
              arrange(nMAPv)
    vMAPc <- as.numeric(factor(vMAPc$cond))

    vMAPp <- distinct(data.frame(hgvs_exp = data_exp$hgvs_exp, nMAPv = nMAPv, type = data_exp$type, position = data_exp$position)) %>%
            arrange(nMAPv)
    vMAPp$position <- as.numeric(vMAPp$position)
    # group position w/ < 10 variants with next position
    positions <- sort(unique(vMAPp$position))
    print("grouping positions")
    print(head(positions))
    for (i in 1:length(positions)) {
      if (i != length(positions)) {
        # print(i)
        # print(sum(vMAPp$position == positions[i]))
        if (sum(vMAPp$position == positions[i]) < 10) {
          vMAPp[vMAPp$position == positions[i],]$position <- positions[i+1]
        }
      }
    }
    # ensure end also has 10 variants
    while (sum(vMAPp$position == positions[i]) < 10) {
      vMAPp[vMAPp$position == positions[i],]$position <- positions[i-1]
      i <- i - 1
    }
    print("setting syn position")
    print(head(vMAPp))
    print(vMAPp[vMAPp$type=="synonymous",]$position)
    vMAPp[vMAPp$type=="synonymous",]$position <- 0 # encode synonymous as diff position
    S <- sum(vMAPp$type=="synonymous")
    vMAPs <- rep(-1, nrow(vMAPp)) # -1 for non syn, o.w. each syn gets own index
    vMAPs[vMAPp$type=="synonymous"] <- 1:S


    vMAPp <- as.numeric(as.factor(vMAPp$position))
    P <- length(unique(vMAPp))
    
    R <- length(unique(data_exp$rep))
    nMAPr <- as.numeric(factor(data_exp$rep))
    sMAPr <- nMAPr[data_exp$type=="synonymous"]
    N_syn <- length(sMAPr)

    y <- as.matrix(data_exp %>% dplyr::select(starts_with("c_")))
    y_syn <- y[data_exp$type=="synonymous",]
    K <- ncol(y)

    # compute count groups
    data_exp <- data_exp %>% group_by(hgvs_exp) %>% mutate(total_counts=sum(n_counts))
    data_exp <- data_exp %>% group_by(hgvs_exp) %>% mutate(mean_counts=mean(n_counts))

    # nMAPg -- make count groups within each replicate
    nMAPg <- rep(NA, nrow(data_exp))
    group_index_offset <- 0
    for (r in 1:R) {
      ngroups_in_rep <- as.integer(sum(nMAPr == r) / G_size)
      within_replicate_group_index <- as.numeric(santoku::chop_equally(data_exp[nMAPr == r,]$n_counts, ngroups_in_rep))
      nMAPg[nMAPr == r] <- within_replicate_group_index + group_index_offset
      group_index_offset <- max(nMAPg[!is.na(nMAPg)])
    }
    # sMAPg -- map synonymous to count groups 
    sMAPg <- nMAPg[data_exp$type=="synonymous"]
    print(table(nMAPg))
    G <- max(nMAPg)

    init <- 1
    if (basename(model_file) == "FACS_double_sample_repq.stan") {
      input <- list(V = V, S = S, N = N, N_syn=N_syn, P=P, R=R, nMAPv = nMAPv, vMAPs = vMAPs, nMAPr = nMAPr, vMAPp = vMAPp, sMAPr=sMAPr, 
                n_counts=n_counts, n_syn_counts=n_syn_counts, K = K, y = y, y_syn = y_syn)
      init <- list(list(sigma=rep(1, P-1)), list(sigma=rep(1, P-1)), list(sigma=rep(1, P-1)), list(sigma=rep(1, P-1)))
    } else if (basename(model_file) == "FACS_double_sample_repq_nopos.stan") {
      input <- list(V = V, S = S, N = N, N_syn=N_syn, P=P, R=R, nMAPv = nMAPv, vMAPs = vMAPs, nMAPr = nMAPr, vMAPp = vMAPp, sMAPr=sMAPr, 
                n_counts=n_counts, n_syn_counts=n_syn_counts, K = K, y = y, y_syn = y_syn)
      init <- list(list(sigma=rep(1, V)), list(sigma=rep(1, V)), list(sigma=rep(1, V)), list(sigma=rep(1, V)))
    } else if (basename(model_file) == "FACS_double_sample_count_grouped.stan") {
      input <- list(V = V, S = S, N = N, N_syn=N_syn, P=P, R=R, G=G, nMAPv = nMAPv, vMAPs = vMAPs, nMAPg = nMAPg, nMAPr = nMAPr, vMAPp = vMAPp, sMAPr=sMAPr, sMAPg=sMAPg,
                n_counts=n_counts, n_syn_counts=n_syn_counts, K = K, y = y, y_syn = y_syn)
      init <- list(list(sigma=rep(1, P-1)), list(sigma=rep(1, P-1)), list(sigma=rep(1, P-1)), list(sigma=rep(1, P-1)))
    } else {
      input <- list(V = V, N = N, C = C, nMAPv = nMAPv, vMAPc = vMAPc,
                K = K, y = y)
    }
    if (noinit) {
      init <- NULL
    }
    mod <- cmdstan_model(model_file)
    
    fit <- mod$sample(
    data = input,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    seed = seed,
    init=init
    )
    saveRDS(input, file = input_file)
    print("Summarizing posteriors")
    # summarize
    mu <- fit$draws(variables = "mu", format = "matrix")
    mu_quantile <- apply(mu, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    mu_mean <- apply(mu, 2, mean)
    mu_sd <- apply(mu, 2, sd)
    mu_lfsr <- apply(mu, 2, function(samples) {
      p_less <- sum(samples < 0) / length(samples)
      p_great <- sum(samples > 0) / length(samples)
      min(c(p_less, p_great))
    })
    
    # compute syn recalibrated scores by sampling from random synonymous
    syn_mu <- as.matrix(mu[,unique(input$nMAPv[data_exp$type=="synonymous"])])
    syn_indices <- t(replicate(nrow(mu), sample.int(ncol(syn_mu), ncol(mu), replace=T)))
    mu2 <- matrix(NA, nrow=nrow(mu), ncol=ncol(mu))
    for (i in 1:nrow(mu)) {
      mu2[i,] <- as.numeric(mu[i,]) - as.numeric(syn_mu[i, syn_indices[i,]])
    }
    syn_recalibrated_mu_mean <- apply(mu2, 2, mean)
    syn_recalibrated_mu_sd <- apply(mu2, 2, sd)
    syn_recalibrated_mu_lfsr <- apply(mu2, 2, function(samples) {
      p_less <- sum(samples <= 0) / length(samples)
      p_great <- sum(samples >= 0) / length(samples)
      min(c(p_less, p_great))
    })

    mu_stats <- data.frame(t(rbind(mu_mean, mu_sd, mu_lfsr, 
                syn_recalibrated_mu_mean, syn_recalibrated_mu_sd, syn_recalibrated_mu_lfsr)))
    mu_stats$param <- rownames(mu_stats)
    mu_stats <- mu_stats %>% filter(substr(param, 1, 2) == "mu")
    
    modelname <- tools::file_path_sans_ext(basename(model_file))
    colnames(mu_stats) <- paste(modelname, colnames(mu_stats), sep="_")
    mu_stats$map <- 1:input$V

    # add position mean + sigma to df
    if (grepl("syngrouped", model_file)) {
      theta <- fit$draws(variables = "theta_pos", format = "matrix")
      sigma <- fit$draws(variables = "sigma_pos", format = "matrix")
    } else {
      theta <- fit$draws(variables = "theta", format = "df")
      sigma <- fit$draws(variables = "sigma", format = "df")
    }
    pos_mean <- apply(theta, 2, mean)
    pos_sd <- apply(theta, 2, sd)
    sigma_mean <- apply(sigma, 2, mean)
    sigma_sd <- apply(sigma, 2, sd)
    pos_lfsr <- apply(theta, 2, function(samples) {
      p_less <- sum(samples <= 0) / length(samples)
      p_great <- sum(samples >= 0) / length(samples)
      min(c(p_less, p_great))
    })
    if (!grepl("nopos", model_file)) {
      mu_stats$p_map <- input$vMAPp
      pos_stats <- data.frame(t(rbind(pos_mean, pos_sd, sigma_mean, sigma_sd, pos_lfsr)))
      pos_stats$pos_param <- rownames(pos_stats)
      pos_stats <- pos_stats %>% filter(substr(pos_param, 1, 5) == "theta")
      colnames(pos_stats) <- paste(modelname, colnames(pos_stats), sep="_")
      # if (grepl("syngrouped", model_file) || grepl("syndist", model_file)) {
      pos_stats$p_map <- 2:input$P
      mu_stats <- distinct(mu_stats %>% left_join(pos_stats))
    }

    if (!is.null(input$nMAPg)) {
      n_pos <- data.frame(hgvs_exp = data_exp$hgvs_exp,
                        rep = data_exp$rep,
                        map = input$nMAPv,
                        gmap = input$nMAPg)
    } else {
      n_pos <- data.frame(hgvs_exp = data_exp$hgvs_exp,
                        rep = data_exp$rep,
                        map = input$nMAPv)
    }

    n_pos <- distinct(n_pos %>% left_join(mu_stats))
    # remove pseudocount
    if (pseudocount) {
      data_exp  <- data_exp %>% mutate(across(starts_with("c_"), ~ .x - 1))
    }
    data_exp <- data_exp %>% left_join(n_pos)

    if (is.null(input$nMAPg)) {
      a <- fit$draws(variables = "a", format = "matrix")
      b <- fit$draws(variables = "b", format = "matrix")
      a_mean <- apply(a, 2, mean)
      a_sd <- apply(a, 2, sd)
      b_mean <- apply(b, 2, mean)
      b_sd <- apply(b, 2, sd)
      # phi <- fit$draws(variables = "phi", format = "matrix")
      # phi_mean <- apply(phi, 2, mean)
      # phi_sd <- apply(phi, 2, sd)
      phi_stats <- data.frame(t(rbind(a_mean, a_sd, b_mean, b_sd)))
      phi_stats$param <- rownames(phi_stats)
      print(phi_stats)
      phi_stats <- phi_stats %>% filter(substr(param, 1, 1) == "a" | substr(param, 1, 1) == "b")
      colnames(phi_stats) <- paste(modelname, colnames(phi_stats), sep="_")
      readr::write_tsv(phi_stats, file = paste0(dirname(res_file), "/", modelname, "_phi_stats.tsv"))

      reps <- unique(data_exp$rep)
      data_exp[[paste0(modelname, "_phi_mean")]] <- NA
      for (i in 1:length(reps)) {
        r <- reps[i]
        data_exp[data_exp$rep==r,paste0(modelname, "_phi_mean")] <- a_mean[i] * data_exp[data_exp$rep==r,]$n_counts + b_mean[i]
        data_exp[data_exp$rep==r,paste0(modelname, "_phi_sd")] <- sqrt(a_sd[i]^2 * data_exp[data_exp$rep==r,]$n_counts^2)
      }
    }
    saveRDS(data_exp, file = res_file)
    print("Reading/saving stanfit")
    samples <- fit$draws(format = "df")
    sampler_diagnostics <- fit$diagnostic_summary()
    rm(fit)
    saveRDS(samples, file = stanfit_file)
    return(sampler_diagnostics)
}


compare_syn_posterior <- function(data_control, G_vec, model_file, param, path, dataname) {
  range_p_list <- vector('list', length(G_vec))
  param_sample_list <- vector('list', length(G_vec))
  for (i in 1:length(G_vec)) {
    G_size <- G_vec[i]
    if (is.na(G_size)) {
      G_size <- length(data_control$hgvs_exp)
    }
    nMAPg <- as.numeric(cut_number(data_control$total_counts, as.integer(length(data_control$hgvs_exp) / G_size)))
    G <- max(nMAPg)
    # plot group count ranges
    G_ranges <- matrix(NA, G, 3)
    for(g in 1:G) {
      obs_counts <- data_control[nMAPg==g,]$total_counts
      G_ranges[g,] <- c(g, range(obs_counts))
    }
    p <- ggplot(as.data.frame(G_ranges), aes(x=V1)) + geom_linerange(aes(ymin=V2,ymax=V3),linetype=2,color="blue")+
      geom_point(aes(y=V2),size=0.1,color="red") +
      geom_point(aes(y=V3),size=0.1,color="red") +
      theme_bw() + labs(title=paste("Groups of", G_size), x="Group ID", y="Count Range")
    range_p_list[[i]] <- p
    
    G_param_stanfit <- fit_baseline(data_control %>% ungroup(), model_file, paste0(path, "/", "baseline_summary.tsv"), G_size=G_size)
    sim_G <- rstan::extract(G_param_stanfit)
    param_sample_list[[i]] <- sim_G[[param]]
  }
  p_ranges <- gridExtra::grid.arrange(grobs=range_p_list,nrow=2, top = "Variant total count ranges")
  ggsave(paste0(path, "/", dataname, "_syn_group_ranges", ".png"), p_ranges)
  
  phi_plot_list <- lapply(1:length(param_sample_list), function(i) {
    phi_G_samples <- param_sample_list[[i]]
    plot_df <- as.data.frame(phi_G_samples) %>% pivot_longer(everything())
    plot_df$name <- factor(plot_df$name, levels=paste0("V", 1:ncol(phi_G_samples)))
    ggplot(plot_df, aes(x=name, y=value)) + geom_violin() + stat_summary(fun.y=mean, geom="point", color="red", size=0.2) + theme_bw() +
      labs(title=paste("Groups of", G_vec[i]), x="Group ID", y=param)
  })
  p_exp <- gridExtra::grid.arrange(grobs=phi_plot_list, top = paste("Phi posterior for", dataname))
  ggsave(paste0(path, "/", dataname, "_syn_group_", param, ".png"), p_exp)
  return(param_sample_list)
}

# fit other methods
shep_moments <- function(counts, breaks) {
  breaks_lower <- breaks[1:(length(breaks)-1)]
  breaks_upper <- breaks[2:(length(breaks))]
  breaks_mid <- (breaks_upper + breaks_lower)/2
  
  n <- sum(counts)
  mu <- sum(breaks_mid * counts) / n
  sigma2 <- (sum(breaks_mid^2 * counts) - n * mu^2) / (n-1)
  return(c(mu, sigma2))
}

cdf_lik_log <- function(theta, counts, breaks, cdf=pnorm) {
  breaks_lower <- breaks[1:(length(breaks)-1)]
  breaks_upper <- breaks[2:(length(breaks))]
  # moments of distribution
  mu <- theta[1] 
  sigma <- exp(theta[2])
  # binned NLL
  return(-sum(sapply(1:length(counts), function(i) {
    counts[i] * log(cdf(breaks_upper[i], mu, sigma) - cdf(breaks_lower[i], mu, sigma))
  })))
}

cdf_lik_log_multi <- function(theta, counts_mat, breaks, cdf=pnorm) {
  breaks_lower <- breaks[1:(length(breaks)-1)]
  breaks_upper <- breaks[2:(length(breaks))]
  # moments of distribution
  mu <- theta[1] 
  sigma <- exp(theta[2])
  # binned NLL
  NLL <- 0
  if (is.null(nrow(counts_mat))) {
    counts_mat <- t(as.matrix(counts_mat))
  } else {
    counts_mat <- as.matrix(counts_mat)
  }
  for (i in 1:nrow(counts_mat)) {
    counts <- counts_mat[i,]
    NLL <- NLL - sum(sapply(1:length(counts), function(i) {
    counts[i] * log(cdf(breaks_upper[i], mu, sigma) - cdf(breaks_lower[i], mu, sigma))
    }))
  }
  return(NLL)
}

dirmult_NLL <- function(alpha, counts) {
  phi <- sum(alpha)
  n <- rowSums(counts)
  left <- sum(lgamma(phi) + lgamma(n+1) - lgamma(n+phi))
  right <- sum(sweep(counts, 2, alpha, function(c,al) lgamma(c + al) - lgamma(al) - lgamma(c+1)))
  return(-(left+right))
}

dirmult_NLL_fixed_alpha <- function(phi, q, counts) {
  alpha <- phi*q
  n <- rowSums(counts)
  left <- sum(lgamma(phi) + lgamma(n+1) - lgamma(n+phi))
  right <- sum(sweep(counts, 2, alpha, function(c,al) lgamma(c + al) - lgamma(al) - lgamma(c+1)))
  return(-(left+right))
}

dirmult_NLL_constrained <- function(alpha, counts) {
  alpha <- exp(alpha)
  phi <- sum(alpha)
  n <- rowSums(counts)
  left <- sum(lgamma(phi) + lgamma(n+1) - lgamma(n+phi))
  right <- sum(sweep(counts, 2, alpha, function(c,al) lgamma(c + al) - lgamma(al) - lgamma(c+1)))
  return(-(left+right))
}

dirmult_NLL_onlyphi <- function(phi, num_variants, n_counts) {
  alpha <- phi*rep(1/num_variants, num_variants)
  n <- sum(n_counts)
  left <- sum(lgamma(phi) + lgamma(n+1) - lgamma(n+phi))
  right <- sum(lgamma(n_counts + alpha) - lgamma(alpha) - lgamma(n_counts+1))
  return(-(left+right))
}