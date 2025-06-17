library(ggplot2)
library(dplyr)
source('../src/build_sim.R')
source('../src/fit_func.R')
source('../src/methods/method_utils.R')
source('../src/methods/diagnostics.R')
source('../src/plot_func.R')


args <- commandArgs(trailingOnly = TRUE)
modelfile <- args[1]
modelname <- tools::file_path_sans_ext(basename(modelfile))
draws <- readRDS(args[2])
dataname <- tools::file_path_sans_ext(basename(args[3]))
# data <- readRDS(args[3])
output_df_file <- args[3]
data <- readRDS(output_df_file)
plot_dir <- args[4]


K <- sum(startsWith(colnames(data), "c_")) # num bins


baseline_file <- paste0(tools::file_path_sans_ext(output_df_file), '_baseline.tsv')
param <- read_tsv(baseline_file)
q <- param[startsWith(param$param, "q"),2][[1]]
q <- q/sum(q)

sim <- posterior::as_draws_rvars(draws)

q_fit <- posterior::draws_of(sim$q)
saveRDS(q_fit, file=paste0(plot_dir, "/q_posterior.RData"))
if (!is.null(q_fit)) {
  q_means <- as.data.frame(colMeans(q_fit))
  colnames(q_means) <- "q_means"
  write_tsv(q_means, file = paste0(plot_dir, "/q_means.tsv"))
}

# ======= posterior histograms =========
print("Making plots...")

single_summaries <- list()
# sigma posterior
if (!is.null(posterior::draws_of(sim$sigma0))) {
  p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$sigma0)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of sigma0 for", dataname), x="sigma0", y="Freq")
  ggsave(paste0(plot_dir, "/sigma0_posterior.png"), p)
  single_summaries <- append(single_summaries, c(sigma0_mean=mean(posterior::draws_of(sim$sigma0)), sigma0_sd=sd(posterior::draws_of(sim$sigma0))))
}

if (!is.null(posterior::draws_of(sim$sigma))) {
  if (is.na(ncol(posterior::draws_of(sim$sigma)))) {
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$sigma)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of sigma for", dataname), x="sigma", y="Freq")
    ggsave(paste0(plot_dir, "/sigma_posterior.png"), p)
  } else {
    plot_df <- as.data.frame(posterior::draws_of(sim$sigma)) %>% pivot_longer(everything())
    plot_df$name <- factor(plot_df$name, levels=paste0("V", 1:ncol(posterior::draws_of(sim$sigma))))
    p <- ggplot(plot_df, aes(x=name, y=value)) + geom_violin() + stat_summary(fun=mean, geom="point", color="red", size=0.2) + theme_bw() +
        labs(title=paste("sigma posteriors for", dataname), x="Group ID", y="sigma")
    ggsave(paste0(plot_dir, "/sigma_posterior.png"), p)
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$sigma)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of sigma for", dataname), x="sigma", y="Freq")
    ggsave(paste0(plot_dir, "/all_sigma_posterior.png"), p)
    rm(plot_df)
  }
  rm(p)
}

if (!is.null(posterior::draws_of(sim$sigma_syn))) {
  if (is.na(ncol(posterior::draws_of(sim$sigma_syn)))) {
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$sigma_syn)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of sigma_syn for", dataname), x="sigma_syn", y="Freq")
    ggsave(paste0(plot_dir, "/sigma_syn_posterior.png"), p)
  } else {
    plot_df <- as.data.frame(posterior::draws_of(sim$sigma_syn)) %>% pivot_longer(everything())
    plot_df$name <- factor(plot_df$name, levels=paste0("V", 1:ncol(posterior::draws_of(sim$sigma_syn))))
    p <- ggplot(plot_df, aes(x=name, y=value)) + geom_violin() + stat_summary(fun=mean, geom="point", color="red", size=0.2) + theme_bw() +
        labs(title=paste("sigma_syn posteriors for", dataname), x="Group ID", y="sigma_syn")
    ggsave(paste0(plot_dir, "/sigma_syn_posterior.png"), p)
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$sigma_syn)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of sigma_syn for", dataname), x="sigma_syn", y="Freq")
    ggsave(paste0(plot_dir, "/all_sigma_syn_posterior.png"), p)
    rm(plot_df)
  }
  rm(p)
}

if (!is.null(posterior::draws_of(sim$sigma_pos))) {
  if (is.na(ncol(posterior::draws_of(sim$sigma_pos)))) {
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$sigma_pos)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of sigma_pos for", dataname), x="sigma_pos", y="Freq")
    ggsave(paste0(plot_dir, "/sigma_pos_posterior.png"), p)
  } else {
    plot_df <- as.data.frame(posterior::draws_of(sim$sigma_pos)) %>% pivot_longer(everything())
    plot_df$name <- factor(plot_df$name, levels=paste0("V", 1:ncol(posterior::draws_of(sim$sigma_pos))))
    p <- ggplot(plot_df, aes(x=name, y=value)) + geom_violin() + stat_summary(fun=mean, geom="point", color="red", size=0.2) + theme_bw() +
        labs(title=paste("sigma_pos posteriors for", dataname), x="Group ID", y="sigma_pos")
    ggsave(paste0(plot_dir, "/sigma_pos_posterior.png"), p)
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$sigma_pos)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of sigma_pos for", dataname), x="sigma_pos", y="Freq")
    ggsave(paste0(plot_dir, "/all_sigma_pos_posterior.png"), p)
    rm(plot_df)
  }
  rm(p)
}

# if (!is.null(posterior::draws_of(sim$tau))) { 
#     plot_df <- as.data.frame(posterior::draws_of(sim$tau)) %>% pivot_longer(everything())
#     plot_df$name <- factor(plot_df$name, levels=paste0("V", 1:ncol(posterior::draws_of(sim$phi))))
#     p <- ggplot(plot_df, aes(x=name, y=value)) + geom_violin() + stat_summary(fun=mean, geom="point", color="red", size=0.2) + theme_bw() +
#         labs(title=paste("tau posteriors for", dataname), x="Syn ID", y="tau")
#     ggsave(paste0(plot_dir, "/tau_posterior.png"), p)
#     rm(plot_df)
#     rm(p)
# }
if (!is.null(posterior::draws_of(sim$nu))) {
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$nu)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of nu for", dataname), x="sigma", y="Freq")
    ggsave(paste0(plot_dir, "/nu_posterior.png"), p)
    rm(p)
}
# phi posterior
if (grepl("nfunc", modelfile) && grepl("rep", modelfile)) {
  p <- ggplot(data, aes(x=n_counts, y=get(paste0(modelname, "_phi_mean")), color=rep)) + geom_line() + 
    geom_ribbon(aes(ymin=get(paste0(modelname, "_phi_mean")) - 2*get(paste0(modelname, "_phi_sd")), 
        ymax=get(paste0(modelname, "_phi_mean")) + 2*get(paste0(modelname, "_phi_sd")), fill=rep), alpha=0.1, linetype = 0) + 
    xlab("Observation counts") + ylab("Phi posterior mean") + ggtitle("Overdispersion vs observation counts") + theme_bw()
  ggsave(paste0(plot_dir, "/phi_posterior_lines.png"), p)
  rm(p)
}
if (!is.null(posterior::draws_of(sim$phi))) {
  if (is.na(ncol(posterior::draws_of(sim$phi)))) {
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$sigma)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("Posterior of phi for", dataname), x="phi", y="Freq")
  } else {
    if (grepl("nfunc", modelfile)) {
      p <- ggplot() + geom_histogram(aes(colMeans(posterior::draws_of(sim$phi))), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
          labs(title=paste("Posterior means of phi for", dataname), x="mean phi", y="freq")
    } else {
      plot_df <- as.data.frame(posterior::draws_of(sim$phi)) %>% pivot_longer(everything())
      plot_df$name <- factor(plot_df$name, levels=paste0("V", 1:ncol(posterior::draws_of(sim$phi))))
      p <- ggplot(plot_df, aes(x=name, y=value)) + geom_violin() + stat_summary(fun=mean, geom="point", color="red", size=0.2) + theme_bw() +
          labs(title=paste("phi posteriors for", dataname), x="Group ID", y="phi")
      rm(plot_df)
    }
  }
  ggsave(paste0(plot_dir, "/phi_posterior.png"), p)
  # random sample sigma-phi pairwise correlation histogram
  # indices <- sample.int(length(posterior::draws_of(sim$sigma), max(c(5000, length(posterior::draws_of(sim$sigma))))
  print("plotting cor mat")
  indices <- sample.int(4000, 100)

  if (!is.null(posterior::draws_of(sim$sigma))) {
    if (is.na(ncol(posterior::draws_of(sim$sigma))) && is.na(ncol(posterior::draws_of(sim$phi)))) {
        cor_mat <- cor(posterior::draws_of(sim$sigma)[indices], posterior::draws_of(sim$phi)[indices])
    } else if (is.na(ncol(posterior::draws_of(sim$sigma))) && !is.na(ncol(posterior::draws_of(sim$phi)))) {
        phi_indices <- sample.int(ncol(posterior::draws_of(sim$phi)), min(c(100, ncol(posterior::draws_of(sim$phi)))))
        print("phi indices")
        cor_mat <- cor(posterior::draws_of(sim$sigma)[indices], posterior::draws_of(sim$phi)[indices,phi_indices])
    } else if (!is.na(ncol(posterior::draws_of(sim$sigma))) && is.na(ncol(posterior::draws_of(sim$phi)))) {
        sigma_indices <- sample.int(ncol(posterior::draws_of(sim$sigma)), min(c(20, ncol(posterior::draws_of(sim$sigma)))))
        print("sigma indices")
        cor_mat <- cor(posterior::draws_of(sim$sigma)[indices,sigma_indices], posterior::draws_of(sim$phi)[indices])
    } else {
        print("phi+sigma indices")
        phi_indices <- sample.int(ncol(posterior::draws_of(sim$phi)), min(c(100, ncol(posterior::draws_of(sim$phi)))))
        sigma_indices <- sample.int(ncol(posterior::draws_of(sim$sigma)), min(c(20, ncol(posterior::draws_of(sim$sigma)))))
        cor_mat <- cor(posterior::draws_of(sim$sigma)[indices,sigma_indices], posterior::draws_of(sim$phi)[indices,phi_indices])
    }
    print(dim(cor_mat))
    cor_vals <- cor_mat[upper.tri(cor_mat)]
    p <- ggplot() + geom_histogram(aes(cor_vals), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("sigma-phi pairwise correlation values for", dataname), x="r^2", y="Freq")
    ggsave(paste0(plot_dir, "/sigma-phi_corr.png"), p)
    rm(cor_mat)
    rm(cor_vals)
    rm(p)
  }
}

if (!is.null(posterior::draws_of(sim$phi_groups))) {
  saveRDS(posterior::draws_of(sim$phi_groups), file=paste0(plot_dir, "/phi_groups_posterior.RData"))
}

if (!is.null(posterior::draws_of(sim$a))) {
  for (i in 1:ncol(posterior::draws_of(sim$a))) {
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$a)[,i]), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("phi coefficient", dataname), x="a", y="Freq")
    ggsave(paste0(plot_dir, "/a_hist", i, ".png"), p)
    rm(p)
    single_summaries[[paste0("a_mean", i)]] <- mean(posterior::draws_of(sim$a)[,i])
    single_summaries[[paste0("a_sd", i)]] <- sd(posterior::draws_of(sim$a)[,i])
    # single_summaries <- append(single_summaries, c(a_mean=, a_sd=sd(posterior::draws_of(sim$a)[,i])))
  }
  saveRDS(posterior::draws_of(sim$a), file=paste0(plot_dir, "/a_posterior.RData"))
}
if (!is.null(posterior::draws_of(sim$b))) {
  for (i in 1:ncol(posterior::draws_of(sim$b))) {
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$b)[,i]), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("phi intercept", dataname), x="b", y="Freq")
    ggsave(paste0(plot_dir, "/b_hist", i, ".png"), p)
    rm(p)
    single_summaries[[paste0("b_mean", i)]] <- mean(posterior::draws_of(sim$b)[,i])
    single_summaries[[paste0("b_sd", i)]] <- sd(posterior::draws_of(sim$b)[,i])
    # single_summaries <- append(single_summaries, c(b_mean=mean(posterior::draws_of(sim$b)[,i]), b_sd=sd(posterior::draws_of(sim$b)[,i])))
  }
  saveRDS(posterior::draws_of(sim$b), file=paste0(plot_dir, "/b_posterior.RData"))
}
print("c sim")
if (!is.null(posterior::draws_of(sim$c))) {
  for (i in 1:ncol(posterior::draws_of(sim$c))) {
    p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$c)[,i]), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
    labs(title=paste("quadratic phi coefficient", dataname), x="c", y="Freq")
    ggsave(paste0(plot_dir, "/c_hist", i, ".png"), p)
    rm(p)
    single_summaries[[paste0("c_mean", i)]] <- mean(posterior::draws_of(sim$c)[,i])
    single_summaries[[paste0("c_sd", i)]] <- sd(posterior::draws_of(sim$c)[,i])
    # single_summaries <- append(single_summaries, c(c_mean=mean(posterior::draws_of(sim$c)[,i]), c_sd=sd(posterior::draws_of(sim$c)[,i])))
  }
  saveRDS(posterior::draws_of(sim$c), file=paste0(plot_dir, "/c_posterior.RData"))
}


# means of mu posteriors
p <- ggplot() + geom_histogram(aes(posterior::draws_of(sim$mu)), bins=50, color="black", fill="#76B7B2") + theme_bw() + 
  labs(title=paste("Posterior means of mu for", dataname), x="mu", y="Freq")
ggsave(paste0(plot_dir, "/mu_means.png"), p)
rm(p)
saveRDS(posterior::draws_of(sim$mu), file=paste0(plot_dir, "/mu_posterior.RData"))
saveRDS(posterior::draws_of(sim$theta), file=paste0(plot_dir, "/theta_posterior.RData"))
saveRDS(posterior::draws_of(sim$theta_syn), file=paste0(plot_dir, "/theta_syn_posterior.Rdata"))


# total counts by position plot
pos_coverage_df <- data %>% group_by(position) %>% summarize(mean_count=mean(n_counts), sd_count=sd(n_counts))

p <- ggplot(pos_coverage_df, aes(x=position, y=mean_count)) + 
  geom_point(color="red", size=0.5) +
  geom_linerange(aes(ymin=mean_count-sd_count, ymax=mean_count+sd_count), size=0.1, color="darkgreen") +
  labs(title="n mean counts by position", x="pos", y="mutant mean counts") + theme_bw()
ggsave(paste0(plot_dir, "/counts_by_pos.png"), p)
rm(p)
rm(pos_coverage_df)

yaml::write_yaml(single_summaries, paste0(plot_dir, "/other_param_summaries.yaml"))
# =========== posterior predictive ==========
print("Computing posterior predictive CI scores...")
# posterior predictive by variant by bin credible interval based on model (TODO: modularize and make based on which model)
# for each variant
  # for each replicate
    # for each bin
      # check if true bin count in credible interval from posterior predictive (+1 if not)
  # score variant

n_sims <- 1000
scores <- posterior_samples(n_sims, data, sim, modelname, K=K) # score each variant obs
if (length(scores) > 0) {
  data$scores <- scores
  data <- data %>% group_by(hgvs_exp) %>% mutate(variant_score=sum(scores)) # add scores across replicates
  p <- ggplot(data, aes(x=variant_score)) + geom_bar(fill="#76B7B2") + theme_bw() +
    labs(title="variant bin count outside 95 pos. pre. CI across reps", x="Num bins CI", y="Num variants")
  ggsave(paste0(plot_dir, "/postpred_CI.png"), p)
  rm(p)
}
saveRDS(scores, paste0(tools::file_path_sans_ext(output_df_file), '_posterior_pred_scores.RData'))
rm(scores)

# effect matrix plots
df <- readRDS(output_df_file)
modelname <- tools::file_path_sans_ext(basename(modelfile))

if (!grepl("nopos", modelname)) {
  df2tsv <- df[c("hgvs", "exp", "wildtype", "position", "mutation", "type", 
                  paste0(modelname, c("_syn_recalibrated_mu_mean", "_syn_recalibrated_mu_sd", "_syn_recalibrated_mu_lfsr",
                  "_pos_mean", "_pos_sd", "_sigma_mean", "_sigma_sd")))]
  colnames(df2tsv) <- c("hgvs", "exp", "wildtype", "position", "mutation", "type", 
                  "effect", "effect_se", "lfsr", 
                  "pos_effect", "pos_effect_se", "pos_sd", "pos_sd_se")
} else {
  df2tsv <- df[c("hgvs", "exp", "wildtype", "position", "mutation", "type", 
                  paste0(modelname, c("_syn_recalibrated_mu_mean", "_syn_recalibrated_mu_sd", "_syn_recalibrated_mu_lfsr")))]
  colnames(df2tsv) <- c("hgvs", "exp", "wildtype", "position", "mutation", "type", 
                  "effect", "effect_se", "lfsr")
}

df2tsv <- distinct(df2tsv)
df2tsv$position <- as.numeric(df2tsv$position)
df2tsv <- df2tsv[order(df2tsv$position, df2tsv$type),]
write_tsv(df2tsv, paste0(plot_dir, "/variant_scores.tsv"))

df2tsv <- df2tsv %>% tidyr::complete(position = min(df2tsv$position):max(df2tsv$position))
df2tsv[is.na(df2tsv$hgvs),]$mutation <- "A" # TODO: fix this to be completely blank
df2tsv[is.na(df2tsv$hgvs),]$type <- "-"
df2tsv <- df2tsv[order(df2tsv$position),]

df2tsv <- df2tsv %>% arrange(str_length(mutation) > 1, mutation)

scoreHeatmap(df2tsv, plot_dir, score.col="effect", name="score_heatmap", x.text=2.5, seq.text=0.75, y.text=2)
scoreDensity(df2tsv, plot_dir, score.col="effect", name="score_histogram", hist=T, scale.free=T)

# plot 3D data colored by score (if K=4)
if (K == 4 && any(df[[paste0(modelname, "_mu_mean")]] < 0) && any(df[[paste0(modelname, "_mu_mean")]] > 0)) {
  print("plotting 3D")
  simplex <- function(n) {
    qr.Q(qr(matrix(1, nrow=n)) ,complete = TRUE)[,-1]
  }
  tetra <- simplex(4)
  counts <- df  %>% ungroup() %>% select(c("c_0", "c_1", "c_2", "c_3"))
  type <- df$type
  hgvs <- df$hgvs
  props <- cbind(counts / df$n_counts, type, hgvs)
  print("converting coordinates")
  props_3D <- geometry::bary2cart(tetra, as.matrix(props[,1:4]))
  props_3D_df <- as.data.frame(props_3D)
  props_3D_df$type <- df$type
  props_3D_df$score <- df[[paste0(modelname, "_mu_mean")]]
  props_3D_df$mean_counts <- df$mean_counts
  props_3D_df$position <- df$position

  q_coord <- geometry::bary2cart(tetra, q)
  # create simplex 3d plot
  library(plotly)
  print("Making colorscale")
  colorlength <- 100
  zeroval <- 0 - min(props_3D_df$score) / (max(props_3D_df$score) - min(props_3D_df$score)) # Note: breaks if no values below zero
  border <- as.integer(zeroval * colorlength)
  colorscale <- as.list(1:colorlength)
  #colorscale below zero
  s <- scales::seq_gradient_pal("#2166ac", "#FFFFFF", "Lab")(seq(0,1,length.out=border))
  for (i in 1:border) {
      colorscale[[i]] <- c((i - 1) / colorlength, s[i])
  }
  #colorscale above zero
  s <- scales::seq_gradient_pal("#FFFFFF", "#b2182b", "Lab")(seq(0,1,length.out=colorlength - border))
  for (i in 1:(colorlength - border)) {
    colorscale[[i + border]] <- c((i + border) / colorlength, s[i])
  }

  axx <- list(
    showticklabels=F,
    title = " "
  )
  axy <- list(
    showticklabels=F,
    title = " "
  )
  axz <- list(
    showticklabels=F,
    title = " "
  )
  print("Making 3D plot")
  print(nrow(props_3D_df))
  fig <- plot_ly() %>% add_trace(data=props_3D_df, x=~V1, y=~V2, z=~V3, color=~score,
                               type = "scatter3d", mode = "markers", marker=list(size=1.5), showlegend=F,
                               hoverinfo="text",
                               hovertext = paste(paste0('<b>', df$hgvs, '</b>'),
                                                 paste0('<br> score: ', props_3D_df$score),
                                                 paste0('<br><i> c_0 </i>: ', df$c_0),
                                                 paste0('<br><i> c_1 </i>: ', df$c_1),
                                                 paste0('<br><i> c_2 </i>: ', df$c_2),
                                                 paste0('<br><i> c_3 </i>: ', df$c_3)
                                                     )) %>%
                      add_trace(data=as.data.frame(tetra[c(1,2,3,4,1,3,1,2,4),]), 
                                x=~V1, y=~V2, z=~V3, type="scatter3d", mode="lines", line=list(color="darkgray", width=3), showlegend=F, hoverinfo='none') %>% 
                      add_trace(x=tetra[,1], y=tetra[,2], z=tetra[,3], type="scatter3d", mode="text", 
                                text=c("c_0 <br> (1,0,0,0)", "c_1 <br> (0,1,0,0)", "c_2 <br> (0,0,1,0)", "c_3 <br> (0,0,0,1)"), 
                                inherit=F, showlegend=F, hoverinfo='none') %>%
                      add_trace(x=q_coord[1], y=q_coord[2], z=q_coord[3], type="scatter3d", mode="markers", name = 'Synonymous baseline', hoverinfo="text",
                                hovertext = paste("Synonymous baseline", "<br> <i>", paste(q, collapse=" "), "</i>"), marker=list(size=10, color="pink", opacity=0.8)) %>%
                      layout(title="Variants in proportion space", showlegend=T, scene = list(xaxis=axx,yaxis=axy,zaxis=axz)) %>%
                      add_annotations(
                        text = "Axes on triangle represent proportion for that bin",
                        x = 1,
                        y = 1,
                        yref = "paper",
                        xref = "paper",
                        xanchor = "middle",
                        yanchor = "top",
                        showarrow = FALSE,
                        font = list(size = 15)
                      )
  print("Saving plot...")
  htmlwidgets::saveWidget(as_widget(fig), paste0(plot_dir, "/variant_props_by_score.html"))
} else {
  file.create(paste0(plot_dir, "/variant_props_by_score.html")) # create empty file for the sake of snakemake
}