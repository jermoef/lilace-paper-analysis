# Functions for plotting
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggpubr)
library(utils)
library(readr)
library(stringr)
library(patchwork)
library(cowplot)
library(ggthemes)

build_plot_effects_df <- function(data_exp, loc_shift_mat, c) {
    dd_est <- data.frame()
    for (i in 1:nrow(loc_shift_mat)) {
        for (j in 1:ncol(loc_shift_mat)) {
            # mut == i+1 bc i=1 is synonymous
            dd_est <- rbind(dd_est, data_exp[data_exp$position==j & 
                                    data_exp$mutation==(i+1) & 
                                    data_exp$exp==c,] %>% 
                                    select('position', 'mutation', 'mu_mean') %>%
                                    distinct())
        }
    }
    names(dd_est) <- c('x', 'y', 'col')
    dd_est$x <- as.integer(dd_est$x)
    dd_est$y <- as.integer(dd_est$y)
    dd_est$col <- dd_est$col
    return(dd_est)
}


#' Makes grid plot of effect sizes by position and mutation
#'  first col is position, second is mutation, third is effect size
plot_effects <- function(df, color_breaks) {
    names(df) <- c('x', 'y', 'col')
    colors <- c("maroon", "white", "steelblue")
    p <- ggplot(df, aes(x = x, y = y)) + 
        geom_tile(aes(fill = col),colour = "white") +
        scale_fill_gradientn(
            limits  = range(df$col),
            colors = colors[c(1, seq_along(colors), length(colors))],
            values  = c(0, scales::rescale(color_breaks, from = range(df$col)), 1)) +
        theme_bw() + 
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_discrete(expand = c(0, 0))
    return(p)
}


# fxns lifted from rosace / modified slightly
scoreHeatmap <- function(data,
                         savedir,
                         ctrl.name = "synonymous",
                         pos.col = "position",
                         wt.col = "wildtype",
                         mut.col = "mutation",
                         type.col = "type",
                         score.col = "mean",
                         aa.order = NA,
                         npos = 100,
                         ncol = 1,
                         pos.step = 5,
                         x.text = 6,
                         y.text = 3,
                         seq.text = 1.1,
                         c.pallete = 'RdBu',
                         c.scale = list(),
                         ht = 11,
                         wd = 8.5,
                         name = "Heatmap",
                         savepdf = TRUE,
                         savesvg = TRUE,
                         show = FALSE,
                         factor_score =FALSE,
                         discovery_score=FALSE,
                         compare_score=FALSE,
                         category_score=FALSE,
                         cat_name="discovery",
                         flipped=FALSE
){
  
  # determine certain plot properties
  if (is.na(aa.order)) {
    aa.order = unique(na.omit(data[[mut.col]]))
  }
  pos.order <- levels(factor(data[[pos.col]]))
  # pos.order <- min(data[[pos.col]]):max(data[[pos.col]])
  npanel <- length(pos.order) %/% npos +1
  nrow <- npanel / ncol
  if (flipped) {
    c.default <- list(palette = c.pallete, mid = 0, rev=FALSE, na.value = '#E5E5E5')
  } else {
    c.default <- list(palette = c.pallete, mid = 0, rev=TRUE, na.value = '#E5E5E5')
  }
  c.args <- modifyList(c.default, c.scale)
  
  # parse the positions
  starts <- seq(1, length(pos.order), by = npos)
  if (length(starts) > 1) {
    ends <- c((starts - 1)[2:length(starts)], length(pos.order))
  } else {
    ends <- length(pos.order)
  }
  
  plot_list <- lapply(1:length(starts), function(i) {
    
    start <- starts[i]
    end <- ends[i]
    legend = ifelse(i == length(starts), 'right', 'none')
    
    # subset data for positions in range
    sub_data <- data[data[[pos.col]] %in% pos.order[start:end], ]
    sub_pos.order <- levels(factor(sub_data[[pos.col]]))
    
    # create subplot heatmaps
    if (factor_score) {
      colorscale <- c("1" = "#5778a4", "2" = "#e49444", "3" = "#6a9f58", "4" = "#967662", "5" = "#a87c9f", "6" = "#d62727")
       p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]), 
        y = factor(.data[[mut.col]], levels = aa.order), fill = factor(.data[[score.col]]))) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) + 
        scale_fill_manual(values=colorscale, drop=FALSE, name = cat_name)
    } else if (discovery_score) {
      colorscale <- c("LOF" = "#e49444", "Not significant" = "lightgray", "GOF" = "#5778a4", "4" = "#967662", "5" = "#a87c9f", "6" = "#d62727")
      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]), 
        y = factor(.data[[mut.col]], levels = aa.order), fill = .data[[score.col]])) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        scale_fill_manual(values=colorscale, drop=FALSE, name = cat_name)
    } else if (compare_score) {
      # colorscale <- c("both" = "darkgray", "Lilace" = "navy", "original approach" = "maroon", "neither" = "lightgray", "5" = "#a87c9f", "6" = "#d62727")
      colorscale <- c("agree" = "darkgray", "+GOF" = "blue", "-GOF" = "navy", "-LOF" = "maroon", "+LOF" = "red", "=LOF" = "pink", "=GOF" = "lightblue")

      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]), 
        y = factor(.data[[mut.col]], levels = aa.order), fill = .data[[score.col]])) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        scale_fill_manual(values=colorscale, drop=FALSE, name = cat_name)
    } else if (category_score) {
      # colorscale <- c("" = "#5778a4", "abundance" = "#e49444", "abundance, surface_abundance " = "#6a9f58", "4" = "#967662", "5" = "#a87c9f", "6" = "#d62727")
      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]), 
        y = factor(.data[[mut.col]], levels = aa.order), fill = as.factor(.data[[score.col]]))) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) + 
        scale_fill_brewer(palette = "Set3")
        # geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        # scale_fill_manual(values=colorscale, drop=FALSE, name = cat_name)
    }
    else {
      p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]), 
        y = factor(.data[[mut.col]], levels = aa.order), fill = .data[[score.col]])) +
        geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
        do.call(scale_fill_continuous_divergingx, c.args)
    }
    p <- p + 
    # p <- ggplot(data = sub_data, aes(x = factor(.data[[pos.col]]), 
    #   y = factor(.data[[mut.col]], levels = aa.order), fill = .data[[score.col]])) +
    #   geom_tile(aes(color = factor(.data[[type.col]] == ctrl.name)), linewidth = 0.2) +
    #   do.call(scale_fill_continuous_divergingx, c.args) +  
      scale_color_manual(values = c(NA,'lightgreen'), name = ctrl.name, guide = "none") +
      # scale_x_discrete(breaks = sub_pos.order[seq(1, length(sub_pos.order), by = pos.step)]) +
      scale_x_discrete(breaks = sub_pos.order[seq(1, length(sub_pos.order), by = pos.step)]) +
      coord_fixed(ratio = 1, clip = "off") +
      theme_cowplot() +
      theme(
        plot.margin = unit(c(1, 0.25, 0.25, 0.25), "lines"),
        axis.title.y = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(angle = 45, hjust = 1, size = x.text),
        axis.text.y = element_text(size = y.text),
        axis.ticks = element_blank(),
        line = element_blank()
        # legend.position= legend,
        # legend.box = "horizontal"
      ) +
      geom_text(data = sub_data, 
                aes(x = factor(.data[[pos.col]], levels=sort(pos.order)), label = factor(.data[[wt.col]]), y = Inf), 
                vjust = -1, check_overlap = TRUE, size = seq.text) +
      labs(y = "Mutation", x = "Position")
    
    
    return(p)
  })
  
  p_all <- do.call(ggarrange, c(plot_list, list(nrow = nrow, ncol = ncol, common.legend=T, legend="right")))
  
  # save plots
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }
  
  ggsave(file.path(savedir, paste0(name,".png")), plot = p_all, height = ht, width = wd)
  if (savepdf) {
    ggsave(file.path(savedir, paste0(name,".pdf")), plot = p_all, height = ht, width = wd)
  } 
  if (savesvg) {
    ggsave(file.path(savedir, paste0(name,".svg")), plot = p_all, height = ht, width = wd)
  }
  
  if (show) {
    cat(paste("Showing the first", as.character(npos) , "positions.",
    "Full figure can be found in the saved directory."))
    print(plot_list[[1]])
  }
  return(p_all)
}


scoreDensity <- function(data,
                         savedir,
                         type.col = "type",
                         score.col = "mean",
                         hist = FALSE,
                         nbins = 30,
                         c.fill = c('#FF7575', 'lightgreen', "#7298BF"),
                         alpha = 0.5,
                         x.text = 10,
                         y.text = 10,
                         scale.free = FALSE,
                         space.free = FALSE,
                         ht = 10,
                         wd = 8,
                         name = "DensityPlot",
                         savepdf = TRUE,
                         savesvg = TRUE,
                         show = TRUE
){
  
  sc <- ifelse(scale.free, "free_y", "fixed")
  sp <- ifelse(space.free, "free_y", "fixed")
  if (length(c.fill) != length(levels(factor(data[[type.col]])))) {
    c.fill = c('#FF7575', "orange", 'lightgreen', "#7298BF", "#8828a8")
    warning("Length of color vector does not match the number of mutation types.")
  }
  
  p <- ggplot(data, aes(x = .data[[score.col]], fill = .data[[type.col]])) +
    facet_grid(.data[[type.col]] ~ ., scales = sc, space = sp) +
    theme_minimal() +
    scale_fill_manual(values = c.fill) +
    theme(
      strip.text.y = element_blank(),  
      strip.background = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      axis.text.x = element_text(size = x.text),
      axis.text.y = element_text(size = y.text)
    )
  
  if (hist) {
    p <- p + geom_histogram(aes(alpha = 0.5), bins = nbins, position="dodge") + 
      labs(x = "Score", y = "Count") +
      scale_alpha(guide = "none")
  }
  else{
    p <- p +geom_density(alpha = alpha) + labs(x = "Score", y = "Density")
  }
  
  # save plots
  if (!dir.exists(savedir)) {
    dir.create(savedir, recursive = TRUE)
  }
  
  ggsave(file.path(savedir, paste0(name,".png")), plot = p, height = ht, width = wd)
  if (savepdf) {
    ggsave(file.path(savedir, paste0(name,".pdf")), plot = p, height = ht, width = wd)
  } 
  if (savesvg) {
    ggsave(file.path(savedir, paste0(name,".svg")), plot = p, height = ht, width = wd)
  } 
  
  if (show) {
    print(p)
  }
}

plot_method_scatter <- function(df_called, method1, method2, m1_name, m2_name, contour=F, title="") {
  effect_col <- rep(NA, nrow(df_called))
  call1 <- df_called[[paste0(method1, "_disc")]]
  call2 <- df_called[[paste0(method2, "_disc")]]
  effect_col[call1 & call2] <- "both"
  effect_col[call1 & !call2] <- method1
  effect_col[!call1 & call2] <- method2
  effect_col[!call1 & !call2] <- "neither"
  print(table(effect_col, useNA="ifany"))
  plot_df <- df_called
  plot_df$effect_col <- effect_col
  plot_df$effect_col[plot_df$effect_col=="effect"] <- "Lilace"
  plot_df$effect_col[plot_df$effect_col=="weight_mean"] <- "Mean bin"
  plot_df$effect_col[plot_df$effect_col=="enrich"] <- "Enrich2"
  
  cat_counts <- plot_df[!is.na(plot_df$effect_col),] %>% group_by(effect_col) %>% summarise(count = n())
  cat_counts <- cat_counts %>% arrange(desc(effect_col))
  labels <- paste0(cat_counts$effect_col, " (", cat_counts$count, ")")
  # print(labels)
  plot_df$effect_col <- factor(plot_df$effect_col, levels=cat_counts$effect_col)
  plot_df <- plot_df[order(plot_df$effect_col),]
  

  call_colors <- setNames(c('#4e79a7', '#76b7b2', '#f28e2b', '#f28e2b', 
                          '#e15759', '#e15759', '#e15759', '#e15759', '#e15759', '#e15759'),
                          c("both", "neither", "Lilace", "effect_sep",
                          "effect_joint", "Mean bin", "weight_effect_mean_syn_sd", "ML_effect_mean", "ML_effect_mean_syn_sd", "Enrich2"))
  
  p <- ggplot(plot_df[!is.na(plot_df$effect_col),] %>% arrange(effect_col), aes(x=get(method1), y=get(method2), col=as.factor(effect_col))) + geom_point(alpha=0.3, size=2) +
    labs(x=m1_name, y=m2_name) +
    # scale_color_discrete(labels = labels) +
    scale_color_manual(labels = labels, values=call_colors) +
    # scale_color_tableau(labels=labels) +
    theme_cowplot() + 
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    guides(color=guide_legend(title="Discovery calls", override.aes = list(alpha=1)))
      
  dens1 <- ggplot(plot_df[!is.na(plot_df$effect_col),], aes(x = get(method1), fill = as.factor(effect_col))) + 
    geom_density(alpha = 0.4) + 
    theme_void() +
    # scale_fill_discrete(labels = labels) +
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    scale_fill_manual(labels = labels, values=call_colors) +
    guides(fill="none")
  
  dens2 <- ggplot(plot_df[!is.na(plot_df$effect_col),], aes(x = get(method2), fill = as.factor(effect_col))) + 
    geom_density(alpha = 0.4) + 
    theme(text = element_text(size = 20), axis.text = element_text(size = 20)) +
    theme_void() +
    coord_flip() +
    # scale_fill_discrete(labels = labels) +
    scale_fill_manual(labels = labels, values=call_colors) +
    guides(fill="none")
  
  # p_full <- dens1 + plot_spacer() + p + dens2 + 
  #   plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4), guides = "collect") +
  #   plot_annotation(paste(m1_name, "vs", m2_name, "effects"),theme=theme(plot.title=element_text(hjust=0.5)))
  # ggsave(paste0(plot_dir, "/", method1, "_", method2, ".png"), p_full)

  dens11 <- ggplot(plot_df[!is.na(plot_df$effect_col),], aes(x = get(method1), fill = type)) + 
    geom_density(alpha = 0.4) + 
    scale_color_manual(values = c(deletion="#b07aa1", insertion="#326191", missense="#f28e2b", synonymous="#59a14f"))+
    scale_fill_manual(values = c(deletion="#b07aa1", insertion="#326191", missense="#f28e2b", synonymous="#59a14f")) +
    theme_void() + 
    theme(text = element_text(size = 20), axis.text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
    guides(fill=guide_legend(title="Mutation types"))
  
  dens22 <- ggplot(plot_df[!is.na(plot_df$effect_col),], aes(x = get(method2), fill = type)) + 
    geom_density(alpha = 0.4) + 
    theme_void() + 
    theme(text = element_text(size = 20), axis.text = element_text(size = 20), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
    scale_color_manual(values = c(deletion="#b07aa1", insertion="#326191", missense="#f28e2b", synonymous="#59a14f"))+
    scale_fill_manual(values = c(deletion="#b07aa1", insertion="#326191", missense="#f28e2b", synonymous="#59a14f")) +
    coord_flip() +
    guides(fill=guide_legend(title="Mutation types"))

  if (stringr::str_length(title) == 0) {
    title <- paste(m2_name, "vs", m1_name, "effects")
  }
  
  # p <- dens11 + plot_spacer() + plot_spacer() + dens1 + plot_spacer() + plot_spacer() + p + dens2 + dens22 +
  #   plot_layout(ncol = 3, nrow = 3, widths = c(4, 1, 1), heights = c(1, 1, 4), guides = "collect") +
  #   plot_annotation(title,theme=theme(plot.title=element_text(hjust=0.5)))

  p <- dens11 + plot_spacer() + p + dens22 +
    plot_annotation(title, theme=theme(plot.title=element_text(hjust=0.5, size=20), legend.text=element_text(size=20))) +
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4), guides = "collect")
    # ggsave(paste0(plot_dir, "/", method1, "_", method2, ".png"), p)
 
  if (contour) {
    p <- ggplot(mapping=aes(x=get(method1), y=get(method2), col=as.factor(effect_col))) + 
    lapply(split(plot_df, plot_df$effect_col), function(x) stat_density_2d(data=x, alpha=0.8)) +
    theme_cowplot() + labs(x=m1_name, y=m2_name) +
    # scale_color_discrete(labels = labels) +
    scale_color_manual(labels = labels, values=call_colors) +
    guides(color=guide_legend(title="Discovery calls"))

    p <- dens11 + plot_spacer() + p + dens22 +
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4), guides = "collect")
    # plot_annotation(title,theme=theme(plot.title=element_text(hjust=0.5)))
  }
  return(p)
}


make_cor_plots <- function(yaml, cor_df, outdir) {
  params <- unique(cor_df$param)
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(params)) {
        for (rescal in c(T)) {
          cur_param <- params[i]
          cur_cor_df <- cor_df[cor_df$dataset == dataset & cor_df$n_pos == n_pos & cor_df$param==cur_param & cor_df$is_rescaled==rescal,]
          cor_df_summary <- cur_cor_df %>% group_by(param_val, method) %>% summarize(mean_corr = mean(corr), sd_corr = sd(corr))
          cor_df_summary <- cor_df_summary[grepl("Lilace|ML\n2sd|mean\ bin\n2sd|enrich2|approx", cor_df_summary$method),]
          cor_df_summary$sample_ord <- factor(cor_df_summary$method, levels = c("Lilace", "Lilace\n(unrecalibrated)", "mean bin\n2sd", "mean bin\nnorm approx", "ML\n2sd", "ML\nnorm approx", "enrich2\n2sd", "enrich2\np-value"))
          if (cur_param == "default") {
             p_corr <- ggplot(cor_df_summary, aes(x=sample_ord, y=mean_corr)) + 
              geom_pointrange(aes(ymin = mean_corr - sd_corr, ymax = mean_corr + sd_corr), color="orange", size=1) +
              # geom_boxplot(width = 0.5, position = position_dodge(0.6), fill="orange") + 
              expand_limits(y=0) + 
              theme_cowplot() + xlab("method") + ylab("True Effect Corr") + guides(color=guide_legend(title=cur_param)) +
              scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
          } else {
            p_corr <- ggplot(cor_df_summary, aes(x=sample_ord, y=mean_corr, color=as.factor(param_val))) + 
              expand_limits(y=0) + 
              geom_pointrange(aes(ymin = mean_corr - sd_corr, ymax = mean_corr + sd_corr), position = position_dodge2(0.5, preserve = "single"), size=0.75) +
              # geom_boxplot(width = 0.5, position = position_dodge(0.6)) + 
              theme_cowplot() + xlab("method") + ylab("True Effect Corr") + guides(color=guide_legend(title=cur_param)) +
              # scale_colour_ptol() +
              scale_color_manual(values=c("#E69F00", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#505050")) +
              scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
          }
          p <- p_corr + plot_layout(ncol=1, guides = "collect", widths = c(3, 3, 1)) & theme(text = element_text(size = 20), axis.text = element_text(size = 20))
          ggsave(paste0(outdir, "/", dataset, n_pos, cur_param, rescal, "_boxes.png"), p, height=3, width=13)
        }
      }
    }
  }
}

# make cross param FDR plots
make_FDR_plots <- function(yaml, FDR_df, outdir, plot_mean_effect=T, axis_text_size=16, filter_methods=T) {
  params <- unique(FDR_df$param)
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(params)) {
        for (rescal in c(T)) {
          cur_param <- params[i]
          cur_FDR_df <- FDR_df[FDR_df$dataset == dataset & FDR_df$n_pos == n_pos & FDR_df$param==cur_param & FDR_df$is_rescaled==rescal,]
          if (!plot_mean_effect) {
            cur_FDR_df <- cur_FDR_df %>% filter(!stringr::str_detect(method, "syn_sd"))
          }
          FDR_df_summary <- cur_FDR_df %>% group_by(param_val, method) %>% summarize(mean_FDR = mean(FDR), mean_sens = mean(1-FNR), mean_FPR = mean(FPR),
                                                                                  sd_FDR = sd(FDR), sd_sens = sd(1-FNR), sd_FPR = sd(FPR))
          
          if (filter_methods) {
            FDR_df_summary <- FDR_df_summary[grepl("Lilace|ML\n2sd|mean\ bin\n2sd|enrich2|approx", FDR_df_summary$method),]
            FDR_df_summary$sample_ord <- factor(FDR_df_summary$method, levels = c("Lilace", "Lilace\n(unrecalibrated)", "mean bin\n2sd", "mean bin\nnorm approx", "ML\n2sd", "ML\nnorm approx", "enrich2\n2sd", "enrich2\np-value"))
          } else {
            FDR_df_summary$sample_ord <- factor(FDR_df_summary$method)
          }
          # FDR_df_summary <- FDR_df_summary[grepl("Lilace|ML\n2sd|mean\ bin\n2sd|enrich2", FDR_df_summary$method),]
          # FDR_df_summary$sample_ord <- factor(FDR_df_summary$method, levels = c("Lilace", "Lilace\n(unrecalibrated)", "ML\n2sd", "mean bin\n2sd", "enrich2\n2sd", "enrich2\np-value"))
          if (cur_param == "default") {
             p_FDR <- ggplot(FDR_df_summary, aes(x=sample_ord, y=mean_FDR)) + 
              geom_pointrange(aes(ymin = mean_FDR - sd_FDR, ymax = mean_FDR + sd_FDR), color="orange", size=1) +
              # geom_boxplot(width = 0.5, position = position_dodge(0.6), fill="orange") + 
              expand_limits(y=0) + 
              theme_cowplot() + xlab("method") + ylab("FDR") + guides(color=guide_legend(title=cur_param)) +
              scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
              # theme(axis.text.x = element_text(angle = 70, vjust = 0.3, hjust=1), legend.position = 'bottom')

             p_FNR <- ggplot(FDR_df_summary, aes(x=sample_ord, y=mean_sens)) + 
             geom_pointrange(aes(ymin = mean_sens - sd_sens, ymax = mean_sens + sd_sens), color="orange", size=1) +
              # geom_boxplot(width = 0.5, position = position_dodge(0.6), fill="orange") + 
              expand_limits(y=0) + 
              theme_cowplot() + xlab("method") + ylab("Sensitivity") + guides(color=guide_legend(title=cur_param)) +
              scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
              # theme(axis.text.x = element_text(angle = 70, vjust = 0.3, hjust=1), legend.position = 'bottom')

            p_FPR <- ggplot(FDR_df_summary, aes(x=sample_ord, y=mean_FPR)) + 
              geom_pointrange(aes(mean_FPR - sd_FPR, ymax = mean_FPR + sd_FPR), color="orange", size=1) +
              # geom_boxplot(width = 0.5, position = position_dodge(0.6), fill="orange") + 
              expand_limits(y=0) + 
              theme_cowplot() + xlab("method") + ylab("FPR") + guides(color=guide_legend(title=cur_param)) +
              scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
              # theme(axis.text.x = element_text(angle = 70, vjust = 0.3, hjust=1), legend.position = 'bottom')
          } else {
            p_FDR <- ggplot(FDR_df_summary, aes(x=sample_ord, y=mean_FDR, color=as.factor(param_val))) + 
              expand_limits(y=0) + 
              geom_pointrange(aes(ymin = mean_FDR - sd_FDR, ymax = mean_FDR + sd_FDR), position = position_dodge2(0.5, preserve = "single"), size=0.75) +
              # geom_boxplot(width = 0.5, position = position_dodge(0.6)) + 
              theme_cowplot() + xlab("method") + ylab("FDR") + guides(color=guide_legend(title=cur_param)) +
              # scale_colour_ptol() +
              scale_color_manual(values=c("#E69F00", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#505050")) +
              scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
              # theme(axis.text.x = element_text(angle = 70, vjust = 0.3, hjust=1), legend.position = 'bottom')

            p_FNR <- ggplot(FDR_df_summary, aes(x=sample_ord, y=mean_sens, color=as.factor(param_val))) + 
              expand_limits(y=0) + 
              geom_pointrange(aes(ymin = mean_sens - sd_sens, ymax = mean_sens + sd_sens), position = position_dodge2(0.5, preserve = "single"), size=0.75) +
              # geom_boxplot(width = 0.5, position = position_dodge(0.6)) + 
              theme_cowplot() + xlab("method") + ylab("Sensitivity") + guides(color=guide_legend(title=cur_param)) +
              scale_color_manual(values=c("#E69F00", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#505050")) +
              scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
              # theme(axis.text.x = element_text(angle = 70, vjust = 0.3, hjust=1), legend.position = 'bottom')

            p_FPR <- ggplot(FDR_df_summary, aes(x=sample_ord, y=mean_FPR, color=as.factor(param_val))) + 
              expand_limits(y=0) + 
              geom_pointrange(aes(mean_FPR - sd_FPR, ymax = mean_FPR + sd_FPR), position = position_dodge2(0.5, preserve = "single"), size=0.75) +
              # geom_boxplot(width = 0.5, position = position_dodge(0.6)) + 
              theme_cowplot() + xlab("method") + ylab("FPR") + guides(color=guide_legend(title=cur_param)) +
              scale_color_manual(values=c("#E69F00", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442", "#505050")) +
              scale_y_continuous(labels = scales::label_number(accuracy = 0.01))
              # theme(axis.text.x = element_text(angle = 70, vjust = 0.3, hjust=1), legend.position = 'bottom')
          }

          # p <- p_FDR + p_FNR + p_FPR + plot_layout(ncol=3, guides = "collect") + plot_annotation(title=paste(cur_param, ": ", dataset, n_pos, "pos", "rescaled:", rescal)) & theme(legend.position = 'bottom')
          p <- p_FDR + p_FNR + guide_area() + plot_layout(ncol=3, guides = "collect", widths = c(3, 3, 1)) & theme(text = element_text(size = 20), axis.text = element_text(size = axis_text_size))
          ggsave(paste0(outdir, "/", dataset, n_pos, cur_param, rescal, "_boxes.png"), p, height=3, width=13)

        }
      }
    }
  }
}


# plot cross sig metrics (FDR-FNR, precision-recall, FDR-sens)
make_cross_sig_plots <- function(yaml, sig_FDR_df, outdir, plot_mean_effect=T) {
  param <- unique(sig_FDR_df$param)
  param_val <- unique(sig_FDR_df$param_val)
  rescal_vals <- unique(sig_FDR_df$is_rescaled)
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (rescal in rescal_vals) {
        cur_FDR_FNR_df <- sig_FDR_df[sig_FDR_df$dataset == dataset & sig_FDR_df$n_pos == n_pos & 
          sig_FDR_df$param==param & sig_FDR_df$param_val==param_val & sig_FDR_df$is_rescaled==rescal,]
        if (!plot_mean_effect) {
          cur_FDR_FNR_df <- cur_FDR_FNR_df %>% filter(!stringr::str_detect(method, "syn_sd"))
        }
        
        FDR_FNR_df_summary <- cur_FDR_FNR_df %>% group_by(method, sig_cutoff) %>% summarize(mean_FDR = mean(FDR), mean_FNR = mean(FNR), mean_FPR = mean(FPR), sd_FDR = sd(FDR), sd_FNR = sd(FNR), sd_FPR = sd(FPR), mean_sens = mean(sens), sd_sens=sd(sens))
        # FDR-FNR plot
        p_FDR_FNR <- ggplot(FDR_FNR_df_summary, aes(x=mean_FDR, y=mean_FNR, col=method)) + geom_point() + geom_line(lwd=1) + 
          ggtitle(paste("FDR-FNR: ", dataset, n_pos, "pos", "rescaled:", rescal)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_cowplot() + scale_color_tableau()
        print(p_FDR_FNR)
        ggsave(paste0(outdir, "/", dataset, n_pos, param, param_val, rescal, "FDR-FNR", ".png"), p_FDR_FNR)
        
        # precision -recall (1-FDR and sens)
        p_prec_recall <- ggplot(FDR_FNR_df_summary, aes(x=1-mean_FDR, y=mean_sens, col=method)) + geom_point() + geom_line(lwd=1) + 
          ggtitle(paste("prec-recall: ", dataset, n_pos, "pos", "rescaled:", rescal)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_cowplot() + scale_color_tableau()
        print(p_prec_recall)
        ggsave(paste0(outdir, "/", dataset, n_pos, param, param_val, rescal, "FDR-recall", ".png"), p_prec_recall)
        
        # FDR-sensitivity plot
        p_FDR_sens <- ggplot(FDR_FNR_df_summary, aes(x=mean_FDR, y=mean_sens, col=method)) + geom_point() + geom_line(lwd=1) + 
          # ggtitle(paste("FDR-sensitivity: ", dataset)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_cowplot() + scale_color_tableau() + theme(text = element_text(size = 14), axis.text = element_text(size = 12)) +
          xlab("Mean FDR") + ylab("Mean Sensitivity")
        print(p_FDR_sens)
        
        ggsave(paste0(outdir, "/", dataset, n_pos, param, param_val, rescal, "FDR-sens", ".png"), p_FDR_sens, height=6, width=11)
      }
    }
  }
}

# plot metrics by rank
make_rank_FDR_plots <- function(yaml, rank_FDR_df, outdir, plot_mean_effect=T) {
  param <- unique(rank_FDR_df$param)
  param_val <- unique(rank_FDR_df$param_val)
  rescal_vals <- unique(rank_FDR_df$is_rescaled)
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (rescal in rescal_vals) {
        cur_rank_df <- rank_FDR_df[rank_FDR_df$dataset==dataset & rank_FDR_df$n_pos == n_pos & rank_FDR_df$is_rescaled == rescal &
          rank_FDR_df$param == param & rank_FDR_df$param_val == param_val,]
        if (!plot_mean_effect) {
          cur_rank_df <- cur_rank_df %>% filter(!stringr::str_detect(method, "syn_sd"))
        }
        cur_rank_df_sum <- cur_rank_df %>% group_by(method, rank_cutoff) %>% 
          summarize(mean_FDR = mean(FDR), mean_FNR = mean(FNR), mean_FPR = mean(FPR),
                    sd_FDR = sd(FDR), sd_FNR = sd(FNR), sd_FPR = sd(FPR))

        p <- ggplot(cur_rank_df_sum, aes(x=rank_cutoff, y=mean_FDR, col=method)) + geom_point() + geom_line(lwd=1) + 
          ggtitle(paste("FDR: ", dataset, n_pos, "pos", "rescaled:", rescal)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_bw()

        ggsave(paste0(outdir, "/", dataset, n_pos, param, param_val, rescal, "rank_FDR", ".png"), p)
        
        p <- ggplot(cur_rank_df_sum, aes(x=rank_cutoff, y=mean_FNR, col=method)) + geom_point() + geom_line(lwd=1) + 
          # geom_ribbon(aes(ymin = mean_FNR - 2*sd_FNR, ymax = mean_FNR + 2*sd_FNR, color=method, fill=method), alpha=0.1) +
          ggtitle(paste("FNR: ", dataset, n_pos, "pos", "rescaled:", rescal)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_bw()
        ggsave(paste0(outdir, "/", dataset, n_pos, param, param_val, rescal, "rank_FNR", ".png"), p)
        
        p <- ggplot(cur_rank_df_sum, aes(x=rank_cutoff, y=mean_FPR, col=method)) + geom_point() + geom_line(lwd=1) + 
          # geom_ribbon(aes(ymin = mean_FPR - 2*sd_FPR, ymax = mean_FPR + 2*sd_FPR, color=method, fill=method), alpha=0.1) +
          ggtitle(paste("FPR: ", dataset, n_pos, "pos", "rescaled:", rescal)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_bw()
        ggsave(paste0(outdir, "/", dataset, n_pos, param, param_val, rescal, "rank_FPR", ".png"), p)
      }
    }
  }
}


# plot discovery boxplots
make_disc_FDR_boxplots <- function(yaml, disc_FDR_df, n_discoveries, outdir, plot_mean_effect=T) {
  param <- unique(disc_FDR_df$param)
  param_val <- unique(disc_FDR_df$param_val)
  rescal_vals <- unique(disc_FDR_df$is_rescaled)
  for (dataset in yaml$datasets) {
    for (n_pos in yaml$n_pos_vals) {
      for (i in 1:length(yaml$params)) {
        for (rescal in rescal_vals) {
          cur_FDR_FNR_df <- disc_FDR_df[disc_FDR_df$dataset == dataset & disc_FDR_df$n_pos == n_pos & 
            disc_FDR_df$param==param & disc_FDR_df$param_val==param_val & disc_FDR_df$is_rescaled==rescal,]
          if (!plot_mean_effect) {
            cur_FDR_FNR_df <- cur_FDR_FNR_df %>% filter(!stringr::str_detect(method, "syn_sd"))
          }
          FDR_FNR_df_summary <- cur_FDR_FNR_df %>% group_by(method) %>%
            summarize(mean_FDR = mean(FDR), mean_FNR = mean(FNR), mean_FPR = mean(FPR),
                      sd_FDR = sd(FDR), sd_FNR = sd(FNR), sd_FPR = sd(FPR),
                      mean_sens = mean(sens), sd_sens=sd(sens))
          cur_FDR_FNR_df <- merge(cur_FDR_FNR_df, FDR_FNR_df_summary)
          p <- ggplot(cur_FDR_FNR_df, aes(x=mean_FDR, y=sens, col=method)) + geom_boxplot(width=0.0075, position = position_dodge(preserve = "single")) + theme_cowplot() + xlab(paste0("mean FDR at ", n_discoveries, " discoveries")) + ylab("sensitivity") + ggtitle(paste(dataset, n_pos, "positions"))
          ggsave(paste0(outdir, "/", dataset, n_pos, param, param_val, rescal, "n_discoveries", n_discoveries, "_full.png"), p)
        }
      }
    }
  }
}

make_masked_FDR_plots <- function(rank_FDR_df, dataset, outdir, plot_mean_effect=T) {
  cur_rank_df <- rank_FDR_df
  if (!plot_mean_effect) {
    cur_rank_df <- cur_rank_df %>% filter(!stringr::str_detect(method, "syn_sd"))
  }
  cur_rank_df_sum <- cur_rank_df %>% group_by(method, sig_cutoff) %>% 
    summarize(N_syn_called = mean(N_syn_called))
  
  p <- ggplot(cur_rank_df_sum, aes(x=sig_cutoff, y=N_syn_called, col=method)) + geom_point() + geom_line(lwd=1) +
    ggtitle(paste("Proportion masked synonymous called incorrectly by significance: ", dataset)) +
    scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
    "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
    theme_bw()
  ggsave(paste0(outdir, "/", dataset, "_masked_syn", ".png"), p)
}

make_ranked_masked_FDR_plots <- function(rank_FDR_df, dataset, outdir, plot_mean_effect=T) {
  cur_rank_df <- rank_FDR_df
  if (!plot_mean_effect) {
    cur_rank_df <- cur_rank_df %>% filter(!stringr::str_detect(method, "syn_sd"))
  }
  cur_rank_df_sum <- cur_rank_df %>% group_by(method, rank_cutoff) %>% 
    summarize(N_syn_called = mean(N_syn_called))
  
  p <- ggplot(cur_rank_df_sum, aes(x=rank_cutoff, y=N_syn_called, col=method)) + geom_point() + geom_line(lwd=1) +
    ggtitle(paste("Proportion masked synonymous called incorrectly by rank: ", dataset)) +
    scale_color_manual(values=c("#505050", "#E69F00", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#F0E442")) +
    theme_bw()
  ggsave(paste0(outdir, "/", dataset, "_masked_syn_ranked", ".png"), p)
}

make_ranked_nonsense_FNR_plots <- function(rank_FDR_df, dataset, outdir, plot_mean_effect=T) {
  cur_rank_df <- rank_FDR_df
  if (!plot_mean_effect) {
    cur_rank_df <- cur_rank_df %>% filter(!stringr::str_detect(method, "syn_sd"))
  }
  cur_rank_df_sum <- cur_rank_df %>% group_by(method, rank_cutoff) %>% 
    summarize(N_nonsense_called = mean(N_nonsense_called))
  
  p <- ggplot(cur_rank_df_sum, aes(x=rank_cutoff, y=N_nonsense_called, col=method)) + geom_point() + geom_line(lwd=1) +
    ggtitle(paste("Proportion nonsense called correctly by rank: ", dataset)) +
    scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
    "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
    theme_bw()
  ggsave(paste0(outdir, "/", dataset, "_nonsense_ranked", ".png"), p)
}

make_am_FNR_plots <- function(rank_FDR_df, dataset, outdir, plot_mean_effect=T) {
  cur_rank_df <- rank_FDR_df
  if (!plot_mean_effect) {
    cur_rank_df <- cur_rank_df %>% filter(!stringr::str_detect(method, "syn_sd"))
  }
  cur_rank_df_sum <- cur_rank_df %>% group_by(method, sig_cutoff) %>% 
    summarize(mean_FDR = mean(FDR), mean_FNR = mean(FNR), mean_FPR = mean(FPR),
              sd_FDR = sd(FDR), sd_FNR = sd(FNR), sd_FPR = sd(FPR))
  
  p <- ggplot(cur_rank_df_sum, aes(x=sig_cutoff, y=mean_FNR, col=method)) + geom_point() + geom_line(lwd=1) + 
          geom_ribbon(aes(ymin = mean_FNR - 2*sd_FNR, ymax = mean_FNR + 2*sd_FNR, color=method, fill=method), alpha=0.1) +
          ggtitle(paste("alphamissense-based FNR: ", dataset)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_bw()
  ggsave(paste0(outdir, "/", dataset, "_am_FNR", ".png"), p)
}

make_nonsense_FNR_plots <- function(rank_FDR_df, dataset, outdir, plot_mean_effect=T) {
  cur_rank_df <- rank_FDR_df
  if (!plot_mean_effect) {
    cur_rank_df <- cur_rank_df %>% filter(!stringr::str_detect(method, "syn_sd"))
  }
  cur_rank_df_sum <- cur_rank_df %>% group_by(method, sig_cutoff) %>% 
    summarize(mean_FDR = mean(FDR), mean_FNR = mean(FNR), mean_FPR = mean(FPR),
              sd_FDR = sd(FDR), sd_FNR = sd(FNR), sd_FPR = sd(FPR))
  
  p <- ggplot(cur_rank_df_sum, aes(x=sig_cutoff, y=mean_FNR, col=method)) + geom_point() + geom_line(lwd=1) + 
          geom_ribbon(aes(ymin = mean_FNR - 2*sd_FNR, ymax = mean_FNR + 2*sd_FNR, color=method, fill=method), alpha=0.1) +
          ggtitle(paste("nonsense-based FNR: ", dataset)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_bw()
  ggsave(paste0(outdir, "/", dataset, "_nonsense_FNR", ".png"), p)
}

make_clinvar_FNR_plots <- function(rank_FDR_df, dataset, outdir, plot_mean_effect=T) {
  cur_rank_df <- rank_FDR_df
  if (!plot_mean_effect) {
    cur_rank_df <- cur_rank_df %>% filter(!stringr::str_detect(method, "syn_sd"))
  }
  cur_rank_df_sum <- cur_rank_df %>% group_by(method, sig_cutoff) %>% 
    summarize(mean_FDR = mean(FDR), mean_FNR = mean(FNR), mean_FPR = mean(FPR),
              sd_FDR = sd(FDR), sd_FNR = sd(FNR), sd_FPR = sd(FPR))
  
  p <- ggplot(cur_rank_df_sum, aes(x=sig_cutoff, y=mean_FNR, col=method)) + geom_point() + geom_line(lwd=1) + 
          geom_ribbon(aes(ymin = mean_FNR - 2*sd_FNR, ymax = mean_FNR + 2*sd_FNR, color=method, fill=method), alpha=0.1) +
          ggtitle(paste("clinvar-based FNR: ", dataset)) +
          scale_color_manual(values=c("#967662", "#85b6b2",  "#5778a4", "#a87c9f", "#6a9f58", 
          "#d1615d", "#f1a2a9", "#b8b0ac", "#fc7d0b", "#3a4032", "#47817f", "#ed5e93", "#adbe77", "#cda211")) +
          theme_bw()
  ggsave(paste0(outdir, "/", dataset, "_clinvar_FNR", ".png"), p)
}


make_3d_plot <- function(df, plot_dir, modelname="FACS_double_sample_syn_recalibrated") {
  print("plotting 3D")
  simplex <- function(n) {
    qr.Q(qr(matrix(1, nrow=n)) ,complete = TRUE)[,-1]
  }
  tetra <- simplex(4)
  counts <- df %>% ungroup() %>% select(c("c_0", "c_1", "c_2", "c_3"))
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
}