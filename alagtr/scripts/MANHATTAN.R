
library(tidyverse)


manhat_plot <- function(mod, outliers) {
    # Make and get tidy data frames for plotting
    #snp_scores <- vegan::scores(mod, choices = 1:ncol(mod$CCA$v), display = "species", scaling = "none")
    # TAB_snps <- data.frame(names = row.names(snp_scores), snp_scores)
    # TAB_snps$type <- "Non-outlier"
    # TAB_snps$type[TAB_snps$label %in% rdadapt$locus] <- "Outlier"
    # TAB_snps$type <- factor(TAB_snps$type, levels = c("Non-outlier", "Outlier"))

    TAB_manhattan <- data.frame(
        pos = 1:nrow(rdadapt),
        pvalues = rdadapt$p.values,
        #chrom = rdadapt$scaff,
        locus = rdadapt$locus)

    TAB_manhattan <- TAB_manhattan %>% 
        separate(locus, sep = "_", c("chrom", "chrom_pos", "ref", "alt"))

    TAB_manhattan <- TAB_manhattan[order(TAB_manhattan$pos), ]
    TAB_manhattan$chrom <- factor(TAB_manhattan$chrom, levels = (unique(TAB_manhattan$chrom)))

    ylim <- TAB_manhattan %>% 
        dplyr::filter(pvalues == min(pvalues)) %>% 
        mutate(ylim = abs(floor(log10(pvalues))) + 2) %>% 
        pull(ylim)

    axis_set <- TAB_manhattan %>% 
        group_by(chrom) %>% 
        summarize(center = mean(pos))

    TAB_manhattan <- TAB_manhattan %>% 
        mutate(type = ifelse(-log10(pvalues) > sig, "Outlier", "Non-outlier"))

    plt_manhat <-
        ggplot2::ggplot() +
        ggplot2::geom_point(data = TAB_manhattan %>% dplyr::filter(type == "Outlier"), 
                            ggplot2::aes(x = pos, y = -log10(pvalues),), col = "orange", size = 1.4, alpha = 0.75) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab("-log10(p)") +
        ggplot2::geom_hline(yintercept = -log10(sig), linetype = "dashed", color = "black", linewidth = 0.6) +
        scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
        ggplot2::geom_point(data = TAB_manhattan %>% dplyr::filter(type == "Non-outlier"), 
                            ggplot2::aes(x = pos, y = -log10(pvalues), col = chrom), size = 1.4, alpha = 0.75) +
        ggplot2::scale_color_manual(values = rep(c("#276FBF","#183059"), 
                                                ceiling(length(unique(TAB_manhattan$chrom))/2))[1:length(unique(TAB_manhattan$chrom))]) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            plot.background = ggplot2::element_blank(),
            axis.text.x = element_text(angle = 60, size = 4, vjust = 0.5)
        )


    chromosomes <- unique(TAB_manhattan$chrom)
        
    pdf(snakemake@output[["per_chrom"]], height = 8.5, width = 11)
    
    for (chr in chromosomes) {
        
        plots_list <- list()
        
        df_chr <- TAB_manhattan %>% 
        filter(chrom == chr)
        
        
        chr_plot <- ggplot2::ggplot() +
                ggplot2::geom_point(data = df_chr, 
                                    ggplot2::aes(x = pos, y = -log10(pvalues),), size = 1.4, alpha = 0.75) +
                #ggplot2::geom_point(data = df_chr %>% dplyr::filter(type == "Outlier"), 
                #                    ggplot2::aes(x = pos, y = -log10(pvalues),), col = "orange", size = 1.4, alpha = 0.75) +
                ggplot2::xlab(NULL) +
                ggplot2::ylab("-log10(p)") +
                ggplot2::geom_hline(yintercept =sig, linetype = "dashed", color = "black", linewidth = 1) +
                scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
                scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
                #ggplot2::geom_point(data = df_chr %>% dplyr::filter(type == "Non-outlier"), 
                #                    ggplot2::aes(x = pos, y = -log10(pvalues), col = "black"), size = 1.4, alpha = 0.75) +
                ggplot2::theme_bw(base_size = 11) +
                ggplot2::theme(
                    legend.position = "none",
                    panel.grid = ggplot2::element_blank(),
                    plot.background = ggplot2::element_blank(),
                    axis.text.x = element_text(size = 16, vjust = 0.5)
                )
        
        print(chr_plot)
        #plots_list[[chr]] <- roh_plot
        
    }
    
    dev.off()

    return(plt_manhat)

}

rdadapt <- snakemake@input[["rdadapt_output"]]
rdadapt <- read_csv(rdadapt)
sig <- snakemake@params[["sig"]]
mantplot <- manhat_plot(rdadapt, sig)
ggsave(snakemake@output[["manhattan_plot"]], mantplot, width = 12, height = 5, dpi = 300, bg = "white")


#' Build Manhattan plots of LFMM results; modified from algatr's `lfmm_manhattanplot()` function
#'
#' @param results LFMM df output
#' @param sig significance threshold
#' @param color_by_contig whether to color by contig or not (defaults to FALSE)
#'
#' @return
#' @export
manhat_plot_lfmm <- function(results, sig = 0.05, color_by_contig = FALSE) {
  df <- data.frame(results$df)
  # if (!is.null(var)) df <- df[df$var %in% var, ]
  df$type[df$adjusted.pvalue < sig] <- "Outlier"
  df$type[!(df$adjusted.pvalue < sig)] <- "Non-outlier"
  df$index <- 1:length(unique(df$snp))

  if (color_by_contig == FALSE) {
    plt_manhat <-
      ggplot2::ggplot(df, ggplot2::aes(x = index, y = -log10(adjusted.pvalue))) +
      ggplot2::geom_point(alpha = 0.75, pch = 16, ggplot2::aes(col = type)) +
      ggplot2::scale_color_manual(values = c("Non-outlier" = rgb(0.7, 0.7, 0.7, 0.5), "Outlier" = "#F9A242FF"), na.translate = F) +
      ggplot2::xlab("SNPs") +
      ggplot2::ylab("-log10(p)") +
      ggplot2::guides(color = ggplot2::guide_legend(title = "SNP type")) +
      ggplot2::facet_wrap(~var, nrow = length(unique(df$var))) +
      ggplot2::xlab("Position") +
      ggplot2::ylab("-log10(p)") +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position = "right",
        legend.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = ggplot2::rel(.8)),
        strip.text = ggplot2::element_text(size = 11))
  }
  
  if (color_by_contig == TRUE) {
    TAB_manhattan <- TAB_manhattan %>% 
      tibble::rownames_to_column(var = "name") %>% 
      tidyr::separate_wider_delim(cols = name,
                                  delim = "_",
                                  names = c("chrom", "site"))
    
    TAB_manhattan <- TAB_manhattan[order(TAB_manhattan$pos), ]
    TAB_manhattan$chrom <- factor(TAB_manhattan$chrom, levels = (unique(TAB_manhattan$chrom)))
    
    axis_set <- TAB_manhattan %>% 
      group_by(chrom) %>% 
      summarize(center = mean(pos))
    
    ylim <- TAB_manhattan %>%
      dplyr::filter(pvalues == min(pvalues)) %>%
      mutate(ylim = abs(floor(log10(pvalues))) + 2) %>%
      pull(ylim)
    
    plt_manhat <-
      ggplot2::ggplot() +
      ggplot2::geom_point(data = df %>% dplyr::filter(type == "Outlier"), 
                          ggplot2::aes(x = index, y = -log10(adjusted.pvalue),), col = "orange", size = 1.4, alpha = 0.75) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab("-log10(p)") +
      ggplot2::geom_hline(yintercept = -log10(sig), linetype = "dashed", color = "black", linewidth = 0.6) +
      # scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center) +
      # scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
      ggplot2::geom_point(data = df %>% dplyr::filter(type == "Non-outlier"), 
                          ggplot2::aes(x = pos, y = -log10(adjusted.pvalue), col = chrom), size = 1.4, alpha = 0.75) +
      ggplot2::scale_color_manual(values = rep(c("#276FBF","#183059"), 
                                               ceiling(length(unique(TAB_manhattan$chrom))/2))[1:length(unique(TAB_manhattan$chrom))]) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        axis.text.x = element_text(angle = 60, size = 4, vjust = 0.5))
  }
  
  return(plt_manhat)
}