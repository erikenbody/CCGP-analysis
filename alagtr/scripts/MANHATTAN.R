
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

    return(plt_manhat)
}

rdadapt <- snakemake@input[["rdadapt_output"]]
rdadapt <- read_csv(rdadapt)
sig <- snakemake@params[["sig"]]
mantplot <- manhat_plot(rdadapt, sig)
ggsave(snakemake@output[["manhattan_plot"]], mantplot, width = 12, height = 5, dpi = 300, bg = "white")