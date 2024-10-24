suppressMessages({
  library(algatr)
  library(vcfR)
  library(LEA)
  library(tidyverse)
  library(wingen)
  library(data.table)
  library(lfmm)
})

# example call: `Rscript LFMM_rancor.R "59-Ursus" "~/../../media/WangLab/WangLab/CCGP_raw_data/" "outputs/LFMM/" FALSE FALSE "simple" "ridge" 3 0.05 "fdr"`
# if (!require("algatr", character.only = TRUE)) {
#   # Install the package if not installed
#   devtools::install_github("TheWangLab/algatr", quiet = T)
# }

# set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

#!/usr/bin/env Rscript # leave line commented
species = snakemake@params[[1]]
data_path = snakemake@params[[2]]
output_path = snakemake@params[[3]]
rmislands = as.logical(snakemake@params[[4]])
pruned = as.logical(snakemake@params[[5]])
impute = snakemake@params[[6]] # "simple"
lfmm_method = snakemake@params[[7]] # "ridge"
manual_k_assignment = snakemake@params[[8]] # from admixture/TESS
sig = snakemake@params[[9]] # 0.05
p_adj = snakemake@params[[10]] # "fdr"
intervals = as.logical(snakemake@params[[11]])
scaff = as.character(snakemake@params[[12]])
color_by_contig = as.character(snakemake@params[[13]]) # whether to colorize manhat plot non-outliers by contig

source(paste0(snakemake@scriptdir, "/general_functions.R"))


# Import and process data -------------------------------------------------

peakRAM_imp <-
  peakRAM::peakRAM(
    dat <- get_input_objects(species = species, 
                             data_path = data_path,
                             analysis = "vcf",
                             pruned = pruned,
                             impute = impute,
                             save_impute = FALSE, # save this for later
                             intervals = intervals,
                             scaff = scaff,
                             vcf_path = snakemake@input[["vcf"]],
                             kvals = NULL,
                             rmislands = rmislands)
  )
print("Imported LFMM data")

# Extract and standardize environmental variables and make into dataframe
env <- raster::extract(dat$envlayers, dat$coords)
env <- scale(env, center = TRUE, scale = TRUE)
env <- data.frame(env)
print("Extracted environmental values")


# Define functions --------------------------------------------------------

#' Run LFMM; modified version of algatr function `lfmm_run()`
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR` object
#' @param env dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted
#' @param K_value number of latent factors (if left as NULL (default), K value selection will be conducted)
#' @param lfmm_method lfmm method (either \code{"ridge"} (default) or \code{"lasso"})
#' @param p_adj method to use for p-value correction (defaults to "fdr"); other options can be found in \code{\link{p.adjust}}
#' @param sig alpha level for determining candidate SNPs (defaults to 0.05)
#' @param calibrate a character string, "gif" or "median+MAD". If the "gif" option is set (default), 
#' significance values are calibrated by using the genomic control method. Genomic control uses a robust 
#' estimate of the variance of z-scores called "genomic inflation factor". If the "median+MAD" option is set, 
#' the pvalues are calibrated by computing the median and MAD of the zscores. If NULL, the pvalues are not 
#' calibrated.
#' 
#' @return
#' @export
lfmm_run_ccgp <- function(gen, env, manual_k_assignment, lfmm_method = "ridge", p_adj = "fdr", sig = 0.05, calibrate = "gif") {
  # Check that env var names don't match coord names
  if (any(colnames(coords) %in% colnames(env))) {
    colnames(env) <- paste(colnames(env), "_env", sep = "")
    warning("env names should differ from x and y. Appending 'env' to env names")
  }
  
  # Handle NA values -----------------------------------------------------
  if (any(is.na(gen))) {
    stop("Missing values found in gen data")
  }
  
  if (any(is.na(env))) {
    warning("Missing values found in env data, removing rows with NAs")
    na_env <- env
    gen <- gen[complete.cases(na_env), ]
    # Must come last
    env <- env[complete.cases(na_env), ]
    if (!is.null(coords)) coords <- coords[complete.cases(na_env), ]
  }
  
  # gen matrix
  genmat <- as.matrix(gen)
  # env matrix
  envmat <- as.matrix(env)
  
  # Run model
  if (lfmm_method == "ridge") {
    lfmm_mod <- lfmm::lfmm_ridge(genmat, envmat, K = manual_k_assignment)
  }
  if (lfmm_method == "lasso") {
    lfmm_mod <- lfmm::lfmm_lasso(genmat, envmat, K = manual_k_assignment)
  }
  
  # Perform association testing using the fitted model:
  lfmm_test_result <- lfmm::lfmm_test(
    Y = genmat,
    X = envmat,
    lfmm = lfmm_mod,
    calibrate = calibrate
  )
  
  # If p_adj method is specified, perform p-value correction by column (by env variable)
  lfmm_test_result$adjusted.pvalue <- apply(dplyr::as_tibble(lfmm_test_result$calibrated.pvalue), 2, stats::p.adjust, method = p_adj)
  
  # Stop if all p-values are NA
  if (all(is.na(lfmm_test_result$adjusted.pvalue))) stop("All p-values are NA")
  
  # Transfer column names
  colnames(lfmm_test_result$adjusted.pvalue) <- colnames(envmat)
  
  # Transfer rownames
  rownames(lfmm_test_result$adjusted.pvalue) <- colnames(genmat)
  
  # Make tidy dataframe of results
  result_df <- algatr::lfmm_df(lfmm_test_result)
  
  # Subset out candidate SNPs
  lfmm_snps <- result_df %>% dplyr::filter(adjusted.pvalue < 0.05)
  
  return(list(lfmm_snps = lfmm_snps, 
              df = result_df, 
              model = lfmm_mod, 
              lfmm_test_result = lfmm_test_result, 
              K = manual_k_assignment))
}

#' Build Manhattan plots of results; modified from algatr's `lfmm_manhattanplot()` function
#'
#' @param results 
#' @param sig 
#' @param color_by_contig 
#'
#' @return
#' @export
manhat_plot <- function(results, sig = 0.05, color_by_contig = FALSE) {
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


# Run LFMM ----------------------------------------------------------------

peakRAM_run <-
  peakRAM::peakRAM(
    results <- lfmm_run_ccgp(gen = dat$gen, 
                             env = env, 
                             K_value = K_value, 
                             lfmm_method = lfmm_method, 
                             p_adj = p_adj, 
                             sig = sig, 
                             calibrate = "gif")
  )
print("LFMM has run!")


# Export results ----------------------------------------------------------

export_lfmm <- function(results) {
  # LFMM model results
  saveRDS(results$model, file = paste0(output_path, species, "_LFMM_model.RDS"))
  
  # Association testing results
  saveRDS(results$lfmm_test_result, file = paste0(output_path, species, "_LFMM_assoctest.RDS"))
  
  # All SNPs and their associations
  readr::write_csv(results$df,
                   file = paste0(output_path, species, "_LFMM_pvalues.csv"),
                   col_names = TRUE)
  
  # Build and save Manhattan plot to file
  plt_manhat <- manhat_plot(results, sig = sig, color_by_contig = color_by_contig) # color_by_contig = color_by_contig
  plt_manhat
  ggsave(paste0(output_path, species, "_LFMM_manhatplot.png"), width = 10, height = 8, bg = "white")
}

# Export results
peakRAM_exp <-
  peakRAM::peakRAM(
    export_lfmm(results)
  )


# Export RAM usage --------------------------------------------------------

RAM <- dplyr::bind_rows(peakRAM_imp,
                        peakRAM_run, 
                        peakRAM_exp) %>% 
  dplyr::mutate(fxn = c("import", "run", "export"))

readr::write_csv(RAM,
                 file = paste0(output_path, species, "_LFMM_peakRAM.csv"))
