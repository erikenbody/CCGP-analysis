suppressMessages({
  library(ggplot2)
  library(terra)
  library(raster)
  library(fields)
  library(rworldmap)
  library(automap)
  library(cowplot)
  library(here)
  library(LEA)
  library(dplyr)
  library(rgdal)
  library(peakRAM)
  library(devtools)
})

if (!require("algatr", character.only = TRUE)) {
  # Install the package if not installed
  devtools::install_github("TheWangLab/algatr", quiet = T)
}

if (!require("wingen", character.only = TRUE)) {
  # Install the package if not installed
  devtools::install_github("AnushaPB/wingen", quiet = T)
}

if (!require("tess3r", character.only = TRUE)) {
  # Install the package if not installed
  devtools::install_github("bcm-uga/TESS3_encho_sen", quiet = T)
}

library(algatr)
library(wingen)
library(tess3r)

#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

# example call: `Rscript TESS_rancor.R "59-Ursus" "~/../../media/WangLab/WangLab/CCGP_raw_data/" TRUE 1:5 "outputs/TESS/"`

#!/usr/bin/env Rscript # leave line commented

species = snakemake@params[[1]]
data_path = snakemake@params[[2]]
rmislands = as.logical(snakemake@params[[3]])
kvals = snakemake@params[[4]]
output_path = snakemake@params[[5]]
incl_env = as.logical(snakemake@params[[6]])
coords = snakemake@input[["coords"]]
vcf = snakemake@input[["vcf"]]
# species <- "5-Mirounga"
# data_path <- "projects/5-Mirounga/results/GCA_029215605.1/"
# rmislands = TRUE
# kvals <- "1:10"
# output_path <- "projects/5-Mirounga/results/GCA_029215605.1/algatr"
# incl_env <- TRUE
# source("/scratch2/erik/CCGP-reruns/alagtr/scripts/general_functions.R")

#need to interpret the kvals string as an expression
kvals <- try(eval(parse(text = kvals)), silent = TRUE)

source(paste0(snakemake@scriptdir, "/general_functions.R"))

# Import and process data -------------------------------------------------
peakRAM_imp <-
  peakRAM::peakRAM(
    dat <- get_input_objects(species = species, 
                             data_path = data_path,
                             analysis = "vcf",
                             pruned = TRUE,
                             impute = "none",
                             rmislands = rmislands,
                             incl_env = incl_env,
                             vcf_path = vcf,
                             coords = coords)
  )


# Convert to dosage matrix ------------------------------------------------

peakRAM_dos <-
  peakRAM::peakRAM(
    gen <- wingen::vcf_to_dosage(dat$gen)
  )


# Run K selection test and TESS -------------------------------------------

run_tess <- function(dat) {
  tess3_result <- algatr::tess_ktest(gen, 
                                     dat$coords,
                                     Kvals = kvals, 
                                     ploidy = 2, 
                                     K_selection = "auto")
  # Extract TESS object and best K
  tess3_obj <- tess3_result$tess3_obj
  bestK <- tess3_result[["K"]]

  # Get Q matrix
  qmat <- tess3r::qmatrix(tess3_obj, K = bestK)
  
  # Get minK value
  ce <- list()
  for (k in kvals) ce[[k]] <- tess3_obj[[k]]$crossentropy
  ce.K <- c()
  for (k in kvals) ce.K[k] <- min(ce[[k]])
  
  ce_results <- as.data.frame(ce.K) %>% 
    dplyr::mutate(Kval = kvals) %>% 
    dplyr::filter(ce.K == min(ce.K))
  
  minK <- ce_results$Kval

  return(list(tess3_obj = tess3_obj, bestK = bestK, minK = minK, qmat = qmat, krig_raster = tess3_result$grid))
}

peakRAM_tess <-
  peakRAM::peakRAM(
    results <- run_tess(dat)
  )


# Krige Q values ----------------------------------------------------------

# Uncomment if you want kriging
# fails for K=1, so removing

# if (!is.null(incl_env) & length(kvals) > 1) {
#   grid <- raster::aggregate(dat$envlayers[[1]], fact = 6)
#   reproj <- reproject(coords = dat$coords, env = grid)
  
#   peakRAM_krig <-
#     peakRAM::peakRAM(
#       krig_admix <- algatr::tess_krig(results$qmat, 
#                                       reproj$coords_proj, 
#                                       reproj$env_proj)
#     )
# }


# Export results ----------------------------------------------------------

export_TESS <- function(dat, results) {
  if (length(kvals) > 1) {
    ks <- 1:length(kvals)
    qvals <- 
      ks %>%
      lapply(function(x) {
        pops <- algatr::pops_helper(gen = gen, tess3_obj = results$tess3_obj, K = x)
        pops <- pops %>% 
          tidyr::pivot_longer(cols = "K1":paste0("K", x), names_to = "kval", values_to = "qvalue") %>% 
          dplyr::mutate(total_K = paste0(x))
        return(pops)
      }) %>% 
      dplyr::bind_rows()
    
    qvals$best_k = unique(results$bestK)
    qvals$min_k = unique(results$minK)
    
    readr::write_csv(qvals,
                     file = paste0(output_path, species, "_TESS_qmatrix.csv"),
                     col_names = TRUE)

    # Export raster of kriged Q values
    # comment out if you make and want to save rasters

    # if (!is.null(incl_env)) {
    #   terra::writeRaster(krig_admix,
    #                      paste0(output_path, species, "_TESS_bestK_krigadmix.tif"),
    #                      overwrite = TRUE)
    # }
  }
  
  if (length(kvals) == 1) {
    pops <- algatr::pops_helper(gen = gen, tess3_obj = results$tess3_obj, K = kvals)
    readr::write_csv(pops, 
                     file = paste0(output_path, species, "_TESS_qmatrix.csv"),
                     col_names = TRUE)
  }

  # Export cross-entropy values for all K values
  x <- results$tess3_obj
  
  xent <- seq_along(x)
  rmse <- seq_along(x)
  K <- seq_along(x)
  
  if (length(x) > 1) {
    for (i in seq_along(x)) {
      K[i] <- x[[i]]$K
      xent[i] <- x[[i]]$crossentropy
      rmse[i] <- x[[i]]$rmse
    }
  }
  
  if (length(x) == 1) {
    K <- x$K
    xent <- x$crossentropy
    rmse <- x$rmse
  }
  
  xval <- dplyr::bind_cols(as.data.frame(K), as.data.frame(xent), as.data.frame(rmse))
  readr::write_csv(xval,
                   file = paste0(output_path, species, "_TESS_xval.csv"),
                   col_names = TRUE)

  # Export plot of CV error
  if (length(x) > 1) {
    xval %>% 
      ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x = K, y = xent)) +
      ggplot2::geom_line(ggplot2::aes(x = K, y = xent)) +
      ggplot2::theme_bw() +
      ggplot2::geom_vline(xintercept = qvals$min_k, 
                          color = "red", 
                          linetype = "dashed") +
      ggplot2::annotate("text", x = qvals$min_k-0.2, y = mean(xval$xent), label = "Minimum K", size = 4, color = "red", angle = 90) +
      ggplot2::geom_vline(xintercept = qvals$best_k, 
                          color = "blue", 
                          linetype = "dashed") +
      ggplot2::annotate("text", x = qvals$best_k-0.2, y = mean(xval$xent), label = "Best K", size = 4, color = "blue", angle = 90) +
      ggplot2::ggtitle(paste0(species)) +
      ggplot2::ylab("Cross-entropy value") +
      ggplot2::scale_x_continuous(breaks = kvals)
    
    ggplot2::ggsave(paste0(output_path, species, "_TESS_xval.png"), 
           width = 20, height = 15, units = "cm")
  }

}

peakRAM_exp <-
  peakRAM::peakRAM(
    export_TESS(dat, results)
  )


# Export RAM usage --------------------------------------------------------

RAM <- dplyr::bind_rows(as.data.frame(peakRAM_imp),
                        as.data.frame(peakRAM_dos),
                        as.data.frame(peakRAM_tess),
                        as.data.frame(peakRAM_exp)) %>% 
  dplyr::mutate(fxn = c("import", "dosage", "run", "export"))

readr::write_csv(RAM,
                 file = paste0(output_path, species, "_TESS_peakRAM.csv"))