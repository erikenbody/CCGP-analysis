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
                             incl_env = incl_env)
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
  
  return(list(tess3_obj = tess3_obj, bestK = bestK, qmat = qmat, krig_raster = tess3_result$grid))
}

peakRAM_tess <-
  peakRAM::peakRAM(
    results <- run_tess(dat)
  )


# Krige Q values ----------------------------------------------------------

if (!is.null(incl_env)) {
  grid <- raster::aggregate(dat$envlayers[[1]], fact = 6)
  reproj <- reproject(coords = dat$coords, env = grid)
  
  peakRAM_krig <-
    peakRAM::peakRAM(
      krig_admix <- algatr::tess_krig(results$qmat, 
                                      reproj$coords_proj, 
                                      reproj$env_proj)
    )
}


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
    readr::write_csv(qvals,
                     file = paste0(output_path, species, "_TESS_qmatrix.csv"),
                     col_names = TRUE)
  }
  
  if (length(kvals) == 1) {
    pops <- algatr::pops_helper(gen = gen, tess3_obj = results$tess3_obj, K = kvals)
    readr::write_csv(pops, 
                     file = paste0(output_path, species, "_TESS_qmatrix.csv"),
                     col_names = TRUE)
  }
  
  # Export raster of kriged Q values
  if (!is.null(incl_env)) {
    terra::writeRaster(krig_admix,
                       paste0(output_path, species, "_TESS_bestK_krigadmix.tif"),
                       overwrite = TRUE)
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
}

peakRAM_exp <-
  peakRAM::peakRAM(
    export_TESS(dat, results)
  )


# Export RAM usage --------------------------------------------------------

if (!is.null(incl_env)) {
  RAM <- dplyr::bind_rows(as.data.frame(peakRAM_imp),
                          as.data.frame(peakRAM_dos),
                          as.data.frame(peakRAM_tess),
                          as.data.frame(peakRAM_krig),
                          as.data.frame(peakRAM_exp)) %>% 
    dplyr::mutate(fxn = c("import", "dosage", "run", "krig", "export"))
}

if (is.null(incl_env)) {
  RAM <- dplyr::bind_rows(as.data.frame(peakRAM_imp),
                          as.data.frame(peakRAM_dos),
                          as.data.frame(peakRAM_tess),
                          as.data.frame(peakRAM_exp)) %>% 
    dplyr::mutate(fxn = c("import", "dosage", "run", "export"))
}

readr::write_csv(RAM,
                 file = paste0(output_path, species, "_TESS_peakRAM.csv"))