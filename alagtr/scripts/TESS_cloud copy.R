suppressMessages(library(ggplot2))
suppressMessages(library(terra))
suppressMessages(library(raster))
suppressMessages(library(fields))
suppressMessages(library(rworldmap))
suppressMessages(library(automap))
suppressMessages(library(cowplot))
suppressMessages(library(here))
suppressMessages(library(LEA))
suppressMessages(library(dplyr))
suppressMessages(library(rgdal))

# example call: `Rscript TESS_rancor.R "59-Ursus" "~/../../media/WangLab/WangLab/CCGP_raw_data/" TRUE 1:5 "outputs/TESS/"`

devtools::install_github("bcm-uga/TESS3_encho_sen")
devtools::install_github("AnushaPB/wingen")
devtools::install_github("TheWangLab/algatr")

suppressMessages(library(algatr))
suppressMessages(library(wingen))
suppressMessages(library(tess3r))

#!/usr/bin/env Rscript # leave line commented

species = snakemake@params[[1]]
data_path = snakemake@params[[2]]
rmislands = snakemake@params[[3]]
kvals = snakemake@params[[4]]
output_path = snakemake@params[[5]]

#need to interpret the kvals string as an expression
#kvals <- try(eval(parse(text = kvals)), silent = TRUE)
kvals <- 1:3
#similar problem with logical true/false
if (rmislands == "true") {
  rmislands <- TRUE
} else if (rmislands == "false") {
  rmislands <- FALSE
} else {
  stop("Invalid logical value provided as an argument. Please use 'True' or 'False'.")
}

source(paste0(snakemake@scriptdir, "/general_functions.R"))

# Import and process data -------------------------------------------------

dat <- get_input_objects(species = species, 
                         data_path = data_path,
                         analysis = "vcf",
                         pruned = TRUE,
                         impute = "none",
                         kvals = kvals,
                         rmislands = rmislands)

# Run TESS ----------------------------------------------------------------

results <- algatr::tess_do_everything(dat$gen, 
                                      dat$coords, 
                                      grid = raster::aggregate(dat$envlayers[[1]], fact = 6), 
                                      Kvals = kvals, 
                                      K_selection = "auto",
                                      plot_method = "maxQ", 
                                      col_breaks = 20,
                                      minQ = 0.10,
                                      tess_method = "projected.ls", 
                                      lambda = 1, 
                                      ploidy = 2, 
                                      correct_kriged_Q = TRUE,
                                      quiet = TRUE)

# Export results ----------------------------------------------------------

# Export Q values with population assignments
export_TESS <- function(dat, results) {
  
  ks <- 1:length(results$Kvals)

  if (length(results$Kvals) > 1) {
    dat <- 
      ks %>%
      lapply(function(x) {
        pops <- algatr::pops_helper(gen = dat$gen, tess3_obj = results$tess_results, K = x)
        pops <- pops %>% 
          tidyr::pivot_longer(cols = "K1":paste0("K", x), names_to = "kval", values_to = "qvalue") %>% 
          dplyr::mutate(total_K = paste0(x))
        return(pops)
        }) %>% 
      dplyr::bind_rows()
    
    readr::write_csv(dat,
                     file = paste0(output_path, species, "_TESS_qmatrix.csv"),
                     col_names = TRUE)
  }
  
  if (length(results$Kvals) == 1) {
    readr::write_csv(results$pops, 
                     file = paste0(output_path, species, "_TESS_qmatrix.csv"),
                     col_names = TRUE)
  }
  
  # Export raster of kriged Q values
  terra::writeRaster(results$krig_admix,
                     paste0(output_path, species, "_TESS_bestK_krigadmix.tif"),
                     overwrite = TRUE)
  
  # Export cross-entropy values for all K values
  x <- results$tess_results
  
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

export_TESS(dat, results)
