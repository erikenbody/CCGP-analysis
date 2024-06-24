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

#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

output = snakemake@output[[1]]

if (!require("algatr", character.only = TRUE)) {
  # Install the package if not installed
  print("installing algatr")
  devtools::install_github("TheWangLab/algatr", quiet = F)
}

if (!require("wingen", character.only = TRUE)) {
  print("installing wingen")
  # Install the package if not installed
  devtools::install_github("AnushaPB/wingen", quiet = F)
}

if (!require("tess3r", character.only = TRUE)) {
  print("installing tess3r")
  # Install the package if not installed
  devtools::install_github("bcm-uga/TESS3_encho_sen", quiet = T)
}

if (!require("algatr", character.only = TRUE)) {
  print("algatr did not get installed nad I dont know why")
}

file.create(output)
