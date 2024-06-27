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
})

output = snakemake@output[[1]]

if (!require("algatr", character.only = TRUE)) {
  # Install the package if not installed
  devtools::install_github("TheWangLab/algatr", quiet = T)
}

if (!require("wingen", character.only = TRUE)) {
  # Install the package if not installed
  install.packages("wingen")
}

if (!require("tess3r", character.only = TRUE)) {
  # Install the package if not installed
  devtools::install_github("bcm-uga/TESS3_encho_sen", quiet = T)
}

file.create(output)
