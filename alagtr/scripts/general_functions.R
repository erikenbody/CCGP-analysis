#!/usr/bin/env Rscript

#' Import, process, and check gen objects and coordinates
#' 
#' TODO: eventually coords will be pulled from the CCGP-module directory
#' TODO add land_type to this function to specify env layers to gather for analyses
#' 
#' @param species species, including project ID (e.g., "41-Chamaea)
#' @param data_path path to raw data files
#' @param analysis type of gen object required as input into given analysis: "vcf" (default; for wingen, RDA, LFMM, TESS) or "gendist" (for GDM, MMRR)
#' @param impute if `analysis = "vcf"`, whether to impute missing values based on simple median-based imputation ("simple") or no imputation ("none"; default)
#' @param rmislands whether to remove islands (TRUE; default) or not (FALSE) depending on project's sampling
#' @param intervals boolean whether the analysis is running per contig or for the genome
#' @param scaff if intervals, which scaffold is being run
#' @param save_impute if `impute`, export imputed genotypes (defaults to FALSE)
#' @param incl_env whether to retrieve envlayers (defaults to TRUE)
#' @param env_var_type if `incl_env = TRUE`, type of environmental variables to gather; options are "rasterpcs" (default) or "bio1ndvi"
#' @param vcf_path path to vcf or gendist file
#' @param coords path to sampling coordinates
#' @param land_type if `env_var_type = "rasterpcs"`, landscape type for env layers; options are "terrestrial" (default) or "marine"
#' @param shape_path if `rmislands = TRUE`, path to shapefile that excludes Channel Islands
#' 
#'
#' @return list with genetic data and coordinates (in same order and matched samples)
#' @export
get_input_objects <- function(species, data_path, analysis = "gendist", impute = "none", 
                              rmislands = TRUE, intervals = FALSE, scaff = NA, save_impute = FALSE, 
                              incl_env = TRUE, env_var_type = "rasterpcs", vcf_path, coords, 
                              land_type = "terrestrial", shape_path) {
  # Get coords --------------------------------------------------------------
  coords <- readr::read_tsv(coords, col_names = c("INDV", "x", "y"))

  # Get input data ----------------------------------------------------------
  if (analysis == "vcf") {
    gen <- vcfR::read.vcfR(vcf_path)
    if (impute == "simple") {
      gen <- algatr::vcf_to_dosage(gen)
      gen <- algatr::simple_impute(x = gen, FUN = median)
    }
  }
  if (analysis == "gendist") {
    dist_id = paste0(vcf_path, ".id")
    # Process dists using algatr
    gen <- algatr::gen_dist(dist_type = "plink", plink_file = vcf_path, plink_id_file = dist_id)
  }
  
  # Check coords and gendist IDs --------------------------------------------
  dat <- check_ccgp_data(gen = gen, coords = coords, filetype = analysis)
  
  # Get envlayers -----------------------------------------------------------
  if (incl_env) {
    envlayers <- get_envlayers(env_var_type = env_var_type, 
                         shape_path = shape_path, 
                         rmislands = rmislands,
                         land_type = land_type)
    # envlayers <- get_envlayers(env_path = snakemake@params[[21]], 
    #                      shape_path = snakemake@params[[23]], 
    #                      layers = snakemake@params[[22]], 
    #                      rmislands = rmislands,
    #                      land_type = land_type,
    #                      shape_path = shape_path)
    # envlayers <- get_envlayers(env_path = "/scratch2/erik/CCGP-reruns/data/", rmislands = rmislands)
  } else {
    envlayers <- NULL
  }
  return(list(gen = gen, coords = dat$coords, sampleIDs = dat$sampleIDs, envlayers = envlayers))
}

#' Checks sample IDs from coordinates against those of genetic data file; also
#' checks ordering of samples within each of those. If genetic data are missing
#' for coordinates, removes those samples and outputs new coords file that's corrected.
#' If missing coords for samples, function will stop.
#' 
#'
#' @param gen genetic data
#' @param coords sampling coordinates
#' @param filetype file type to compare coords to: "gendist" (default) or "vcf"
#'
#' @export
check_ccgp_data <- function(gen, coords, filetype = "gendist"){
  # Get IDs from gen data
  if(filetype == "gendist") genID <- colnames(gen)
  
  if(filetype == "vcf") {
    if (inherits(gen, "vcfR")) {
      # colnames(gen@gt) <- stringr::str_replace_all(colnames(gen@gt), "[^[:alnum:]]", "_")
      genID <- colnames(gen@gt[, -1])
    }
    # If structure-based imputation was performed, output from str_impute is a dosage matrix: 
    if (!inherits(gen, "vcfR")) genID <- rownames(gen)
  }
  
  # Check dimensions match
  if(length(coords$INDV) != length(genID)) {warning(paste0("Number of individuals in coords (", length(coords$INDV) ,") and genetic data (", length(genID), ") do not match"))}
  
  # Check overlap
  overlap <- coords$INDV %in% genID
  if(!all(overlap)){warning("Missing genetic data for: ", paste(coords$INDV[!overlap]), ", removing coordinate data for these individuals...")}
  coordsF <- coords[overlap,]

  overlap <- genID %in% coords$INDV
  if(!all(overlap)){stop("Missing coordinate data for: ", paste(genID[!overlap]))}
  
  # Check order
  coordsF <- coordsF[match(genID, coordsF$INDV),]
  if(!all(coordsF$INDV == genID)){warning("Order of samples in coordinates and genetic data do not match")}
  
  # note: first column must be included because it is the format column
  # vcfF <- vcf[,c(TRUE, overlap)]
  # reorder and confirm order is the same
  # coordsF <- coordsF[match(colnames(vcfF@gt[,-1]), coordsF$INDV),]
  # if(!all(coordsF$INDV == colnames(vcfF@gt[,-1]))){warning("order of coords and vcf do not match")}
  final_coords <- coordsF %>% dplyr::select(-INDV)
  
  return(list(coords = final_coords, sampleIDs = genID))
}

#' Get PC envlayers and remove islands; requires env_path/PC_layers and env_path/CA_State
#' 
#' 
#' @param env_var_type type of env vars desired; options are "rasterpcs" or "bio1ndvi"
#' @param rmislands whether to remove islands or not
#' @param shape_path if `rmislands = TRUE`, path to shapefile
#' @param land_type if `env_var_type = "rasterpcs"`, landscape type for env layers; options are "terrestrial" (default) or "marine"
#'
#' @return envlayers
#' @export
get_envlayers <- function(env_var_type, rmislands, shape_path, land_type) {
  if (env_var_type == "rasterpcs") {
    if (land_type == "terrestrial") envlayers <- terra::rast(paste0(snakemake@scriptdir, "/../../data/california_chelsa_bioclim_1981-2010_V.2.1_pca.tif"))
    # if (land_type == "marine") envlayers <- terra::rast(paste0(snakemake@scriptdir, "/../../data/XXX_pca.tif"))
  }

  if (env_var_type == "bio1ndvi") {
    all_bio <- terra::rast(paste0(snakemake@scriptdir, "/../../data/california_chelsa_bioclim_1981-2010_V.2.1.tif"))
    # Extract only BIO1
    bio1 <- terra::subset(all_bio, "CHELSA_bio1_1981-2010_V.2.1")
    # Import NDVI layer
    ndvi <- terra::rast(paste0(snakemake@scriptdir, "/../../data/california_ndvi_mean_2000_2020.tif"))
    ndvi <- terra::resample(ndvi, bio1)
    envlayers <- c(bio1, ndvi)
    # Can do checks like crs(), ext(), or res() to ensure layers match
  }

  print(summary(envlayers))
  
  # Remove islands; requires rmapshaper package
  # TODO get rid of rgdal
  if(rmislands == TRUE) {
    # Shape file of region of interest
    # states <- tigris::states(cb = TRUE)
    # ca <- states[states$STUSPS == "CA", "STUSPS"]
    # ca <- sf::st_transform(ca, sf::st_crs(4326))
    # ca_noislands <- sf::st_read(here("data", "CA_State_2024_noChannelIslands", "CA_State_2024.shp"))
    ca_noislands <- sf::st_read(paste0(snakemake@scriptdir, shape_path))
    ca_noislands <- sf::st_as_sf(ca_noislands, crs = sf::st_crs(envlayers))
    ca_noislands <- sf::st_transform(ca_noislands, sf::st_crs(envlayers))
    envlayers <- terra::mask(envlayers, ca_noislands)
  }
  
  return(envlayers)
}

#' Convert from data frame to formatted sp
#'
#' @param coords sf object or matrix representing coordinates
#'
#' @return converted coords in sp format
#' @export
coords_to_df <- function(coords) {
  if (inherits(coords, "sf")) coords <- sf::st_coordinates(coords)
  if (is.matrix(coords)) coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")
  return(coords)
}

#' Reproject sampling coordinates and raster to same CRS
#'
#' @param coords sampling coordinates, in latitude/longitude degrees
#' @param env RasterStack of envlayers to be reprojected
#' @param newcrs CRS to reproject data to (defaults to 3310)
#'
#' @return reprojected coords
#' @export
reproject <- function(coords, env, newcrs = 3310) {
  latlong <- sf::st_as_sf(coords, coords = c("x", "y"), crs = "+proj=longlat")
  coords_proj <- sf::st_transform(latlong, crs = newcrs)
  
  # TODO convert to terra package
  # if (!is.null(env)) env_proj <- terra::projectRaster(env, crs = newcrs)
  if (!is.null(env)) env_proj <- raster::projectRaster(env, crs = newcrs)
  if (is.null(env)) env_proj <- NULL
  
  return(list(coords_proj = coords_proj, env_proj = env_proj))
}
