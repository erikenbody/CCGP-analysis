#!/usr/bin/env Rscript

#' Import, process, and check gen objects and coordinates
#' 
#' TODO: eventually coords won't be pulled from the QC dir but rather the spreadsheet
#' 
#' @param species species, including project ID (e.g., "41-Chamaea)
#' @param data_path path to raw data files
#' @param analysis type of gen object required as input into given analysis: "gendist" (default; for GDM, MMRR) or "vcf" (for wingen, RDA, LFMM, TESS)
#' @param pruned whether to retrieve LD-pruned data (TRUE; default) or not (FALSE)
#' @param impute if `analysis = "vcf"`, whether to impute missing values based on simple median-based imputation ("simple"), structure-based using sNMF ("structure") or no imputation ("none"; default)
#' @param kvals if `analysis = "vcf"` and `impute = "structure"`, integer vector corresponding to the number of ancestral populations for which the sNMF algorithm estimates have to be calculated
#' @param save_impute if `impute`, export imputed genotypes (defaults to FALSE)
#' @param incl_env whether to retrieve envlayers (defaults to TRUE)
#' @param rmislands whether to remove islands (TRUE; default) or not (FALSE) depending on project's sampling
#' @param intervals boolean whether the analysis is running per contig or for the genome
#' @param scaff if intervals, which scaffold is being run
#'
#' @return list with genetic data and coordinates (in same order and matched samples)
#' @export
get_input_objects <- function(species, data_path, analysis = "gendist", pruned = TRUE, impute = "none", 
  kvals = NULL, rmislands = TRUE, intervals = FALSE, scaff = NA, save_impute = FALSE, incl_env = TRUE, vcf_path) {
  
  # Get coords --------------------------------------------------------------

  coords <- readr::read_tsv(paste0(data_path, "algatr/", species, ".coords.txt", sep = ""), col_names = FALSE)
  colnames(coords) <- c("INDV", "x", "y")
  coords$x <- as.numeric(coords$x)
  coords$y <- as.numeric(coords$y)

  # Get input data ----------------------------------------------------------
  if (analysis == "gendist") {
    if (!pruned) {
      dist = paste0(data_path, "CCGP/", species, "_filtered.dist", sep = "")
      dist_id = paste0(data_path, "CCGP/", species, "_filtered.dist.id", sep = "")
    }
    if (pruned) {
      dist = paste0(data_path, "CCGP/", species, "_annotated_pruned_0.6.dist", sep = "")
      dist_id = paste0(data_path, "CCGP/", species, "_annotated_pruned_0.6.dist.id", sep = "")
    }
    # Process dists using algatr
    gen <- algatr::gen_dist(dist_type = "plink", plink_file = dist, plink_id_file = dist_id)
  }
  
  if (analysis == "vcf") {
      gen <- vcfR::read.vcfR(vcf_path)

    # if (!pruned) {
    #   if (!intervals){
    #     gen <- vcfR::read.vcfR(paste0(data_path, "CCGP/", species, "_annotated.vcf.gz"))
    #   }
    #   if (intervals){
    #     #gen <- vcfR::read.vcfR(paste0(data_path, "algatr/subsets/", species, "_", scaff, "_annotated.vcf.gz"))
    #     gen <- vcfR::read.vcfR(paste0(data_path, "algatr/subsets/", species, "_", scaff, "_annotated_pruned_0.6.vcf.gz"))
    #   }
    # }
    # if (pruned) {
    #   #gen <- vcfR::read.vcfR(paste0(data_path, "CCGP/", species, "_annotated.vcf.gz"))
    #   #gen <- vcfR::read.vcfR(paste0(data_path, "CCGP/", species, "_annotated_pruned_0.6.vcf.gz"))
    #   #use this to test out, much smaller file
    #   #gen <- vcfR::read.vcfR(paste0(data_path, "QC/", species, ".pruned.vcf.gz"))
    #   gen <- vcfR::read.vcfR(paste0(data_path, "algatr/", species, "_complete_coords_pruned_0.6.vcf.gz"))
    # }
    if (impute == "simple") {
      gen <- algatr::vcf_to_dosage(gen)
      gen <- algatr::simple_impute(x = gen, FUN = median)
    }
    if (impute == "structure") gen <- algatr::str_impute(gen = gen, K = kvals) # N.B.: output is a dosage matrix, not a vcfR object
  }
  
  # Check coords and gendist IDs --------------------------------------------

  # Filter gen object to include only individuals present in coords$INDV
  # genID <- colnames(gen@gt[, -1])
  # gen <- gen[, genID %in% coords$INDV]
  # print(gen)
  # print(gen[,"_E0044B"])

  dat <- check_ccgp_data(gen = gen, coords = coords, filetype = analysis)
  
  # Get envlayers -----------------------------------------------------------
  if (incl_env) {
    env <- get_envlayers(env_path = snakemake@params[[21]], shape_path = snakemake@params[[23]], layers = snakemake@params[[22]], rmislands = rmislands)
    #env <- get_envlayers(env_path = "/scratch2/erik/CCGP-reruns/data/", rmislands = rmislands)
  } else {
    env <- NULL
  }
  
  return(list(gen = gen, coords = dat$coords, sampleIDs = dat$sampleIDs, envlayers = env))
}

#' Checks sample IDs from coordinates against those of genetic data file; also
#' checks ordering of samples within each of those. If genetic data are missing
#' for coordinates, removes those samples and outputs new coords file that's corrected.
#' If missing coords for samples, function will stop.
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
#' @param env_path path to PC layers
#' @param shape_path path to shapefile
#' @param layers names of layers to be retained in raster stack; options are "all" (for all layers in tif; default) or list of names
#' @param rmislands whether to remove islands or not
#'
#' @return envlayers
#' @export
get_envlayers <- function(env_path, shape_path, layers = "all", rmislands = FALSE){
  # Get env layers
  # env_files <- list.files(paste0(env_path, "PC_layers"), full.names = TRUE)
  envlayers <- raster::stack(paste0(snakemake@scriptdir, env_path))

  # Subset particular layers of env object if so desired
  if (all(layers != "all")) envlayers <- raster::subset(x = envlayers, subset = layers)
  print(summary(envlayers))
  
  # Shape file of region of interest
  spdf <- rgdal::readOGR(paste0(snakemake@scriptdir, shape_path))
  boundary <- sp::spTransform(spdf, raster::crs("+proj=longlat +datum=WGS84 +no_defs"))
  
  # Remove islands
  if(rmislands == TRUE) envlayers <- algatr::rm_islands(envlayers, boundary)
  
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
#'
#' @return reprojected coords
#' @export
reproject <- function(coords, env) {
  latlong <- sf::st_as_sf(coords, coords = c("x", "y"), crs = "+proj=longlat")
  coords_proj <- sf::st_transform(latlong, crs = 3310)
  
  if (!is.null(env)) env_proj <- raster::projectRaster(env, crs = 3310)
  if (is.null(env)) env_proj <- NULL
  
  return(list(coords_proj = coords_proj, env_proj = env_proj))
}
