library(algatr)
library(tidyverse)
library(peakRAM)
algatr::rda_packages()

# example call: `Rscript RDA_cloud_simple.R "59-Ursus" "~/../../media/WangLab/WangLab/CCGP_raw_data/" "outputs/RDA/" TRUE FALSE "simple" "full" NULL "~/../../media/WangLab/WangLab/CCGP_raw_data/59-Ursus/QC/59-Ursus.eigenvec" 3 0.05 1000 TRUE "fdr" 0.01 TRUE`

#!/usr/bin/env Rscript # leave line commented
args = commandArgs(trailingOnly=TRUE)
species = args[1]
data_path = args[2]
output_path = args[3]
rmislands = args[4]
pruned = args[5]
impute = args[6]
model = args[7]
correctGEO = args[8] # NULL or non-null string
correctPC = args[9] # must be path to Plink PCA results or NULL
nPC = args[10]
Pin = args[11]
R2permutations = args[12]
R2scope = args[13]
p_adj = args[14] # for outlier_method = "p"
sig = args[15] # for outlier_method = "p"
save_impute = args[16]
kvals = args[17] # for structure-based imputation

# Exported files are as follows:
#     (1) species_imputed_`impute`.txt # if `save_impute = TRUE`
#     (2) species_RDA_anova_`model`.csv # if `model = "best"`
#     (3) species_RDA_cortest_`model`.csv
#     (4) species_RDA_outliers_`model`_rdadapt.csv
#     (5) species_RDA_outliers_`model`_Zscores.csv
#     (6) species_RDA_peakRAM.csv

source("./general_functions.R")


# Import and process data -------------------------------------------------

peakRAM_imp <-
  peakRAM::peakRAM(
    dat <- get_input_objects(species = species, 
                             data_path = data_path,
                             analysis = "vcf",
                             gendata = gendata, # SWITCH TO PRUNED FOR CLOUD!
                             impute = impute,
                             kvals = kvals,
                             rmislands = rmislands,
                             save_impute = FALSE) # save this for later
  )

# Extract and standardize environmental variables and make into dataframe
env <- raster::extract(dat$envlayers, dat$coords)
env <- scale(env, center = TRUE, scale = TRUE)
env <- data.frame(env)


# Run RDA -----------------------------------------------------------------

#' Redefined algatr `rda_run()` function to take in Plink PC files for pRDA correctPC
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR` object
#' @param env dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted
#' @param coords dataframe with coordinates (only needed if correctGEO = TRUE) or if env is a Raster* from which values should be extracted
#' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); defaults to "full"
#' @param correctGEO whether to condition on geographic coordinates (defaults to FALSE)
#' @param correctPC path to eigenvector file produced from Plink --pca to condition on PCs from PCA of genotypes; defaults to NULL
#' @param nPC if `correctPC` not NULL, number of PCs to use (defaults to 3)
#' @param Pin if `model = "best"`, limits of permutation P-values for adding (`Pin`) a term to the model, or dropping (`Pout`) from the model. Term is added if` P <= Pin`, and removed if `P > Pout` (see \link[vegan]{ordiR2step}) (defaults to 0.05)
#' @param R2permutations if `model = "best"`, number of permutations used in the estimation of adjusted R2 for cca using RsquareAdj (see \link[vegan]{ordiR2step}) (defaults to 1000)
#' @param R2scope if `model = "best"` and set to TRUE (default), use adjusted R2 as the stopping criterion: only models with lower adjusted R2 than scope are accepted (see \link[vegan]{ordiR2step})
#'
#' @return RDA model
#' @export
rda_run_pc <- function(gen, env, coords = NULL, model = "full", correctGEO = FALSE, correctPC = NULL, nPC = 3,
                       Pin = 0.05, R2permutations = 1000, R2scope = T) {
  # Format coordinates ------------------------------------------------------
  if (!is.null(coords)) coords <- coords_to_df(coords)
  
  # Check that env var names don't match coord names
  if(any(colnames(coords) %in% colnames(env))) {
    colnames(env) <- paste(colnames(env), "_env", sep = "")
    warning("env names should differ from x and y. Appending 'env' to env names")
  }
  
  # Handle NA values -----------------------------------------------------
  if (any(is.na(gen))) {
    stop("Missing values found in gen data")
  }
  
  if (any(is.na(env))) {
    warning("Missing values found in env data, removing rows with NAs")
    gen <- gen[complete.cases(env), ]
    if (!is.null(coords)) coords <- coords[complete.cases(env), ]
    # NOTE: this must be last
    env <- env[complete.cases(env), ]
  }
  
  # Set up model ---------------------------------------------------------
  if (is.null(correctPC) & is.null(correctGEO)) {
    moddf <- data.frame(env)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+")))
  }
  
  if (!is.null(correctPC) & is.null(correctGEO)) {
    pc <- read_tsv(paste0(correctPC)) %>% 
      tibble::column_to_rownames(var = "#IID") %>% 
      dplyr::select(tidyselect::all_of(1:nPC))
    
    # Check env var naming ----------------------------------------------------
    if(any(colnames(pc) %in% colnames(env))) {
      colnames(env) <- paste(colnames(env), "_env", sep = "")
      warning("env names should differ from PC1, PC2, etc if correctPC is TRUE. Appending 'env' to env names")
    }
    moddf <- data.frame(env, pc)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(", paste(colnames(pc), collapse = "+"), ")"))
  }
  
  if (!is.null(correctPC) & !is.null(correctGEO)) {
    if (is.null(coords)) stop("Coordinates must be provided if correctGEO is TRUE")
    pc <- read_tsv(paste0(correctPC)) %>% 
      tibble::column_to_rownames(var = "#IID") %>% 
      dplyr::select(tidyselect::all_of(1:nPC))
    
    # Check env var naming ----------------------------------------------------
    if(any(colnames(pc) %in% colnames(env))) {
      colnames(env) <- paste(colnames(env), "_env", sep = "")
      warning("Enviro var names should differ from PC1, PC2, etc if correctPC is TRUE. Appending env to enviro var names")
    }
    moddf <- data.frame(env, coords, pc)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(", paste(colnames(pc), collapse = "+"), "+ x + y)"))
  }
  
  if (is.null(correctPC) & !is.null(correctGEO)) {
    if (is.null(coords)) stop("Coordinates must be provided if correctGEO is TRUE")
    moddf <- data.frame(env, coords)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(x + y)"))
  }
  
  if (model == "best") {
    mod_full <- vegan::rda(f, data = moddf)
    mod_null <- vegan::rda(gen ~ 1, data = moddf)
    mod <- vegan::ordiR2step(mod_null, mod_full, Pin = Pin, R2permutations = R2permutations, R2scope = R2scope)
    if (mod$call == mod_null$call) {
      mod <- NULL
      stop("Best model is NULL model, returning NULL")
    }
  } else {
    mod <- vegan::rda(f, data = moddf)
  }
  
  return(mod)
}

peakRAM_run <-
  peakRAM::peakRAM(
    mod <- rda_run_pc(gen = dat$gen,
                      env = env,
                      coords = dat$coords, 
                      model = model, 
                      correctGEO = correctGEO, 
                      correctPC = correctPC, 
                      nPC = nPC, 
                      Pin = Pin, 
                      R2permutations = R2permutations, 
                      R2scope = R2scope)
  )


# Get outliers and run cortest --------------------------------------------

get_outliers <- function(mod) {
  rda_sig_z <- rda_getoutliers(mod, naxes = "all", outlier_method = "z", z = 3, plot = FALSE)
  rda_sig_p <- rda_getoutliers(mod, naxes = "all", outlier_method = "p", p_adj = p_adj, sig = sig, plot = FALSE)
  
  # Extract genotypes for outlier SNPs
  rda_snps_p <- rda_sig_p$rda_snps
  rda_gen_p <- dat$gen[, rda_snps_p]
  rda_snps_z <- rda_sig_z$rda_snps
  rda_gen_z <- dat$gen[, rda_snps_z]
  
  results <- list(rda_sig_z = rda_sig_z,
                  rda_sig_p = rda_sig_p,
                  rda_gen_p = rda_gen_p,
                  rda_gen_z = rda_gen_z)
  
  return(results)
}

peakRAM_outliers <-
  peakRAM::peakRAM(
    outliers <- get_outliers(mod)
  )

# Run correlation test
run_cortest <- function(rda_gen_p, rda_gen_z) {
  # Run correlation test
  cor_df_p <- algatr::rda_cor(rda_gen_p, env)
  cor_df_z <- algatr::rda_cor(rda_gen_z, env)
  
  results <- list(cor_df_p = cor_df_p, 
                  cor_df_z = cor_df_z)
  
  return(results)
}

peakRAM_cortest <-
  peakRAM::peakRAM(
    cortest <- run_cortest(rda_gen_p = outliers$rda_gen_p, rda_gen_z = outliers$rda_gen_z)
  )


# Export results ----------------------------------------------------------

export_rda <- function(mod, rda_sig_z, rda_sig_p, cor_df_p, cor_df_z, save_impute) {
  outlier_helper <- function(df, outlier) {dat <- df %>% dplyr::mutate(outlier_method = outlier)}
  
  # RDA model results
  if (!is.null(mod$anova)) {
    readr::write_csv(mod$anova,
                     file = paste0(output_path, species, "_RDA_anova_", model, ".csv"),
                     col_names = TRUE)
  }
  
  # Sig results Z-scores
  readr::write_csv(rda_sig_z,
                   file = paste0(output_path, species, "_RDA_outliers_", model, "_Zscores.csv"),
                   col_names = TRUE)
  
  # Sig results p-values
  snps <- rda_sig_p$rdadapt %>% 
    dplyr::mutate(locus = colnames(dat$gen))
  readr::write_csv(snps,
                   file = paste0(output_path, species, "_RDA_outliers_", model, "_rdadapt.csv"),
                   col_names = TRUE)
  
  # Correlation test results
  cor_test <- rbind(outlier_helper(cor_df_p, outlier = "p"),
                    outlier_helper(cor_df_z, outlier = "z"))
  readr::write_csv(cor_test, file = paste0(output_path, species, "_RDA_cortest_", model, ".csv"),
                   col_names = TRUE)
  
  # Save imputed data
  if (save_impute) {
    write.table(dat$gen, file = paste0(output_path, species, "_imputed_", impute, ".txt"),
                sep = " ", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

# Export results
peakRAM_exp <- 
  peakRAM::peakRAM(
    export_rda(mod, outliers$rda_sig_z, outliers$rda_sig_p, cortest$cor_df_p, cortest$cor_df_z, save_impute = save_impute)
  )


# Export RAM usage --------------------------------------------------------

RAM <- dplyr::bind_rows(peakRAM_imp,
                        peakRAM_run, 
                        peakRAM_outliers,
                        peakRAM_cortest,
                        peakRAM_exp) %>% 
  dplyr::mutate(fxn = c("import", "run", "outliers", "cortest", "export"))

readr::write_csv(RAM,
                 file = paste0(output_path, species, "_RDA_peakRAM.csv"))
