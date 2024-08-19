suppressMessages({
  library(ggplot2)
  library(here)
  library(dplyr)
})

#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

# example call: `Rscript admixture_cloud.R "59-Ursus" "outputs/admixture/"`

#!/usr/bin/env Rscript # leave line commented

species = snakemake@params[[1]]
output_path = snakemake@params[[2]]


# Get CV scores -----------------------------------------------------------

admix_cv <- function(output_path) {
  files <- list.files(output_path, pattern = ".log", full.names = TRUE)
  shortfiles <- list.files(output_path, pattern = ".log", full.names = FALSE)
  
  cv_scores <-
    1:length(files) %>% 
    lapply(function(x) {
      cv <- readLines(con = files[[x]])
      cv <- cv[grepl("^CV error", cv)]
      cv <- as.data.frame(cv) %>% dplyr::mutate(filename = shortfiles[[x]])
      return(cv)
    }) %>% 
    dplyr::bind_rows() %>% 
    tidyr::separate(cv, sep = ": ", into = c("temp", "cv")) %>% 
    dplyr::mutate(kval = readr::parse_number(temp)) %>% 
    tidyr::separate(filename, sep = "\\_", into = c("temp2", "temp3", "temp4", "tempseed")) %>% 
    tidyr::separate(tempseed, sep = "\\.", into = c("seed", "temp5")) %>% 
    dplyr::mutate(seed = readr::parse_number(seed)) %>% 
    dplyr::select(-c(temp, temp2, temp3, temp4, temp5)) %>% 
    dplyr::mutate(cv = as.numeric(cv), 
                  kval = as.numeric(kval), 
                  seed = as.numeric(seed)) %>% 
    dplyr::group_by(kval) %>% 
    dplyr::mutate(min_cv_seed = min(cv)) %>% 
    dplyr::ungroup()

  cv_scores_seeds <- cv_scores %>% 
    dplyr::arrange(kval) %>% 
    dplyr::select(kval, min_cv_seed) %>% 
    distinct()
  
  # Get min_k value across all Ks
  min_k <- cv_scores_seeds %>% 
    dplyr::filter(min_cv_seed == min(min_cv_seed)) %>%
    dplyr::pull(kval)
  
  # Get best_k value
  kvals <- cv_scores_seeds$kval
  ce.K <- cv_scores_seeds$min_cv_seed
  diff <- ce.K[-1] - ce.K[-max(kvals)]
  slope <- exp(-diff) - 1
  best_k <- min(which(slope <= quantile(slope)[4]))

  # Evanno's delta K value; adapted from code provided by Katherine Drotos
  delta_k <- cv_scores %>%
    dplyr::arrange(seed, kval) %>%
    dplyr::mutate(secondrate = abs((lead(cv) - cv) - (cv - lag(cv)))) %>%
    dplyr::group_by(kval) %>%
    dplyr::mutate(delta_k = mean(secondrate) / sd(cv)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(kval, delta_k)
  
  evanno_k <- delta_k$kval[which.max(delta_k$delta_k)]
  
  cv_scores <- cv_scores %>% 
    dplyr::mutate(min_k = min_k,
                  best_k = best_k,
                  evanno_k = evanno_k)
  
  return(cv_scores)
}

# Run above function
cv_scores <- admix_cv(output_path)


# Export results ----------------------------------------------------------

export_admix <- function(output_path, species, cv_scores) {
  readr::write_csv(cv_scores,
                   file = paste0(output_path, species, "_best_K.txt"),
                   col_names = TRUE)
  
  # Export plot of CV error
  cv_scores %>%
    dplyr::group_by(kval) %>% 
    dplyr::mutate(mean = mean(cv)) %>% 
    dplyr::mutate(sd = sd(cv)) %>% 
    dplyr::ungroup() %>% 
    ggplot2::ggplot(aes(x = kval, y = mean)) +
    ggplot2::geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2, color = "grey") +
    ggplot2::geom_point() +
    ggplot2::geom_line(aes(group = 1)) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept = unique(cv_scores$min_k),
                        color = "red",
                        linetype = "dashed") +
    ggplot2::annotate("text", x = unique(cv_scores$min_k-0.2), y = mean(cv_scores$cv), label = "Minimum K", size = 4, color = "red", angle = 90) +
    ggplot2::geom_vline(xintercept = unique(cv_scores$best_k),
                        color = "blue",
                        linetype = "dashed") +
    ggplot2::annotate("text", x = unique(cv_scores$best_k-0.2), y = mean(cv_scores$cv), label = "Best K", size = 4, color = "blue", angle = 90) +
    ggplot2::geom_vline(xintercept = unique(cv_scores$evanno_k), 
                        color = "purple", 
                        linetype = "dashed") +
    ggplot2::annotate("text", x = unique(cv_scores$evanno_k-0.2), y = max(cv_scores$cv), label = "ΔK", size = 4, color = "purple", angle = 90) +
    ggplot2::ggtitle(paste0(species, " admixture results")) +
    ggplot2::ylab("Cross-validation score") +
    ggplot2::scale_x_continuous(breaks = unique(cv_scores$kval))

  ggplot2::ggsave(paste0(output_path, species, "_admixture_cv.png"),
                  width = 20, height = 15, units = "cm")
  
  # Export Evanno's delta K plot; adapted from code provided by Katherine Drotos
  cv_scores %>%
    dplyr::arrange(seed, kval) %>%
    dplyr::mutate(secondrate = abs((lead(cv) - cv) - (cv - lag(cv)))) %>%
    dplyr::group_by(kval) %>%
    dplyr::mutate(delta_k = mean(secondrate) / sd(cv)) %>%
    dplyr::ungroup() %>% 
    dplyr::select(kval, delta_k) %>% 
    ggplot(aes(x = kval, y = delta_k)) + 
    ggplot2::geom_vline(aes(xintercept = unique(cv_scores$evanno_k)), color = "purple", linetype = "dashed", linewidth = 1) +
    geom_point() + 
    geom_line(aes(group = 1)) +
    scale_x_continuous(breaks = unique(cv_scores$kval)) +
    ylab("Δk") +
    ggtitle(paste0(species, " delta K"))
  
  ggplot2::ggsave(paste0(output_path, species, "_admixture_deltak.png"),
                  width = 20, height = 15, units = "cm")
}

export_admix(output_path, species, cv_scores)