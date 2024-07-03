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
  
  cv_scores <-
    1:length(files) %>% 
    lapply(function(x) {
      cv <- readLines(con = files[[x]])
      cv <- cv[grepl("^CV error", cv)]
      return(as.data.frame(cv))
    }) %>% 
    dplyr::bind_rows() %>% 
    tidyr::separate(cv, sep = ": ", into = c("temp", "cv")) %>% 
    dplyr::mutate(kval = readr::parse_number(temp)) %>% 
    dplyr::select(-temp) %>% 
    dplyr::mutate(cv = as.numeric(cv))
  
  # Get min_k value
  min_k <- cv_scores %>% 
    filter(cv == min(cv)) %>% 
    pull(kval)
  
  # Get best_k value
  ce.K <- cv_scores$cv
  diff <- ce.K[-1] - ce.K[-max(cv_scores$kval)]
  slope <- exp(-diff) - 1
  # K is selected based on the smallest slope value in the upper quartile
  best_k <- min(which(slope <= quantile(slope)[4]))
  
  cv_scores <- cv_scores %>% 
    dplyr::mutate(min_k = min_k,
                  best_k = best_k)
  
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
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = kval, y = cv)) +
    ggplot2::geom_line(ggplot2::aes(x = kval, y = cv)) +
    ggplot2::theme_bw() +
    ggplot2::geom_vline(xintercept = cv_scores$min_k, 
                        color = "red", 
                        linetype = "dashed") +
    ggplot2::annotate("text", x = cv_scores$min_k-0.2, y = mean(cv_scores$cv), label = "Minimum K", size = 4, color = "red", angle = 90) +
    ggplot2::geom_vline(xintercept = cv_scores$best_k, 
                        color = "blue", 
                        linetype = "dashed") +
    ggplot2::annotate("text", x = cv_scores$best_k-0.2, y = mean(cv_scores$cv), label = "Best K", size = 4, color = "blue", angle = 90) +
    ggplot2::ggtitle(paste0(species, " admixture results")) +
    ggplot2::ylab("Cross-validation score") +
    ggplot2::scale_x_continuous(breaks = cv_scores$kval)
  
  ggplot2::ggsave(paste0(output_path, species, "_admixture_cv.png"), 
                  width = 20, height = 15, units = "cm")
}

export_admix(output_path, species, cv_scores)