suppressMessages({
  library(tidyverse)
  library(patchwork)
  library(maps)
  library(scatterpie)
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

create_plots <- function(coords_file, tess_file, xval_file, cv_errors_file, 
                          eigenvec_file, eigenval_file, admixture_path, output_path) {
  # Load the coordinate data
  coords <- read_tsv(coords_file, col_names = c("Sample", "Longitude", "Latitude"))
  
  # Load the Tess data
  tess <- read_csv(tess_file)
  names(tess) <- c("Sample", "pop_assignment", "kval", "qvalue", "total_K", "best_k", "min_k")
  tess_min_k <- unique(tess$min_k)
  tess_best_k <- unique(tess$best_k)
  
  # Load the cross-validation error data
  cv_errors <- read_delim(cv_errors_file, skip = 2, col_names = FALSE, delim = ": CV error=")

  cv_errors$X1 <- as.numeric(gsub("K=", "", cv_errors$X1))
  names(cv_errors) <- c("K", "CV_Error")
  
  # Load the xval data
  xval <- read_csv(xval_file)
  
  tess <- tess %>% filter(total_K == tess_best_k)
  df_tess <- left_join(tess, coords, by = "Sample")
  
  wide_data <- df_tess %>%
    pivot_wider(names_from = kval, values_from = qvalue)
  wide_data <- wide_data %>% rowwise() %>%
    mutate(pop_name = names(select(., starts_with("K")))[which.max(c_across(starts_with("K")))])
  
  # Load PCA data
  eigenvec <- read_tsv(eigenvec_file, col_names = FALSE, skip = 1)
  eigenval <- read_tsv(eigenval_file, col_names = FALSE)
  colnames(eigenvec) <- c("Sample", paste0("PC", 1:(ncol(eigenvec) - 1)))
  
  eigenvec <- left_join(eigenvec, wide_data, by = "Sample")
  
  # Read the first line of the cv_errors_file
  first_line <- readLines(cv_errors_file, n = 1)
  admix_best_k <- as.numeric(str_extract(first_line, "\\d+"))
  admix_min_k <- cv_errors %>%
    filter(CV_Error == min(CV_Error)) %>%
    pull(K)
  
  # Load the admixture data for the best K value
  num_columns <- as.numeric(admix_best_k)  # Convert admix_best_k to a numeric value
  column_names <- paste0("K", 1:num_columns)  # Generate column names from K1 to K{num_columns}

  # Read the admixture data with dynamic column names
  admixture_data <- read_delim(paste0(admixture_path, ".", admix_best_k, ".Q"), 
                              col_names = column_names, 
                              delim = " ")
  admixture_data$Sample <- eigenvec$Sample
  admixture_data <- inner_join(coords, admixture_data, by = "Sample")
  
  # Plot 1: Outline of California with pie plots for each individual
  california_map <- map_data("state", region = "california")
  k_columns <- names(admixture_data)[startsWith(names(admixture_data), "K")]
  
  p1 <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_scatterpie(data = admixture_data, aes(x = Longitude, y = Latitude), cols = k_columns, color = NA) +
    coord_fixed() +
    labs(title = "Admixture Proportions") +
    theme_minimal() +
    theme(legend.position = c(0.8, 0.8))
  
  # Plot 2: Cross-validation error plot
  p2 <- ggplot(cv_errors, aes(x = K, y = CV_Error)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = admix_min_k, color = "red", linetype = "dashed") +
    geom_vline(xintercept = admix_best_k, color = "blue", linetype = "dashed") +
    annotate("text", x = admix_min_k - 0.2, y = mean(cv_errors$CV_Error), label = "Minimum K", size = 4, color = "red", angle = 90) +
    annotate("text", x = admix_best_k - 0.2, y = mean(cv_errors$CV_Error), label = "Best K", size = 4, color = "blue", angle = 90) +
    labs(title = "Admixture CV Error", x = "K", y = "CV Error") +
    theme_bw() +
    ggplot2::scale_x_continuous(breaks = cv_errors$K)
  
  # Plot 3: TESS Admixture Proportions
  kvals <- unique(df_tess$kval)
  p3 <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_scatterpie(data = wide_data, aes(x = Longitude, y = Latitude), cols = kvals, color = NA) +
    coord_fixed() +
    labs(title = "TESS Admixture Proportions") +
    theme_minimal() +
    theme(legend.position = c(0.8, 0.8))
  
  # Plot 4: Tess cross-entropy error
  p4 <- ggplot(xval, aes(x = K, y = xent)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    geom_vline(xintercept = tess_min_k, color = "red", linetype = "dashed") +
    annotate("text", x = tess_min_k - 0.2, y = mean(xval$xent), label = "Minimum K", size = 4, color = "red", angle = 90) +
    geom_vline(xintercept = tess_best_k, color = "blue", linetype = "dashed") +
    annotate("text", x = tess_best_k - 0.2, y = mean(xval$xent), label = "Best K", size = 4, color = "blue", angle = 90) +
    labs(title = "Tess cross entropy error", x = "K", y = "Cross-entropy value") +
    scale_x_continuous(breaks = cv_errors$K)
  
  # Plot 5: PCA plot
  p6 <- ggplot(eigenvec, aes(x = PC1, y = PC2, color = pop_name)) +
    geom_point() +
    labs(title = "PCA colored by best K Tess", x = paste0("PC1 (", round(100 * eigenval[1, 1] / sum(eigenval), 2), "%)"),
         y = paste0("PC2 (", round(100 * eigenval[2, 1] / sum(eigenval), 2), "%)")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Combine all plots into one using patchwork
  combined_plot <- (p1 | p2) / (p3 | p4) / p6
  
  # Return the combined plot
  ggplot2::ggsave(filename = output_path, combined_plot, 
        width = 8.5, height = 11, units = "in")
}

# Example usage
combined_plot <- create_plots(
  coords_file = snakemake@input[["coords"]],
  tess_file = snakemake@input[["tess"]],
  xval_file = snakemake@input[["xval"]],
  cv_errors_file = snakemake@input[["cv_errors"]],
  eigenvec_file = snakemake@input[["eigenvec"]],
  eigenval_file = snakemake@input[["eigenval"]],
  admixture_path = snakemake@params[["admixture_path"]],
  output_path = snakemake@output[["out"]]
)

#create_plots <- function(coords_file, tess_file, xval_file, cv_errors_file, 
#                          eigenvec_file, eigenval_file, admixture_path, output_path) {
