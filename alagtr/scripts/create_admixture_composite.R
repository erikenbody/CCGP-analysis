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
                         eigenvec_file, eigenval_file, admixture_path, output_path, output_match) {
  # Load the coordinate data
  coords <- read_tsv(coords_file, col_names = c("Sample", "Longitude", "Latitude"))
  
  # Load the Tess data
  tess.in <- read_csv(tess_file)
  names(tess.in) <- c("Sample", "pop_assignment", "kval", "qvalue", "total_K", "best_k", "min_k")
  tess_min_k <- unique(tess.in$min_k)
  tess_best_k <- unique(tess.in$best_k)
  
  # Load the cross-validation error data
  cv_scores <- read_csv(cv_errors_file)
  admix_min_k <- unique(cv_scores$min_k)
  admix_best_k <- unique(cv_scores$best_k)
  admx_evank <- unique(cv_scores$evanno_k)
  
  # List of files
  log_path <- gsub("Q_files", "logs", dirname(admixture_path))
  files <- list.files(path = log_path, pattern = "_K[0-9]+_", full.names = T)
  files <- files[grepl("seed42.log", files)]
  
  # Function to extract Fst matrix from each admixture log file
  extract_fst <- function(file) {
    lines <- readLines(file)
    
    # Extract K value from the file name
    K <- as.numeric(str_match(file, "_K([0-9]+)")[, 2])
    
    # Find the start of the Fst matrix
    fst_start <- grep("Fst divergences between estimated populations:", lines) + 2
    
    # Extract the matrix data
    fst_data <- lines[fst_start:(fst_start + K - 1)]
    
    # Split the lines by tabs and remove empty strings
    numeric_data <- fst_data %>%
      str_split("\t+") %>%                       # Split lines by tabs
      purrr::map(~ .x[.x != ""]) %>%                    # Remove empty strings
      purrr::map(~ as.numeric(.x[-1]))                  # Convert to numeric, excluding the first entry
    
    # Construct data frame
    fst_matrix <- numeric_data[-1] %>%
      purrr::map_dfr(~ setNames(as.list(.x), paste0("Pop", 0:(length(.x) - 1))))
    
    # Add K and Population columns to the data frame
    fst_matrix <- fst_matrix %>%
      mutate(Population = paste0("Pop",c(1:(K-1))))
    
    fst_df <- fst_matrix %>% 
      gather(key = "Comparison", value = "Fst", -Population) %>% 
      filter(!is.na(Fst)) %>% 
      mutate(K = K)
    
    return(fst_df)
  }
  
  # Extract data from all files and combine into a single dataframe
  combined_fst_df <- purrr::map_df(files[-1], extract_fst)
  
  # Load the xval data
  xval <- read_csv(xval_file)
  
  tess <- tess.in %>% filter(total_K == tess_best_k)
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
  
  # Load the admixture data for the best K value
  num_columns <- as.numeric(admix_best_k)  # Convert admix_best_k to a numeric value
  column_names <- paste0("K", 1:num_columns)  # Generate column names from K1 to K{num_columns}
  
  # Read the admixture data with dynamic column names
  
  q_files <- list.files(path = dirname(admixture_path), pattern = "\\.42\\.Q$", full.names = TRUE)
  
  # Function to read and process each .Q file
  process_q_file <- function(file) {
    admix_k <- str_extract(file, "\\d+(?=\\.42\\.Q)")
    column_names <- paste0("K", seq_len(as.numeric(admix_k)))
    admixture_data <- read_delim(file, col_names = column_names, delim = " ", show_col_types = FALSE)
    admixture_data$Sample <- eigenvec$Sample
    admixture_data <- inner_join(coords, admixture_data, by = "Sample")
    admixture_data$K_value <- as.numeric(admix_k)
    return(admixture_data)
  }
  
  all_admixture_data <- map_df(q_files, process_q_file)
  
  #for simplicity, just load admixture_data as the relavent Q. Easier to plot
  admixture_data <- read_delim(paste0(admixture_path, ".", admix_best_k, ".42.Q"), 
                               col_names = column_names, 
                               delim = " ")
  admixture_data$Sample <- eigenvec$Sample
  admixture_data <- inner_join(coords, admixture_data, by = "Sample")
  #
  # Plot 1: Outline of California with pie plots for each individual
  california_map <- map_data("state", region = "california")
  
  k_columns <- names(admixture_data)[startsWith(names(admixture_data), "K")]
  

  # merge tess and admixture ------------------------------------------------
  #get 
  tess_wide <- tess.in %>% 
    pivot_wider(id_cols = c(Sample, pop_assignment, total_K), names_from = kval, values_from = qvalue) %>% 
    dplyr::rename(K_value = total_K)
  
  df_comp <- inner_join(all_admixture_data, tess_wide, by = c("Sample", "K_value"))
  
  #calculate correlation coefficients
  results_list <- list()
  renamed_list <- list() 
  # Iterate over i from 1 to 10
  for (i in 1:10) {
    # Filter and select relevant columns
    df_smaller <- df_comp %>%
      select(-Longitude, -Latitude) %>%
      filter(K_value == i) %>%
      select(where(~ !any(is.na(.))))
    
    # Subset the relevant columns
    k_x <- df_smaller %>% select(starts_with("K") & ends_with(".x"))
    k_y <- df_smaller %>% select(starts_with("K") & ends_with(".y"))
    
    # Calculate the correlation matrix
    cor_matrix <- cor(k_x, k_y)
    if(i == 1){
      cor_matrix[1] <- 1
    }

    print(paste("Correlation matrix for K_value =", i))
    print(cor_matrix)
    
    # Find the highest correlation for each K.x column
    max_cor_indices <- apply(cor_matrix, 1, which.max)
    max_cor_names <- colnames(k_y)[max_cor_indices]
    
    # Combine the K.x and the corresponding best match K.y column names
    corresponding_columns <- tibble(
      K_x = colnames(k_x),
      K_y = max_cor_names
    )
    print(paste("Corresponding columns for K_value =", i))
    print(corresponding_columns)
    
    # Check for duplicates in K_y
    if (anyDuplicated(corresponding_columns$K_y)) {
      # If duplicates are found, set result to NA
      result <- tibble(K_value = i)
      for (j in seq_along(colnames(k_x))) {
        result[[paste0("K", j)]] <- NA
      }
    } else {
      # Proceed with renaming if no duplicates
      rename_vector <- setNames(gsub(".y", ".x", corresponding_columns$K_y), corresponding_columns$K_x)
      
      # Rename the columns
      df_renamed <- df_smaller %>%
        rename_with(~ rename_vector[.x], .cols = starts_with("K") & ends_with(".x"))
      
      # Check that this works
      k_x <- df_renamed %>% select(starts_with("K") & ends_with(".x"))
      k_y <- df_renamed %>% select(starts_with("K") & ends_with(".y"))
      
      cor_matrix <- cor(k_x, k_y)
      print(paste("Reordered correlation matrix for K_value =", i))
      print(cor_matrix)
      
      # Step 1: Remove suffixes from row and column names
      rownames(cor_matrix) <- gsub(".x$", "", rownames(cor_matrix))
      colnames(cor_matrix) <- gsub(".y$", "", colnames(cor_matrix))
      
      # Step 2: Order the rows and columns consistently
      k_cols <- unique(c(rownames(cor_matrix), colnames(cor_matrix)))
      k_cols <- sort(k_cols)
      
      # Step 3: Reorder the correlation matrix
      reordered_matrix <- cor_matrix[k_cols, k_cols, drop = FALSE]
      
      # Step 4: Extract diagonal correlations
      diag_cor <- diag(reordered_matrix)
      
      # Step 5: Create the data frame
      result <- tibble(
        K_value = i
      )
      
      # Add each K column dynamically
      for (j in seq_along(diag_cor)) {
        result[[paste0("K", j)]] <- diag_cor[j]
      }
    }
    
    # Add the result to the list
    results_list[[i]] <- result
    renamed_list[[i]] <- df_renamed
  }
  
  # Combine all results into a single data frame
  final_result <- bind_rows(results_list)
  renamed_dfs <- bind_rows(renamed_list)
  write_csv(renamed_dfs, output_match)

  # Reshape data from wide to long format
  heatmap_data <- final_result %>%
    pivot_longer(-K_value, names_to = "K_column", values_to = "Value") %>% 
    mutate(K_column = as.factor(as.numeric(gsub("K","", K_column))))

  #check pop assignment based on 
  renamed_tess_best_K <- renamed_list[[tess_best_k]] 
  
  temp <- renamed_tess_best_K %>% 
    filter(K_value == tess_best_k) %>% 
    select(ends_with(".x"), -K_value) 
  
  if (nrow(temp)== 0){
    p9 <- ggplot() + geom_blank() + labs(title = "no overlap in K assignments")
  } else {
    admx_assign <- temp %>% 
      mutate(Category = names(temp)[max.col(temp)]) %>% 
      mutate(pop_assignment = gsub("K", "", gsub(".x", "", Category, fixed = T))) %>% 
      mutate(Sample = renamed_tess_best_K$Sample) %>% 
      select(Sample, pop_assignment)
    
    compare_bestk_coefficients <- left_join(admx_assign, renamed_tess_best_K, by = "Sample")
    
    df_matches <- compare_bestk_coefficients %>% 
      filter(pop_assignment.x == pop_assignment.y)
    
    df_matches4plot <- df_matches %>%
      select(-K_value) %>% 
      pivot_longer(cols = starts_with("K"), names_to = "K_type", values_to = "K_value") %>%
      mutate(type = ifelse(grepl("x", K_type), "admixture", "tess")) %>% 
      group_by(Sample, type) %>%
      filter(K_value == max(K_value, na.rm = TRUE)) %>%
      ungroup() %>%
      select(Sample, K_value, type) %>%
      pivot_wider(id_cols = "Sample", values_from = "K_value", names_from = "type") %>% 
      mutate(qvalue_diff = abs(admixture - tess)) %>% 
      arrange(qvalue_diff) 
    
    df_mismatches <- compare_bestk_coefficients %>% 
      filter(pop_assignment.x!=pop_assignment.y)
    
    df_long <- df_mismatches %>%
      select(-K_value) %>% 
      pivot_longer(cols = starts_with("K"), names_to = "K_type", values_to = "K_value") %>%
      mutate(type = ifelse(grepl("x", K_type), "admixture", "tess")) %>% 
      group_by(Sample, type) %>%
      filter(K_value == max(K_value, na.rm = TRUE)) %>%
      ungroup() %>%
      select(Sample, K_value, type) %>%
      arrange(K_value)

    p9_title = paste0("Matching assignments for Tess best K value: K", tess_best_k)
  
    p9.5 <- df_matches4plot %>% 
      ggplot(aes(x = reorder(Sample, qvalue_diff), y = qvalue_diff)) +
        geom_col() +
        labs(x = NULL, y = "Abs. Q value difference", title = p9_title) +
        geom_hline(yintercept = mean(df_matches4plot$qvalue_diff), color = "red") +
        theme_bw() +
        theme(axis.text.x = element_blank()) +
        ylim(c(0,0.5))
    
    p9.25_title = paste0("Mismatching assignments for Tess best K value: K", tess_best_k)
    
    p9.25 <-   
      ggplot(df_long, aes(x = reorder(Sample, K_value), y = K_value, fill = type)) +
      geom_col(position = position_dodge()) +
      labs(x = NULL, y = "Max Q value", title = p9.25_title) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      lims(y = c(0,1)) +
      geom_hline(yintercept = 0.50, color = "red")
    
    p9 <- p9.5/p9.25
  
  }
  # plots -------------------------------------------------------------------

  
  p1 <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_scatterpie(data = admixture_data, aes(x = Longitude, y = Latitude), cols = k_columns, color = NA) +
    coord_fixed() +
    labs(title = "Admixture Proportions", x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.position = c(0.8, 0.8))
  
  # Plot 2: Cross-validation error plot
  p2 <- cv_scores %>%
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
    ggplot2::ggtitle("Admixture CV error") +
    ggplot2::labs(x = "K") +
    ggplot2::ylab("Cross-validation score") +
    ggplot2::scale_x_continuous(breaks = unique(cv_scores$kval))
  
  # Plot 3: TESS Admixture Proportions
  kvals <- unique(df_tess$kval)
  p3 <- ggplot() +
    geom_polygon(data = california_map, aes(x = long, y = lat, group = group), fill = "white", color = "black") +
    geom_scatterpie(data = wide_data, aes(x = Longitude, y = Latitude), cols = kvals, color = NA) +
    coord_fixed() +
    labs(title = "TESS Admixture Proportions", x = NULL, y = NULL) +
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
    scale_x_continuous(breaks = xval$K)
  
  # Plot 5: PCA plot
  p6 <- ggplot(eigenvec, aes(x = PC1, y = PC2, color = pop_name)) +
    geom_point() +
    labs(title = "PCA colored by best K Tess", x = paste0("PC1 (", round(100 * eigenval[1, 1] / sum(eigenval), 2), "%)"),
         y = paste0("PC2 (", round(100 * eigenval[2, 1] / sum(eigenval), 2), "%)")) +
    theme_bw() +
    theme(legend.position = "none")
  
  p7 <- combined_fst_df %>% 
    ggplot(aes(x = K, y = Fst)) +
    geom_boxplot(aes(group = K)) + 
    geom_point() +
    theme_bw() +
    labs(title = "Fst calculated by Admixture") +
    scale_x_continuous(breaks = unique(cv_scores$kval))
  
  
  # Create the heatmap with coefficients
  p8 <- ggplot(heatmap_data, aes(x = K_column, y = as.factor(K_value), fill = Value)) +
    geom_tile() +
    geom_text(aes(label = ifelse(is.na(Value), "", sprintf("%.2f", Value))), color = "black", size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey") +
    theme_minimal() +
    labs(x = "Population", y = "Number of Ks in run", fill = "Correlation") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  layout <- "
            AABB
            AABB
            CCDD
            EEFF
            GGHH
            "
  
  #combined_plot <- (p1 | p3) / (p2 | p4) / (p6 | p7) +
  combined_plot <- p1 + p3 + p2 + p4 + p8 + p9 + p7 + p6 +
    plot_annotation(tag_levels = 'A') +
    plot_layout(design = layout)
  
  ggplot2::ggsave(filename = output_path, combined_plot, 
                  width = 10.5, height = 15, units = "in", device = cairo_pdf)
}

combined_plot <- create_plots(
  coords_file = snakemake@input[["coords"]],
  tess_file = snakemake@input[["tess"]],
  xval_file = snakemake@input[["xval"]],
  cv_errors_file = snakemake@input[["cv_errors"]],
  eigenvec_file = snakemake@input[["eigenvec"]],
  eigenval_file = snakemake@input[["eigenval"]],
  admixture_path = snakemake@params[["admixture_path"]],
  output_path = snakemake@output[["out"]],
  output_match = snakemake@output[["match"]]
)

# combined_plot <- create_plots(
#   coords_file = "data/super_secret_structure_project/41-Cyanocitta.coords.txt",
#   tess_file = "data/super_secret_structure_project/41-Cyanocitta_TESS_qmatrix.csv",
#   xval_file = "data/super_secret_structure_project/41-Cyanocitta_TESS_xval.csv",
#   cv_errors_file = "data/super_secret_structure_project/59-Ursus_best_K.txt",
#   eigenvec_file = "data/super_secret_structure_project/41-Cyanocitta.eigenvec",
#   eigenval_file = "data/super_secret_structure_project/41-Cyanocitta.eigenval",
#   admixture_path = "data/super_secret_structure_project/Q_files/41-Cyanocitta",
#   output_path = "data/super_secret_structure_project/test_output.pdf"
# )

# coords_file = "/scratch2/erik/CCGP-reruns/projects/5-Mirounga/results/GCA_029215605.1/algatr/5-Mirounga.coords.txt"
# tess_file = "/scratch2/erik/CCGP-reruns/projects/5-Mirounga/results/GCA_029215605.1/algatr/5-Mirounga_TESS_qmatrix.csv"
# xval_file = "/scratch2/erik/CCGP-reruns/projects/5-Mirounga/results/GCA_029215605.1/algatr/5-Mirounga_TESS_xval.csv"
# cv_errors_file = "/scratch2/erik/CCGP-reruns/projects/5-Mirounga/results/GCA_029215605.1/algatr/admixture/logs/5-Mirounga_best_K.txt"
# eigenvec_file = "/scratch2/erik/CCGP-reruns/projects/5-Mirounga/results/GCA_029215605.1/algatr/5-Mirounga.eigenvec"
# eigenval_file = "/scratch2/erik/CCGP-reruns/projects/5-Mirounga/results/GCA_029215605.1/algatr/5-Mirounga.eigenval"
# admixture_path = "/scratch2/erik/CCGP-reruns/projects/5-Mirounga/results/GCA_029215605.1/algatr/admixture/Q_files/5-Mirounga"
