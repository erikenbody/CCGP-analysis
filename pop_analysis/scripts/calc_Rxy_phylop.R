#### PREREQUISITES #####
# Load necessary libraries
library(tidyverse)

#### Load MAF Files from Snakemake Input ####
maf_files <- snakemake@input
plot_out <- snakemake@output[["plot_out"]]
rxy_out <- snakemake@output[["rxy_out"]]
maf_files <- list.files("/scratch2/erik/CCGP-reruns/projects/68-Strix/results/GCA_030819815.1/pop_analysis/68-Strix_mutation_load", full.names = TRUE)
maf_files <- grep(".af.txt", maf_files, value = T)


load_maf_data <- function(files) {
  maf_data <- list()

  for (file in files) {
    # Extract variant type and population ID from filename
    #variant_type <- str_extract(basename(file), "^[a-z]+")  # Variant type
    #print(variant_type)
    #pop_id <- str_extract(basename(file), "\\d+")            # Population ID (1 or 2)
    pop_id <- strsplit(basename(file), "_")[[1]][3]  # Extract the 2nd element (1-based indexing)
    print(pop_id)
    # Load the MAF file and store under the corresponding population ID
    maf_data[[pop_id]] <- read_tsv(file, col_names = FALSE)
    names(maf_data[[pop_id]]) <- c("chrom", "pos","alt_freq")
    
  }

  return(maf_data)
}

# Load MAF data into a structured list
maf_data <- load_maf_data(maf_files)

#### Calculate Rxy Statistics for Population Pairs ####
# Function to calculate L(A, B) statistic
calculate_L <- function(pop1, pop2) {
  sum(pop1$alt_freq * (1 - pop2$alt_freq)) / sum(pop2$alt_freq * (1 - pop1$alt_freq))
}

# Function to calculate Rxy for population pairs (for each variant type)
calculate_Rxy <- function(pop1_df, pop2_df) {
  L_AB <- calculate_L(pop1_df, pop2_df)
  L_BA <- calculate_L(pop2_df, pop1_df)
  return(L_AB / L_BA)
}

Rxy <- data.frame(PopPair = character(), Rxy = numeric())
# Get all unque population names (assuming they are the list names in maf_data)
populations <- names(maf_data)
# Generate all unique pairs of populations (combinations without repetition)
pop_pairs <- combn(populations, 2, simplify = FALSE)

# Iterate over each population pair and calculate Rxy
for (pair in pop_pairs) {
  pop1 <- pair[1]
  pop2 <- pair[2]
  
  # Extract data for the two populations
  pop1_data <- maf_data[[pop1]]
  pop2_data <- maf_data[[pop2]]

  # Create a unique index to match variants by chromosome and position
  pop1_data$idx <- paste(pop1_data$chrom, pop1_data$pos, sep = "_")
  pop2_data$idx <- paste(pop2_data$chrom, pop2_data$pos, sep = "_")

  # Inner join to keep only shared variants between the two populations
  df_comb <- inner_join(pop1_data, pop2_data, by = c("chrom", "pos"))
  df_comb$idx <- paste(df_comb$chrom, df_comb$pos, sep = "_")
  # Filter the original data to only include shared variants
  pop1_data <- pop1_data %>% filter(idx %in% df_comb$idx)
  pop2_data <- pop2_data %>% filter(idx %in% df_comb$idx)
  print(pop1_data)
  print(pop2_data)
  # Calculate Rxy for the current population pair
  Rxy_value <- calculate_Rxy(pop1_data, pop2_data)

  # Store the result in the Rxy dataframe
  Rxy <- rbind(Rxy, data.frame(PopPair = paste0(pop1, "-", pop2), Rxy = Rxy_value))
}

t.test(df_comb$alt_freq.x, df_comb$alt_freq.y)

# #temporarily, lets load up the chamaea phylop

# phylop <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/SIFT_merged.tsv")
# phylop <- phylop %>% 
#  mutate(PhyloP = ifelse(PhyloP == "None", NA, as.numeric(PhyloP)),
#         SIFT_PREDICTION = ifelse(grepl("WARNING", SIFT_PREDICTION), "DELETERIOUS(Low)", SIFT_PREDICTION))

# # Compute means and standard deviations by SIFT_PREDICTION and VARIANT_TYPE
# sift_means <- phylop %>%
#   group_by(SIFT_PREDICTION) %>%
#   summarise(mean_phylop = mean(PhyloP, na.rm = TRUE))

# variant_means <- phylop %>%
#   group_by(VARIANT_TYPE) %>%
#   summarise(mean_phylop = mean(PhyloP, na.rm = TRUE))

# # Compute 95% confidence interval for PhyloP mean
# phylop_mean <- mean(phylop$PhyloP, na.rm = TRUE)
# phylop_sd <- sd(phylop$PhyloP, na.rm = TRUE)
# n <- sum(!is.na(phylop$PhyloP))
# ci_95 <- phylop_mean + c(-1, 1) * 1.96 * (phylop_sd / sqrt(n))



# p1 <- ggplot(phylop, aes(x = PhyloP)) +
#   geom_histogram(bins = 50, fill = "lightgrey", color = "black", alpha = 0.7) +
#   geom_vline(data = sift_means, aes(xintercept = mean_phylop, color = SIFT_PREDICTION),
#              linetype = "solid", size = 1) +
#   geom_vline(xintercept = ci_95, linetype = "solid", color = "red", size = 1) +
#   labs(title = "Distribution of PhyloP Scores", x = "PhyloP Score", y = "Count") +
#   theme_bw() +
#   theme(text = element_text(size = 15)) 

# p2 <- ggplot(phylop, aes(x = PhyloP)) +
#   geom_histogram(bins = 50, fill = "lightgrey", color = "black", alpha = 0.7) +
#   geom_vline(data = variant_means, aes(xintercept = mean_phylop, color = VARIANT_TYPE),
#              linetype = "solid", size = 1) +
#   geom_vline(xintercept = ci_95, linetype = "solid", color = "red", size = 1) +
#   labs(title = "Distribution of PhyloP Scores", x = "PhyloP Score", y = "Count") +
#   theme_bw() +
#   theme(text = element_text(size = 15)) 


# library(patchwork)
# p1 + p2
# ggsave("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/freq/phylop_hist.png",
# width = 16, height = 8)

# p1 <- phylop %>%
#   ggplot() + 
#     geom_boxplot(aes(y = PhyloP, x = SIFT_PREDICTION)) +
#     labs(title = "PhyloP vs. SIFT Score", x = "PhyloP Score", y = "SIFT Score") +
#     theme_bw() +
#     theme(text = element_text(size = 15))

# p2 <- phylop %>%
#   ggplot() + 
#     geom_boxplot(aes(y = PhyloP, x = VARIANT_TYPE)) +
#     labs(title = "PhyloP vs. SIFT Score", x = "PhyloP Score", y = "SIFT Score") +
#     theme_bw() +
#     theme(text = element_text(size = 15))

# p1 + p2
# for_test <- phylop %>% 
#   filter(SIFT_PREDICTION %in% c("DELETERIOUS", "TOLERATED")) 

# t.test(for_test$PhyloP ~ for_test$SIFT_PREDICTION)
