#### PREREQUISITES #####
# Load necessary libraries
library(tidyverse)

#### Load MAF Files from Snakemake Input ####
maf_files <- snakemake@input
plot_out <- snakemake@output[["plot_out"]]
rxy_out <- snakemake@output[["rxy_out"]]
maf_files <- list.files("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/freq", full.names = TRUE)
maf_files <- grep("frq", maf_files, value = T)

load_maf_data <- function(files) {
  maf_data <- list()

  for (file in files) {
    # Extract variant type and population ID from filename
    variant_type <- str_extract(basename(file), "^[a-z]+")  # Variant type
    print(variant_type)
    pop_id <- str_extract(basename(file), "\\d+")            # Population ID (1 or 2)
    print(pop_id)
    # Initialize nested list if it doesn't exist yet
    if (!variant_type %in% names(maf_data)) {
      maf_data[[variant_type]] <- list()
    }

    # Load the MAF file and store under the corresponding population ID
    df <- read_tsv(file, skip = 1, col_names = FALSE)
    names(df) <- c("chrom", "pos", "dip","n", "ref_freq", "alt_freq")    
    df <- df %>% separate(alt_freq, into = c(NA, "alt_freq"), sep = ":") %>%
      mutate(alt_freq = as.numeric(alt_freq)) %>%
      select(chrom, pos, alt_freq)
    
    maf_data[[variant_type]][[pop_id]] <- df
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

# Initialize dataframe to store Rxy results
Rxy <- data.frame(Variant = character(), PopPair = character(), Rxy = numeric())

# Iterate over variant types and calculate Rxy for population pairs
for (variant in names(maf_data)) {
  pop1_data <- maf_data[[variant]][["1"]]
  pop2_data <- maf_data[[variant]][["2"]]

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

  # Calculate Rxy for this variant type
  Rxy_value <- calculate_Rxy(pop1_data, pop2_data)

  # Store result in the Rxy dataframe
  Rxy <- rbind(Rxy, data.frame(Variant = variant, PopPair = "1-2", Rxy = Rxy_value))
}

#### Plot the Rxy Results ####
Rxy$Variant <- factor(Rxy$Variant, levels = unique(Rxy$Variant))

ggplot(Rxy, aes(x = Variant, y = Rxy, fill = Variant)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Variant Type", y = "Rxy") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  theme_bw() +
  theme(text = element_text(size = 30), 
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")
ggsave(plot_out, width = 6, height = 8)

write_tsv(Rxy, rxy_out)

#temporarily, lets load up the chamaea phylop

phylop <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/SIFT_merged.tsv")
phylop <- phylop %>% 
 mutate(PhyloP = ifelse(PhyloP == "None", NA, as.numeric(PhyloP)),
        SIFT_PREDICTION = ifelse(grepl("WARNING", SIFT_PREDICTION), "DELETERIOUS(Low)", SIFT_PREDICTION))

# Compute means and standard deviations by SIFT_PREDICTION and VARIANT_TYPE
sift_means <- phylop %>%
  group_by(SIFT_PREDICTION) %>%
  summarise(mean_phylop = mean(PhyloP, na.rm = TRUE))

variant_means <- phylop %>%
  group_by(VARIANT_TYPE) %>%
  summarise(mean_phylop = mean(PhyloP, na.rm = TRUE))

# Compute 95% confidence interval for PhyloP mean
phylop_mean <- mean(phylop$PhyloP, na.rm = TRUE)
phylop_sd <- sd(phylop$PhyloP, na.rm = TRUE)
n <- sum(!is.na(phylop$PhyloP))
ci_95 <- phylop_mean + c(-1, 1) * 1.96 * (phylop_sd / sqrt(n))



p1 <- ggplot(phylop, aes(x = PhyloP)) +
  geom_histogram(bins = 50, fill = "lightgrey", color = "black", alpha = 0.7) +
  geom_vline(data = sift_means, aes(xintercept = mean_phylop, color = SIFT_PREDICTION),
             linetype = "solid", size = 1) +
  geom_vline(xintercept = ci_95, linetype = "solid", color = "red", size = 1) +
  labs(title = "Distribution of PhyloP Scores", x = "PhyloP Score", y = "Count") +
  theme_bw() +
  theme(text = element_text(size = 15)) 

p2 <- ggplot(phylop, aes(x = PhyloP)) +
  geom_histogram(bins = 50, fill = "lightgrey", color = "black", alpha = 0.7) +
  geom_vline(data = variant_means, aes(xintercept = mean_phylop, color = VARIANT_TYPE),
             linetype = "solid", size = 1) +
  geom_vline(xintercept = ci_95, linetype = "solid", color = "red", size = 1) +
  labs(title = "Distribution of PhyloP Scores", x = "PhyloP Score", y = "Count") +
  theme_bw() +
  theme(text = element_text(size = 15)) 


library(patchwork)
p1 + p2
ggsave("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/freq/phylop_hist.png",
width = 16, height = 8)

p1 <- phylop %>%
  filter(!grepl("Low", SIFT_PREDICTION) & !is.na(SIFT_PREDICTION)) %>%
  ggplot() + 
      geom_boxplot(aes(y = PhyloP, x = SIFT_PREDICTION)) +
    labs(title = "PhyloP vs. SIFT Score", y = "PhyloP Score", x = "SIFT Score") +
    theme_bw() +
    theme(text = element_text(size = 15)) +


p2 <- phylop %>%
  filter(!grepl("Low", SIFT_PREDICTION) & !is.na(SIFT_PREDICTION)) %>%
  ggplot() + 
    geom_boxplot(aes(y = PhyloP, x = VARIANT_TYPE)) +
    labs(title = "PhyloP vs. SIFT Score", y = "PhyloP Score", x = "SIFT Score") +
    theme_bw() +
    theme(text = element_text(size = 15))

p1 + p2
for_test <- phylop %>% 
  filter(SIFT_PREDICTION %in% c("DELETERIOUS", "TOLERATED")) 

t.test(for_test$PhyloP ~ for_test$SIFT_PREDICTION)
p1
ggsave("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/freq/phylop_sift_boxplot.png",
width = 6, height = 8)


##
pop1 <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/41-Chamaea_population_1_phylo_p_hom_count.txt")
pop2 <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/41-Chamaea_population_2_phylo_p_hom_count.txt")
mean(pop1$Homozygous_Alternate_Count)
mean(pop2$Homozygous_Alternate_Count)

pop1 <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/population_1_hom_count.txt")
pop2 <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/population_2_hom_count.txt")
mean(pop1$Homozygous_Alternate_Count)
mean(pop2$Homozygous_Alternate_Count)
t.test(pop1$Homozygous_Alternate_Count, pop2$Homozygous_Alternate_Count)
coords <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/algatr/41-Chamaea.coords.txt", col_names = F)
names(coords) <- c("Sample", "lon", "lat")

froh1 <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_merged_roh/41-Chamaea_filtered_population_1.froh", col_names = F)
names(froh1) <- c("Sample", "FROH")
froh1$pop <- "1"
froh2 <- read_tsv("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_merged_roh/41-Chamaea_filtered_population_2.froh", col_names = F)
names(froh2) <- c("Sample", "FROH")
froh2$pop <- "2"
froh <- bind_rows(froh1, froh2)

pop1$pop <- "1"
pop2$pop <- "2"

pops <- bind_rows(pop1, pop2) 
pops <- left_join(pops, coords, by = "Sample")
pops <- left_join(pops, froh, by = "Sample")

pops %>%
  ggplot(aes(x = pop, y = Homozygous_Alternate_Count)) + 
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_bw() +
  labs(x = "Population", y = "No. of homozygous deleterious sites") +
  theme(text = element_text(size = 20))
ggsave("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/freq/pop_boxplot.png")

# Load required libraries
library(sf)      # For handling spatial data
library(maps)    # For state outlines (California)

# Step 1: Load California map data
ca_map <- map_data("state") %>% filter(region == "california")
mean_count <- mean(pops$Homozygous_Alternate_Count, na.rm = TRUE)

# Step 2: Plot the map with points
ggplot() +
  geom_polygon(data = ca_map, aes(x = long, y = lat, group = group), 
               fill = "grey80", color = "black", size = 0.5) +
  geom_point(data = pops, 
             aes(x = lon, y = lat, color = Homozygous_Alternate_Count), 
             size = 10) +
scale_color_gradient2(
    low = "green", mid = "yellow", high = "red", 
    midpoint = mean_count,  # Center the color scale at the mean
    name = "Alt Count"
  ) +
  theme_minimal() +
  labs(color = "Alt Count") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed(1.3)  # Ensure correct aspect ratio
ggsave("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/freq/pop_map.png",)

ggplot() +
  geom_polygon(data = ca_map, aes(x = long, y = lat, group = group), 
               fill = "white", color = "black", size = 0.5) +
  geom_point(data = pops, 
             aes(x = lon, y = lat, color = pop), 
             size = 3, alpha = 0.8) 

ggplot(pops) +
  geom_point(aes(x = FROH, y = Homozygous_Alternate_Count)) +
  theme_bw() +
  labs(x = "FROH", y = "No. homozygous deleterious sites") +
  theme(text = element_text(size = 20))
ggsave("/scratch2/erik/CCGP-reruns/projects/41-Chamaea/results/GCA_029207755.1/pop_analysis/41-Chamaea_mutation_load/freq/froh_vs_hom.png")

