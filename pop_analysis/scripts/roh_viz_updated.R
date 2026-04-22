#!/usr/bin/env Rscript

# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(zoo))
suppressMessages(library(patchwork))
suppressMessages(library(gridExtra))
suppressMessages(library(ggridges))
suppressMessages(library(GenomicRanges))

# --- Get inputs from Snakemake ---
roh_path   <- snakemake@input$roh
top_path   <- snakemake@input$top
pi_paths   <- snakemake@input$pi_files
pop_id     <- snakemake@params$population

# --- Load and Process Data ---

all_pi_data <- purrr::map_dfr(pi_paths, function(file) {
    read_tsv(
        file,
        col_names = TRUE,
        show_col_types = FALSE,
        col_types = cols(
            CHROM = col_character(),
            BIN_START = col_double(),
            BIN_END = col_double(),
            N_VARIANTS = col_integer(),
            PI = col_double()
        )
    ) %>%
    mutate(
        pi_rm = zoo::rollmean(PI, 10, fill = NA, align = "center"),
        sample_id = stringr::str_remove(basename(file), "\\.windowed\\.pi$")
    )
})

df_runs <- read.table(roh_path, header = FALSE, col.names = c("chrom", "from", "to", "id")) %>% 
    mutate(lengthBps = to - from) %>%
    mutate(
        length_bin = case_when(
            lengthBps > 500000 ~ ">500kbp",
            lengthBps > 100000 ~ "<500kbp",
            TRUE               ~ "<100kbp"
        ),
        length_bin = factor(length_bin, levels = c("<100kbp", "<500kbp", ">500kbp"))
    )

top_samples <- read_lines(top_path)

# Check if top_samples is empty (no individuals with ROH > 500kbp)
if (length(top_samples) == 0 || all(is.na(top_samples)) || all(top_samples == "")) {
    warning(paste("Population", pop_id, "has no individuals with ROH > 500kbp. Creating placeholder plots."))

    # Create placeholder plots with explanatory text
    pdf(snakemake@output$roh1, width = 8.5, height = 11)
    plot.new()
    text(0.5, 0.5, paste("Population", pop_id, ":\nNo individuals with ROH > 500kbp\n\nThis population has no detectable\nruns of homozygosity greater than 500kbp"),
         cex = 1.5, col = "darkgray")
    dev.off()

    pdf(snakemake@output$roh2, width = 8.5, height = 8)
    plot.new()
    text(0.5, 0.5, paste("Population", pop_id, ":\nNo individuals with ROH > 500kbp"),
         cex = 1.5, col = "darkgray")
    dev.off()

    pdf(snakemake@output$roh3, width = 8, height = 4.5)
    plot.new()
    text(0.5, 0.5, paste("Population", pop_id, ":\nNo individuals with ROH > 500kbp"),
         cex = 1.5, col = "darkgray")
    dev.off()

    # Exit successfully since this is not an error condition
    quit(save = "no", status = 0)
}

sample1 <- top_samples[1]
sample2 <- if(length(top_samples) > 1) top_samples[2] else NA

pi1 <- filter(all_pi_data, sample_id == sample1)
roh1 <- filter(df_runs, id == sample1)
pi2 <- if(!is.na(sample2)) filter(all_pi_data, sample_id == sample2) else data.frame()
roh2 <- if(!is.na(sample2)) filter(df_runs, id == sample2) else data.frame()

if (nrow(pi1) == 0) {
    error_msg <- paste("Fatal Error: No parsable pi data found for the top F_ROH sample:", sample1,
                       "\nThis can happen if sample names do not match between the .froh and .pi files.",
                       "\nCannot generate plots for population", pop_id)
    stop(error_msg)
}

pi1$CHROM <- as.character(pi1$CHROM)
if(nrow(pi2) > 0) pi2$CHROM <- as.character(pi2$CHROM)
all_pi_data$CHROM <- as.character(all_pi_data$CHROM)
df_runs$chrom <- as.character(df_runs$chrom)

# --- Plot 1: All Individuals ROH + Pi ---
pdf(snakemake@output$roh1, width = 8.5, height = 11)

# List to hold plots for the current page
page_plots <- list()
plots_per_page <- 10

# Get list of large chromosomes from the top sample
large_chr <- pi1 %>%
  group_by(CHROM) %>%
  summarise(max_pos = max(BIN_END, na.rm = TRUE), .groups = 'drop') %>%
  filter(max_pos > 1000000) %>%
  pull(CHROM)

all_pop_samples <- unique(all_pi_data$sample_id)

# Loop through each chromosome and each individual to create plots
for (chr in large_chr) {
  for (indiv in all_pop_samples) {
    
    df_chr <- filter(all_pi_data, sample_id == indiv, CHROM == chr)
    box_filt <- filter(df_runs, id == indiv, chrom == chr)
    
    # Skip if there's no pi data for this specific combination
    if (nrow(df_chr) == 0) next

    pi_range <- (max(df_chr$PI, na.rm = TRUE) - min(df_chr$PI, na.rm = TRUE)) / 50
    
    p <- ggplot() +
      geom_line(data = df_chr, aes(x = BIN_START / 1e6, y = pi_rm), na.rm = TRUE) +
      geom_rect(data = box_filt, aes(xmin = from / 1e6, xmax = to / 1e6, ymin = -pi_range, ymax = 0, fill = length_bin)) +
      theme_bw() +
      labs(x = NULL, y = expression(pi), subtitle = paste("Sample:", indiv, "| Chromosome:", chr)) +
      scale_fill_manual(values = c("<100kbp" = "grey80", "<500kbp" = "blue", ">500kbp" = "red"), name = "ROH Length", drop = FALSE) +
      theme(legend.position = "none")
    
    # Add the generated plot to our page list
    page_plots[[length(page_plots) + 1]] <- p
    
    # If the page is full (10 plots), print it and reset the list
    if (length(page_plots) == plots_per_page) {
      grid.arrange(grobs = page_plots, ncol = 1)
      page_plots <- list() # Clear the list for the next page
    }
  }
}

# After the loops, print any remaining plots on the final page
if (length(page_plots) > 0) {
  grid.arrange(grobs = page_plots, ncol = 1)
}

dev.off()

# --- Plot 2: Top 2 Individuals ---

# --- 💡 FIX 💡 ---
# The previous dplyr pipe (`group_by`, `summarise`, `slice`) was failing due to a package conflict.
# This block is rewritten in base R to bypass the problematic functions.
split_by_chrom <- split(pi1, pi1$CHROM)
max_ends <- sapply(split_by_chrom, function(df) max(df$BIN_END, na.rm = TRUE))
sorted_max_ends <- sort(max_ends, decreasing = TRUE)
biggest_chrom_name <- names(sorted_max_ends)[1]
biggest_chr_df <- data.frame(CHROM = biggest_chrom_name)
# --- End of Fix ---


pi_range1 <- (max(pi1$PI, na.rm = TRUE) - min(pi1$PI, na.rm = TRUE)) / 50
pA <- pi1 %>% filter(CHROM == biggest_chr_df$CHROM) %>%
    ggplot() + geom_line(aes(x = BIN_START/1e6, y = pi_rm), na.rm=TRUE) +
    geom_rect(data = filter(roh1, chrom == biggest_chr_df$CHROM), aes(xmin = from/1e6, xmax = to/1e6, ymin = -pi_range1, ymax = 0, fill = length_bin)) +
    theme_bw() + labs(x = NULL, title = paste(sample1, "(Top F_ROH) | Largest Scaffold"), y = expression(pi)) +
    scale_fill_manual(values = c("<100kbp" = "grey80", "<500kbp" = "blue", ">500kbp" = "red"), drop=FALSE)

pB <- ggplot() + theme_void() 
if(nrow(pi2) > 0){
    pi_range2 <- (max(pi2$PI, na.rm = TRUE) - min(pi2$PI, na.rm = TRUE)) / 50
    pB <- pi2 %>% filter(CHROM == biggest_chr_df$CHROM) %>%
        ggplot() + geom_line(aes(x = BIN_START/1e6, y = pi_rm), na.rm=TRUE) +
        geom_rect(data = filter(roh2, chrom == biggest_chr_df$CHROM), aes(xmin = from/1e6, xmax = to/1e6, ymin = -pi_range2, ymax = 0, fill = length_bin)) +
        theme_bw() + labs(x = paste(biggest_chr_df$CHROM, "(Mbp)"), title = paste(sample2, "(2nd Top F_ROH) | Largest Scaffold"), y = expression(pi)) +
        scale_fill_manual(values = c("<100kbp" = "grey80", "<500kbp" = "blue", ">500kbp" = "red"), drop=FALSE)
}

(pA / pB) + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave(snakemake@output$roh2, width = 8.5, height = 8)


# --- Plot 3: Ridges Plot for Top Individual ---
gr_pi <- with(pi1, GRanges(seqnames = CHROM, ranges = IRanges(BIN_START, BIN_END)))
gr_roh <- with(roh1, GRanges(seqnames = chrom, ranges = IRanges(from, to)))

overlaps <- findOverlaps(gr_pi, gr_roh)
pi1$in_roh <- "Outside ROH"
pi1$in_roh[queryHits(overlaps)] <- "Inside ROH"

ggplot(pi1, aes(x = PI, y = fct_rev(in_roh), fill = fct_rev(in_roh))) +
  geom_density_ridges(alpha = 0.8) +
  theme_bw() +
  scale_fill_viridis_d(guide = "none") +
  labs(y = NULL, x = expression(pi), title = paste("Pi Density for Sample:", sample1))
ggsave(snakemake@output$roh3, width = 8, height = 4.5)