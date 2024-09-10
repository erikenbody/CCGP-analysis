library(ggplot2)
library(dplyr)

plot_gone <- function(gone_files, output_png, ccgp_id, population_id) {
  data_list <- lapply(names(gone_files), function(seed) {
    df <- read.table(gone_files[[seed]], sep = "\t", skip = 2, header = FALSE)
    colnames(df) <- c("Generation", "Ne")
    df$Ne <- as.numeric(df$Ne)
    df$Seed <- seed
    return(df)
  })
  
  combined_df <- bind_rows(data_list)
  
  p <- ggplot(combined_df, aes(x = Generation, y = Ne, color = Seed)) +
    geom_line() +
    labs(title = paste("Effective Population Size (Ne) over Generations -", ccgp_id, "-", population_id),
         x = "Generations",
         y = "Ne") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  ggsave(output_png, plot = p, width = 10, height = 6)
}

args <- commandArgs(trailingOnly = TRUE)
gone_files <- setNames(args[1:10], args[11:20])
output_png <- args[21]
ccgp_id <- args[22]
population_id <- args[23]

plot_gone(gone_files, output_png, ccgp_id, population_id)