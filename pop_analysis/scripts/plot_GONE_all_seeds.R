# Title: Replicate GONE Plotting in R with Additional Outputs
# Description: This script reads output files from GONE, plots the effective 
#              population size (Ne) for multiple seeds, and plots the average Ne 
#              with a confidence interval. It also saves the summary statistics
#              and a table of Ne change at specific generations. It is designed 
#              to integrate with a Snakemake workflow.

# Load necessary libraries
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)


#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

#' Plot Ne, save summary stats, and calculate Ne change at specific generations.
#'
#' @param gone_files A named list where names are seeds and values are paths to GONE output files.
#' @param output_png Path for saving the plot of individual seed Ne trajectories.
#' @param output_avg_png Path for saving the plot of average Ne.
#' @param output_summary_tsv Path for saving the summary statistics table (mean, CI).
#' @param output_change_tsv Path for saving the Ne change table for specific generations.
#' @param ccgp_id The project identifier for plot titles.
#' @param population_id The population identifier for plot titles.
plot_ne <- function(gone_files, output_png, output_avg_png, output_summary_tsv, output_change_tsv, ccgp_id, population_id) {

    # Read and combine all GONE output files into a single dataframe.
    # The `purrr::map_dfr` function iterates through the list of files,
    # reads each one, and binds them together into a single data frame.
    # The `.id = "Seed"` argument automatically creates a column named "Seed"
    # from the names of the `gone_files` list.
    combined_df <- map_dfr(gone_files, ~read_tsv(., skip = 2, col_names = c("Generation", "Ne"), col_types = "id"), .id = "Seed")

    # Plot 1: Ne for all individual seeds
    #-----------------------------------------
    ggplot(combined_df, aes(x = Generation, y = Ne, group = Seed, color = Seed)) +
        geom_line(alpha = 0.8) +
        labs(
            title = paste("GONE -", ccgp_id, "- population", population_id),
            x = "Generations",
            y = "Ne"
        ) +
        theme_bw() +
        theme(panel.grid = element_line(color = "lightgray"))

    ggsave(output_png, width = 10, height = 6, dpi = 300)

    # Calculate mean and standard deviation for each generation
    summary_df <- combined_df %>%
        group_by(Generation) %>%
        summarise(
            mean = mean(Ne, na.rm = TRUE),
            std = sd(Ne, na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        mutate(
            lower_CI = mean - std,
            upper_CI = mean + std
        ) %>%
        mutate(across(c(mean, std, lower_CI, upper_CI), ~round(., 1)))

    write_tsv(summary_df, output_summary_tsv)


    # Plot 2: Average Ne with a confidence interval (mean +/- 1 std)
    #-----------------------------------------------------------------
    ggplot(summary_df, aes(x = Generation)) +
        geom_ribbon(aes(ymin = lower_CI, ymax = upper_CI), fill = "lightgray", alpha = 0.5) +
        geom_line(aes(y = mean, color = "Mean Ne")) +
        scale_color_manual(values = c("Mean Ne" = "blue")) +
        labs(
            title = paste("GONE -", ccgp_id, "- population", population_id),
            x = "Generations",
            y = "Average Ne",
            color = ""
        ) +
        theme_bw() +
        theme(
            panel.grid = element_line(color = "lightgray"),
            legend.position = "none"
        )

    # Save the average plot
    ggsave(output_avg_png, width = 10, height = 6, dpi = 300, bg = "white")

    # --- Extract Ne at specific generations and calculate change ---
    target_generations <- c(10, 50, 100, 150, 200)

    change_df <- summary_df %>%
        filter(Generation %in% target_generations) %>%
        select(Generation, Ne = mean) %>%
        mutate(Ne_change = Ne - lead(Ne, order_by = Generation)) %>%
        mutate(across(c(Ne, Ne_change), ~round(., 1)))

    # Save the resulting change dataframe
    write_tsv(change_df, output_change_tsv)
}
main <- function() {
  
    gone_files <- snakemake@input[["GONE_out_Ne"]]
    names(gone_files) <- snakemake@params[["seeds"]]
    
    # Get output paths from Snakemake object
    output_png <- snakemake@output[["GONE_plot"]]
    output_avg_png <- snakemake@output[["GONE_avg_plot"]]
    output_summary_tsv <- snakemake@output[["GONE_summary_stats"]] # New output
    output_change_tsv <- snakemake@output[["GONE_Ne_change"]]     # New output

    # Get parameters from Snakemake object
    ccgp_id <- snakemake@params[["project_id"]]
    population_id <- snakemake@params[["population_id"]]
    
    # Call the main plotting function with the new arguments
    plot_ne(gone_files, output_png, output_avg_png, output_summary_tsv, output_change_tsv, ccgp_id, population_id)
}

# This mimics the Python `if __name__ == "__main__":` block.
# It ensures that the main() function is called when the script is executed.
if (sys.nframe() == 0) {
    main()
}