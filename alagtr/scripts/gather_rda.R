library(tidyverse)

#set up log file writing
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

log_smk()

zscore_input <- snakemake@output[["zscore_output"]]
rdadapt_input <- snakemake@input[["rdadapt_output"]]
cortest_input <- snakemake@input[["cortest_output"]]
#imputed_input <- snakemake@input[["imputed_output"]]
peakram_input <- snakemake@input[["peakram_output"]]
colsum_input <- snakemake@input[["colsum_output"]]
ybar_input <- snakemake@input[["ybar_output"]]
v_input <- snakemake@input[["v_output"]]
u_input <- snakemake@input[["u_output"]]
wa_input <- snakemake@input[["wa_output"]]
qr_input <- snakemake@input[["qr_output"]]
eig_input <- snakemake@input[["eig_output"]]
biplot_input <- snakemake@input[["biplot_output"]]
qraux_input <- snakemake@input[["qraux_output"]]
envcentre_input <- snakemake@input[["envcentre_output"]]
chi_input <- snakemake@input[["chi_output"]]
scalload_input <- snakemake@input[["scalload_output"]]
unscalload_input <- snakemake@input[["unscalload_output"]]

zscore_output <- snakemake@output[["zscore_output"]]
rdadapt_output <- snakemake@output[["rdadapt_output"]]
cortest_output <- snakemake@output[["cortest_output"]]
#imputed_output <- snakemake@output[["imputed_output"]]
peakram_output <- snakemake@output[["peakram_output"]]
colsum_output <- snakemake@output[["colsum_output"]]
ybar_output <- snakemake@output[["ybar_output"]]
v_output <- snakemake@output[["v_output"]]
u_output <- snakemake@output[["u_output"]]
wa_output <- snakemake@output[["wa_output"]]
qr_output <- snakemake@output[["qr_output"]]
eig_output <- snakemake@output[["eig_output"]]
biplot_output <- snakemake@output[["biplot_output"]]
qraux_output <- snakemake@output[["qraux_output"]]
envcentre_output <- snakemake@output[["envcentre_output"]]
chi_output <- snakemake@output[["chi_output"]]
scalload_output <- snakemake@output[["scalload_output"]]
unscalload_output <- snakemake@output[["unscalload_output"]]

# concat_samples <- function(input_files, output_file){
#     # Read the first file and write its header and content to the output file
#     first_file <- readLines(input_files[1])
#     header <- first_file[1]
#     writeLines(first_file, output_file)

#     # Append the content of the remaining files, skipping their headers
#     for (file in input_files[-1]) {
#         lines <- readLines(file)
#         cat(lines[-1], file = output_file, sep = "\n", append = TRUE)
#     }
# }

concat_samples <- function(input_files, output_file){
    all_data <- list()
    
    add_scaff_col <- function(data, file){
        data <- data %>% mutate(scaff = basename(dirname(file)))
        return(data)
    }

    for (file in input_files) {
        if (file.info(file)$size > 0) {
            print(file)
            data <- read_tsv(file, col_types = cols(.default = "c"))
            data <- data %>% mutate(scaff = basename(dirname(file)))
            all_data <- c(all_data, list(data))  # Use `c()` to concatenate data frames
        } else {
            message(paste("Skipping empty file:", file))
        }
    }
    
    if (length(all_data) > 0) {
        combined_data <- bind_rows(all_data, .id = "path")
        write_tsv(combined_data, output_file)
    } else {
        file.create(output_file)
    }
}

concat_samples(snakemake@input[["zscore_output"]], snakemake@output[["zscore_output"]])
concat_samples(snakemake@input[["rdadapt_output"]], snakemake@output[["rdadapt_output"]])
concat_samples(snakemake@input[["cortest_output"]], snakemake@output[["cortest_output"]])
#concat_samples(snakemake@input[["imputed_output"]], snakemake@output[["imputed_output"]])
concat_samples(snakemake@input[["peakram_output"]], snakemake@output[["peakram_output"]])
concat_samples(snakemake@input[["colsum_output"]], snakemake@output[["colsum_output"]])
#concat_samples(snakemake@input[["ybar_output"]], snakemake@output[["ybar_output"]])
concat_samples(snakemake@input[["v_output"]], snakemake@output[["v_output"]])
concat_samples(snakemake@input[["u_output"]], snakemake@output[["u_output"]])
concat_samples(snakemake@input[["wa_output"]], snakemake@output[["wa_output"]])
concat_samples(snakemake@input[["qr_output"]], snakemake@output[["qr_output"]])
concat_samples(snakemake@input[["eig_output"]], snakemake@output[["eig_output"]])
concat_samples(snakemake@input[["biplot_output"]], snakemake@output[["biplot_output"]])
concat_samples(snakemake@input[["qraux_output"]], snakemake@output[["qraux_output"]])
concat_samples(snakemake@input[["envcentre_output"]], snakemake@output[["envcentre_output"]])
concat_samples(snakemake@input[["chi_output"]], snakemake@output[["chi_output"]])
concat_samples(snakemake@input[["scalload_output"]], snakemake@output[["scalload_output"]])
concat_samples(snakemake@input[["unscalload_output"]], snakemake@output[["unscalload_output"]])