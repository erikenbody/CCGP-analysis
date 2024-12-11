input_files <- snakemake@input[["zscore_output"]]
output_file <- snakemake@output[["zscore_output"]]
# Read the first file and write its header and content to the output file
first_file <- readLines(input_files[1])
header <- first_file[1]
writeLines(first_file, output_file)

# Append the content of the remaining files, skipping their headers
for (file in input_files[-1]) {
  lines <- readLines(file)
  cat(lines[-1], file = output_file, sep = "\n", append = TRUE)
}