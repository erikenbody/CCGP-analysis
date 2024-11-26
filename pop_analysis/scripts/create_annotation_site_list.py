# import pandas as pd
# import sys

# # Accept arguments from Snakemake
# input_file = snakemake.input.annotations
# output_files = snakemake.output.pos

# # Mapping of variant types to their filtering logic
# variant_filters = {
#     "deleterious": {"column": "SIFT_PREDICTION", "value": "DELETERIOUS"},
#     "tolerated": {"column": "SIFT_PREDICTION", "value": "TOLERATED"},
#     "synonymous": {"column": "VARIANT_TYPE", "value": "SYNONYMOUS"},
#     "nonsynonymous": {"column": "VARIANT_TYPE", "value": "NONSYNONYMOUS"}
# }

# def filter_and_save(df, variant_type, output_path):
#     """Filter the DataFrame based on variant type and save the positions."""
#     filter_info = variant_filters[variant_type]
#     filtered_df = df[df[filter_info["column"]] == filter_info["value"]]
#     filtered_df[["CHROM", "POS"]].to_csv(output_path, sep="\t", index=False, header=False)

# def main():
#     """Main function to load data, filter, and create output files."""
#     # Load the annotation file
#     df = pd.read_csv(input_file, sep="\t")

#     # Loop through the variant types and save their corresponding files
#     for variant_type, output_file in zip(variant_filters.keys(), output_files):
#         filter_and_save(df, variant_type, output_file)

# # Run the script
# if __name__ == "__main__":
#     main()

import pandas as pd
from pathlib import Path

# Input and output paths from Snakemake
annotations_file = snakemake.input["annotations"]
output_dir = Path(snakemake.output[0]).parent  # Directory for position files

# Read the annotations file
annotations = pd.read_csv(annotations_file, sep="\t")

# Define variant types to extract
variant_types = ["DELETERIOUS", "TOLERATED"]

# Create position files for each variant type
for var_type in variant_types:
    positions = annotations[annotations['SIFT_PREDICTION'].str.upper() == var_type]
    output_file = output_dir / f"{var_type.lower()}_positions.txt"
    
    # Write positions to file in VCFTOOLS-compatible format (CHROM and POS only)
    positions[['CHROM', 'POS']].to_csv(output_file, sep="\t", index=False, header=False)

# Define variant types to extract
variant_types = ["SYNONYMOUS", "NONSYNONYMOUS"]

# Create position files for each variant type
for var_type in variant_types:
    positions = annotations[annotations['VARIANT_TYPE'].str.upper() == var_type]
    output_file = output_dir / f"{var_type.lower()}_positions.txt"
    
    # Write positions to file in VCFTOOLS-compatible format (CHROM and POS only)
    positions[['CHROM', 'POS']].to_csv(output_file, sep="\t", index=False, header=False)