from cyvcf2 import VCF
import pandas as pd

# Input and output file paths from Snakemake
vcf_file = snakemake.input["vcf"]
annotations_file = snakemake.input["annotations_file"]
output_file = snakemake.output["counts"]

# Load annotations and filter for deleterious variants
annotations = pd.read_csv(annotations_file, sep="\t")
deleterious = annotations[annotations['SIFT_PREDICTION'] == 'DELETERIOUS']

# Create a set of (CHROM, POS) tuples for quick lookup
deleterious_sites = set(zip(deleterious['CHROM'], deleterious['POS']))

# Open the VCF file
vcf = VCF(vcf_file)

# Initialize a dictionary to track homozygous 1/1 counts per individual
sample_names = vcf.samples  # List of sample names from the VCF
homo_alt_counts = {sample: 0 for sample in sample_names}

# Iterate through variants in the VCF
for variant in vcf:
    if (variant.CHROM, variant.POS) in deleterious_sites:
        for i, genotype in enumerate(variant.genotypes):
            # Check if genotype is homozygous alternate (1/1)
            if genotype[0] == genotype[1] == 1:
                sample_name = sample_names[i]
                homo_alt_counts[sample_name] += 1

# Write results to a tab-delimited output file
with open(output_file, "w") as f:
    f.write("Sample\tHomozygous_Alternate_Count\n")
    for sample, count in homo_alt_counts.items():
        f.write(f"{sample}\t{count}\n")