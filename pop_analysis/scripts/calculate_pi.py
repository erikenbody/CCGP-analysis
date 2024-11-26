import allel
import csv
import pandas as pd
import numpy as np

def calculate_mean_pairwise_diff_not_roh(vcf_path, roh_file):
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    
    # Extract sample names and variant positions
    samples = callset['samples']
    chrom = callset['variants/CHROM']
    pos = callset['variants/POS']

    # Read the ROH file into a DataFrame
    roh_df = pd.read_csv(roh_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'sample'])

    # Group ROH intervals by sample for bulk processing
    sample_groups = roh_df.groupby('sample')

    # Iterate through samples
    for sample, group in sample_groups:
        sample_index = np.where(samples == sample)[0]

        if sample_index.size == 0:
            print(f"Sample {sample} not found in genotypes columns")
            continue

        sample_index = sample_index[0]

        # Create a mask for all intervals of the sample in one operation
        mask = np.zeros(len(pos), dtype=bool)
        for _, row in group.iterrows():
            roh_chrom, start, end = row['chrom'], row['start'], row['end']
            mask |= (chrom == str(roh_chrom)) & (pos >= start) & (pos <= end)

        # Apply the mask to the genotype array
        genotypes[mask, sample_index] = -1

    # Calculate allele counts and mean pairwise differences
    ac = genotypes.count_alleles()
    pairwise_diff = allel.mean_pairwise_difference(ac, an=None, fill=0)
    sum_pairwise_diff = sum(pairwise_diff)

    return sum_pairwise_diff

def calculate_mean_pairwise_diff(vcf_path):
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    ac = genotypes.count_alleles()

    pairwise_diff = allel.mean_pairwise_difference(ac, an=None, fill=0)
    sum_pairwise_diff = sum(pairwise_diff)
    return sum_pairwise_diff

def load_bed_file(bed_path):
    # calculate the total length of callable sequence based on a bed file of callable regions
    total_length = 0

    with open(bed_path, 'r') as bed_file:
        for line in bed_file:
            fields = line.strip().split()
            chrom, start, end = fields[:3]
            total_length += int(end) - int(start)

    return total_length

def main(vcf_path, bed_path, k_value, ccgp_project_id, roh_file):

    sum_pairwise_diff = calculate_mean_pairwise_diff(vcf_path)

    total_length = load_bed_file(bed_path)

    pi = sum_pairwise_diff / total_length

    with open(out_pi, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([ccgp_project_id, k_value, pi])

    # calc pi for regions not in ROH
    sum_pairwise_diff = calculate_mean_pairwise_diff_not_roh(vcf_path, roh_file)

    total_length = load_bed_file(bed_path)

    pi = sum_pairwise_diff / total_length

    with open(out_pi_not_roh, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([ccgp_project_id, k_value, pi])



    print(sum_pairwise_diff)
    print(total_length)
    print(pi)

if __name__ == "__main__":

    vcf_input = snakemake.input["pop_vcf"]
    bed_input = snakemake.input["bed"]
    k_value = snakemake.params["k_value"]
    ccgp_project_id = snakemake.params["project_id"]
    roh_file = snakemake.input["roh_file"]
    out_pi = snakemake.output["out_pi"]
    out_pi_not_roh = snakemake.output["out_pi_not_roh"]

    main(vcf_input, bed_input, k_value, ccgp_project_id, roh_file)

# Old version of the function
# def calculate_mean_pairwise_diff_not_roh(vcf_path, roh_file):
#     callset = allel.read_vcf(vcf_path)
#     genotypes = allel.GenotypeArray(callset['calldata/GT'])
#     samples = callset['samples']

#     roh_df = pd.read_csv(roh_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'sample'])

#     chrom = callset['variants/CHROM']
#     pos = callset['variants/POS']

#     for index, row in roh_df.iterrows():
#         sample = row['sample']
#         roh_chrom = str(row['chrom'])
#         start = int(row['start'])
#         end = int(row['end'])

#         # Find the indices in the genotypes object that correspond to the interval
#         mask = (chrom == roh_chrom) & (pos >= start) & (pos <= end)

#         # Find the column index for the sample
#         sample_index = np.where(samples == sample)[0]

#         if sample_index.size > 0:
#             sample_index = sample_index[0]
#             #mask using -1 because its an integer
#             genotypes[mask, sample_index] = -1
#         else:
#             print(f"Sample {sample} not found in genotypes columns")

#     # Calculate allele counts
#     ac = genotypes.count_alleles()
#     # Calculate mean pairwise difference for all variant sites
#     pairwise_diff = allel.mean_pairwise_difference(ac, an=None, fill=0)
#     sum_pairwise_diff = sum(pairwise_diff)
#     return sum_pairwise_diff
