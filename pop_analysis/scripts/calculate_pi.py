import allel
import csv
import pandas as pd
import numpy as np
import pyranges as pr

def calculate_mean_pairwise_diff_not_roh(vcf_path, roh_file):
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    
    # Extract sample names and variant positions
    samples = callset['samples']
    vchrom = callset['variants/CHROM']
    pos = callset['variants/POS']

    # Read the ROH file into a DataFrame
    #roh_df = pd.read_csv(roh_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'sample'])
    
    roh_data = {}
    with open(roh_file, 'r') as f:
        for i, line in enumerate(f):
            chrom, start, end, sample = line.strip().split('\t')
            roh_data[i] = {'chrom': chrom, 'start': int(start), 'end': int(end), 'sample': sample}

    # Create a dictionary for the next ROH index
    next_roh_index = {}
    for i in range(len(roh_data)):
        next_roh_index[i] = i + 1 if i + 1 < len(roh_data) else None
   
 # Initialize counters
    n0 = 0  # Homozygous reference
    n1 = 0  # Other genotypes
    
    # Iterate through samples
    for sample_idx, sample_name in enumerate(samples):
        # Iterate by chromosome and position
        for chrom_name, variant_pos in zip(vchrom, pos):
            # Check if the position is in any ROH region
            for roh in roh_data.values():
                if chrom_name == roh['chrom'] and roh['start'] <= variant_pos <= roh['end']:

                    # Get the genotype for the sample at this position
                    variant_idx = np.where(pos == variant_pos)[0][0]
                    geno = genotypes[variant_idx, sample_idx]              
                    # Check genotype and update counters
                    if (geno[0] == geno[1]):  # Checks if the two alleles are the same (homozygous)
                        n0 += 1
                    else:
                        n1 += 1
                    break  # Stop checking further ROH regions for this position


    # Output the results
    print(f"Homozygous reference genotypes (n0): {n0}")
    print(f"Other genotypes (n1): {n1}")


# def calculate_mean_pairwise_diff_not_roh(vcf_path, roh_file):
#     # Load VCF and initialize genotypes
#     callset = allel.read_vcf(vcf_path)
#     genotypes = allel.GenotypeArray(callset['calldata/GT'])
    
#     # Extract sample names, chromosome, and positions
#     samples = callset['samples']
#     chrom = np.array(callset['variants/CHROM'], dtype=str)  # Ensure string type for comparison
#     pos = callset['variants/POS']
    
#     # Create a DataFrame for VCF positions
#     vcf_df = pd.DataFrame({'Chromosome': chrom, 'Start': pos, 'End': pos + 1})
    
#     # Convert to pyranges object
#     vcf_ranges = pr.PyRanges(vcf_df)
    
#     # Load ROH file into a DataFrame
#     roh_df = pd.read_csv(roh_file, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Sample'])
    
#     # Create a PyRanges object for ROH intervals
#     roh_ranges = pr.PyRanges(roh_df)
    
#     # Copy the genotype array for masking
#     genotypes_filtered = genotypes.copy()
    
#     # Iterate over unique samples in the ROH file
#     for sample in roh_df['Sample'].unique():
#         # Subset ROH intervals for the current sample
#         sample_roh_ranges = roh_ranges[roh_ranges.Sample == sample]
        
#         # Find sample index in the VCF samples
#         sample_index = np.where(samples == sample)[0]
#         if sample_index.size == 0:
#             print(f"Sample {sample} not found in genotypes columns")
#             continue
#         sample_index = sample_index[0]
#         # Intersect ROH intervals with VCF positions
#         sample_roh_ranges = sample_roh_ranges.drop(like="Sample")
#         print(sample_roh_ranges)
#         #overlaps = vcf_ranges.join(sample_roh_ranges, how='right')
#         overlaps = vcf_ranges.intersect(sample_roh_ranges)
        
#         # Mask genotypes for intersecting positions
#         for _, row in overlaps.df.iterrows():
#             chrom_mask = chrom == row['Chromosome']
#             pos_mask = (pos >= row['Start']) & (pos < row['End'])
#             combined_mask = chrom_mask & pos_mask
#             genotypes_filtered[combined_mask, sample_index] = -1
    
#     # Calculate allele counts and mean pairwise differences
#     ac = genotypes_filtered.count_alleles()
#     pairwise_diff = allel.mean_pairwise_difference(ac, an=None, fill=0)
#     sum_pairwise_diff = np.sum(pairwise_diff)
    
#     return sum_pairwise_diff

# def calculate_mean_pairwise_diff_not_roh(vcf_path, roh_file):
#     callset = allel.read_vcf(vcf_path)
#     genotypes = allel.GenotypeArray(callset['calldata/GT'])
    
#     # Extract sample names and variant positions
#     samples = callset['samples']
#     chrom = callset['variants/CHROM']
#     pos = callset['variants/POS']

#     # Read the ROH file into a DataFrame
#     roh_df = pd.read_csv(roh_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'sample'])

#     # Group ROH intervals by sample for bulk processing
#     sample_groups = roh_df.groupby('sample')

#     # Iterate through samples
#     for sample, group in sample_groups:
#         sample_index = np.where(samples == sample)[0]

#         if sample_index.size == 0:
#             print(f"Sample {sample} not found in genotypes columns")
#             continue

#         sample_index = sample_index[0]

#         # Create a mask for all intervals of the sample in one operation
#         mask = np.zeros(len(pos), dtype=bool)
#         for _, row in group.iterrows():
#             roh_chrom, start, end = row['chrom'], row['start'], row['end']
#             mask |= (chrom == str(roh_chrom)) & (pos >= start) & (pos <= end)

#         # Apply the mask to the genotype array
#         genotypes[mask, sample_index] = -1

#     # Calculate allele counts and mean pairwise differences
#     ac = genotypes.count_alleles()
#     pairwise_diff = allel.mean_pairwise_difference(ac, an=None, fill=0)
#     sum_pairwise_diff = sum(pairwise_diff)

#     return sum_pairwise_diff

def calculate_mean_pairwise_diff(vcf_path):
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    ac = genotypes.count_alleles()

    pairwise_diff = allel.mean_pairwise_difference(ac, an=None, fill=0)
    sum_pairwise_diff = sum(pairwise_diff)
    return sum_pairwise_diff

def calculate_sample_heterozygosity(vcf_path, total_length):
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])

    # Compute heterozygosity (Ho) per individual
    is_het = genotypes.is_het()  # A boolean array marking heterozygous calls
    het_counts = is_het.sum(axis=0)  # Sum heterozygous calls per individual
    avg_het = het_counts / total_length  # Average per callable site

    # Get sample names
    sample_names = callset['samples']
    return pd.DataFrame({'sample_name': sample_names, 'avg_heterozygosity': avg_het})

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

    total_length = load_bed_file(bed_path)

    heterozygosity_df = calculate_sample_heterozygosity(vcf_path, total_length)
    heterozygosity_df.to_csv(out_het, index=False)

    sum_pairwise_diff = calculate_mean_pairwise_diff(vcf_path)

    pi = sum_pairwise_diff / total_length

    avg_heterozygosity = heterozygosity_df['avg_heterozygosity'].mean()

    with open(out_pi, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([ccgp_project_id, k_value, pi, avg_heterozygosity])

    # calc pi for regions not in ROH
    # sum_pairwise_diff = calculate_mean_pairwise_diff_not_roh(vcf_path, roh_file)

    # total_length = load_bed_file(bed_path)

    # pi = sum_pairwise_diff / total_length

    # with open(out_pi_not_roh, 'w', newline='') as file:
    #     writer = csv.writer(file)
    #     writer.writerow([ccgp_project_id, k_value, pi])

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
    #out_pi_not_roh = snakemake.output["out_pi_not_roh"]
    out_het = snakemake.output["out_het"]

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
