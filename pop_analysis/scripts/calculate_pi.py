import allel
import csv
import pandas as pd
import numpy as np
import pyranges as pr

def process_and_augment_with_fai_per_sample(roh_file, fai_file):
    # Step 1: Process the ROH file into a dictionary organized by sample
    sample_intervals = {}
    with open(roh_file, 'r') as f:
        for line in f:
            chrom, start, end, sample = line.strip().split('\t')
            if sample not in sample_intervals:
                sample_intervals[sample] = []
            sample_intervals[sample].append({
                'chrom': chrom,
                'start': int(start),
                'end': int(end),
                'state': 1  # Inbred
            })

    # Step 2: Load the .fai file and determine chromosome order
    fai_chrom_order = []
    with open(fai_file, 'r') as f:
        for line in f:
            chrom = line.split('\t')[0]  # Extract chromosome name
            fai_chrom_order.append(chrom)

    # Step 3: Augment intervals for each sample
    augmented_intervals_by_sample = {}
    for sample, intervals in sample_intervals.items():
        augmented_intervals = []
        existing_chroms = {interval['chrom'] for interval in intervals}

        # Iterate through the chromosomes in the order of the .fai file
        for chrom in fai_chrom_order:
            chrom_intervals = [i for i in intervals if i['chrom'] == chrom]
            chrom_intervals.sort(key=lambda x: x['start'])

            # Add outbred intervals between inbred intervals
            prev_end = 0
            for interval in chrom_intervals:
                if interval['start'] > prev_end:
                    augmented_intervals.append({
                        'chrom': chrom,
                        'pos': prev_end,
                        'sample': sample,
                        'state': 0  # Outbred
                    })
                augmented_intervals.append({
                    'chrom': chrom,
                    'pos': interval['start'],
                    'sample': sample,
                    'state': 1  # Inbred
                })
                prev_end = interval['end'] + 1

            # Add an outbred interval after the last inbred interval
            augmented_intervals.append({
                'chrom': chrom,
                'pos': prev_end,
                'sample': sample,
                'state': 0  # Outbred
            })

        # Add outbred intervals for chromosomes not present in this sample
        for chrom in fai_chrom_order:
            if chrom not in existing_chroms:
                augmented_intervals.append({
                    'chrom': chrom,
                    'pos': 0,
                    'sample': sample,
                    'state': 0  # Outbred
                })

        # Sort by chromosome order and position
        chrom_rank = {chrom: rank for rank, chrom in enumerate(fai_chrom_order)}
        augmented_intervals.sort(key=lambda x: (chrom_rank[x['chrom']], x['pos']))

        augmented_intervals_by_sample[sample] = augmented_intervals

    return augmented_intervals_by_sample

def process_genotypes_with_augmented_intervals(vcf_path, augmented_intervals, genome_length):
    import allel
    import numpy as np

    # Read the VCF file
    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    
    # Extract sample names and variant positions
    samples = callset['samples']
    vchrom = callset['variants/CHROM']
    pos = callset['variants/POS']
    
    results = {}  # Dictionary to store n1, inbred length, and other metrics for each sample

    # Iterate through samples
    for sample_idx, sample_name in enumerate(samples):
        # Fetch the intervals for this sample
        intervals = augmented_intervals.get(sample_name, [])
        
        if not intervals:
            results[sample_name] = {
                "heterozygous_count": 0,
                "inbred_length": 0,
                "total_length": genome_length,
                "new_length": genome_length,
                "heterozygosity_metric": 0
            }
            continue
        
        # Sort intervals by chrom and pos
        intervals.sort(key=lambda x: (x['chrom'], x['pos']))

        interval_idx = 0  # Start with the first interval
        current_interval = intervals[interval_idx]

        n1 = 0  # Heterozygous count for this sample
        inbred_length = 0  # Total length of inbred regions

        # Iterate through variants
        for chrom_name, variant_pos in zip(vchrom, pos):
            # Advance the interval index until we find the relevant interval
            while interval_idx < len(intervals) - 1 and (
                current_interval['chrom'] != chrom_name or variant_pos >= intervals[interval_idx + 1]['pos']
            ):
                interval_idx += 1
                current_interval = intervals[interval_idx]
            
            # Determine the state for the current variant
            if current_interval['chrom'] == chrom_name:
                state = current_interval['state']
                
                # Check if the current position is in an outbred region (state = 0)
                if state == 0:
                    #print(current_interval)
                    variant_idx = np.where(pos == variant_pos)[0][0]
                    geno = genotypes[variant_idx, sample_idx]
                    #print(chrom_name, variant_pos, geno)
                    # Count heterozygous sites
                    if geno[0] != geno[1]:  # Heterozygous
                        n1 += 1

        # Calculate inbred length for this sample
        for i in range(len(intervals) - 1):
            if intervals[i]['state'] == 1:
                inbred_length += intervals[i + 1]['pos'] - intervals[i]['pos']

        # Compute new length (outbred length)
        new_length = genome_length - inbred_length

        # Calculate heterozygosity metric
        heterozygosity_metric = n1 / new_length if new_length > 0 else 0

        results[sample_name] = {
            "heterozygous_count": n1,
            "inbred_length": inbred_length,
            "total_length": genome_length,
            "new_length": new_length,
            "heterozygosity_metric": heterozygosity_metric
        }
    
    return results

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
    return pd.DataFrame({'sample_name': sample_names, 'avg_heterozygosity': avg_het, "het_counts": het_counts})

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
    #process_genotypes_with_roh(vcf_path, roh_file)
    #augmented_intervals=process_and_augment_with_fai_per_sample(roh_file, fai)
    
    #adjusted_het = process_genotypes_with_augmented_intervals(vcf_path, augmented_intervals, total_length)

    # with open(out_het_not_roh, 'w', newline='') as f:
    #     writer = csv.writer(f, delimiter='\t')
    #     writer.writerow(["Sample", "Heterozygous_Count", "Inbred_Length", "Total_Length", "New_Length", "Heterozygosity_Metric"])
    #     for sample, metrics in adjusted_het.items():
    #         writer.writerow([
    #             sample,
    #             metrics["heterozygous_count"],
    #             metrics["inbred_length"],
    #             metrics["total_length"],
    #             metrics["new_length"],
    #             metrics["heterozygosity_metric"]
    #         ])

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
    fai = snakemake.input["fai"]
    k_value = snakemake.params["k_value"]
    ccgp_project_id = snakemake.params["project_id"]
    roh_file = snakemake.input["roh_file"]
    out_pi = snakemake.output["out_pi"]
    #out_pi_not_roh = snakemake.output["out_pi_not_roh"]
    #out_het_not_roh = snakemake.output["out_het_not_roh"]
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
