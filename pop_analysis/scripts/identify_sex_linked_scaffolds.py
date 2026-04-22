import os
import sys
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import gzip
import gc
import psutil

# --- Snakemake I/O ---
raw_data_path = snakemake.input.raw_data
vcf_path = snakemake.input.vcf_path_for_lengths
output_csv = snakemake.output.data
report_csv = snakemake.output.report
plot_file = snakemake.output.plot
exclude_file = snakemake.output.exclude_scaffolds

# --- Directory and Logging Setup ---
output_dir = os.path.dirname(report_csv)
os.makedirs(output_dir, exist_ok=True)
print(f"Reading raw data from: {raw_data_path}")
print(f"Reading VCF header from: {vcf_path}")

# --- Helper Functions ---
def get_chromosome_lists(vcf_path, exclude_file, report_csv, output_csv, plot_file):
    """
    Quickly reads the VCF header, returning a dictionary of chromosomes > 5Mb (ref_chroms),
    a dictionary of all chromosomes > 1Mb (long_chroms), and a list of chromosomes < 1Mb.
    """
    long_chroms = {}
    short_chroms = []
    ref_chroms = {}
    with pysam.VariantFile(vcf_path, 'r') as vcf_reader:
        for contig in vcf_reader.header.contigs.values():
            if contig.length: # Ensure length is not None
                if contig.length >= 5_000_000:
                    ref_chroms[contig.name] = contig.length
                if contig.length >= 1_000_000:
                    long_chroms[contig.name] = contig.length
                else:
                    short_chroms.append(contig.name)

    if not long_chroms:
        print("No chromosomes > 1Mb found in VCF header.")
        print("Writing empty output files and exiting successfully.")
        # Write empty exclude file (no chromosomes to exclude since none are large enough)
        with open(exclude_file, 'w') as f:
            pass  # Empty file
        # Write empty report
        pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj', 'chr_len']).to_csv(report_csv, index=False)
        # Write empty data output
        pd.DataFrame(columns=['Sample', 'Chromosome', 'Depth', 'Heterozygosity', 'Normalized_Depth', 'Normalized_Heterozygosity']).to_csv(output_csv, index=False)
        # Write empty plot file
        with open(plot_file, 'w') as f:
            pass  # Empty file
        sys.exit(0)

    if not ref_chroms:
        print("No chromosomes > 5Mb found in VCF header — assembly is too fragmented for sex chromosome normalization.")
        print("Writing empty output files and exiting successfully.")
        with open(exclude_file, 'w') as f:
            pass  # Empty file — no sex chromosomes identified, nothing to exclude
        pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj', 'chr_len']).to_csv(report_csv, index=False)
        pd.DataFrame(columns=['Sample', 'Chromosome', 'Depth', 'Heterozygosity', 'Normalized_Depth', 'Normalized_Heterozygosity']).to_csv(output_csv, index=False)
        with open(plot_file, 'w') as f:
            pass  # Empty file
        sys.exit(0)

    print(f"Found {len(ref_chroms)} reference chromosomes (> 5Mb) for normalization.")
    print(f"Found {len(long_chroms)} total chromosomes > 1Mb for analysis.")
    print(f"Found {len(short_chroms)} chromosomes < 1Mb to be excluded.")
    return ref_chroms, long_chroms, short_chroms

def check_memory_usage(threshold_percent=80):
    """
    Check current memory usage and warn if it exceeds threshold.
    Returns True if memory usage is acceptable, False if critical.
    """
    memory = psutil.virtual_memory()
    percent_used = memory.percent

    if percent_used > threshold_percent:
        print(f"WARNING: High memory usage detected: {percent_used:.1f}% of {memory.total / (1024**3):.1f} GB")
        print(f"  Available: {memory.available / (1024**3):.1f} GB")
        gc.collect()  # Force garbage collection
        memory = psutil.virtual_memory()
        print(f"  After GC: {memory.percent:.1f}% used, {memory.available / (1024**3):.1f} GB available")

        if memory.percent > 90:
            print("CRITICAL: Memory usage above 90%! Stopping to prevent crash.")
            return False

    return True

def process_chromosome_data(chrom_data, chrom_name):
    """
    Process data for a single chromosome and return aggregated summary.
    This is called for each chromosome separately to minimize memory usage.
    """
    if not chrom_data:
        return pd.DataFrame()

    # Convert list of dicts to DataFrame
    df = pd.DataFrame(chrom_data)

    # Calculate heterozygosity
    df['is_het'] = df['Genotype'].str.contains('0[/|]1|1[/|]0', regex=True).astype(int)

    # Aggregate by sample for this chromosome
    summary = df.groupby('Sample').agg(
        Total_Depth=('Depth', 'sum'),
        Variant_Count=('Depth', 'size'),
        Het_Count=('is_het', 'sum')
    ).reset_index()

    summary['Chromosome'] = chrom_name

    # Explicitly delete the DataFrame and force garbage collection
    del df
    gc.collect()

    return summary

# --- MODIFICATION 1: Updated the function signature ---
def process_raw_data(data_path, all_chrom_lengths, ref_chrom_keys):
    """
    Processes raw data to calculate depth and heterozygosity.
    Normalization medians are calculated ONLY from reference chromosomes (>5Mb).
    The analysis is performed on all chromosomes > 1Mb.
    Uses chromosome-by-chromosome processing to minimize memory usage.
    """
    print(f"Processing large file chromosome-by-chromosome to minimize memory usage...")
    print(f"Analyzing {len(all_chrom_lengths)} chromosomes > 1Mb")

    # Get list of chromosomes to process
    chromosomes_to_process = set(all_chrom_lengths.keys())

    # Process file chromosome by chromosome
    print("Reading file and processing by chromosome...")
    chrom_summaries = []
    current_chunk_data = []
    current_chunk_chrom = None
    warned_large_chrom = False  # Track if we've already warned about large chromosome

    with gzip.open(data_path, 'rt') as f:
        for line_num, line in enumerate(f, 1):
            if line_num % 10_000_000 == 0:
                print(f"  Processed {line_num:,} lines...")

            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue

            chrom = parts[0]

            # Skip chromosomes we don't need
            if chrom not in chromosomes_to_process:
                # If we were collecting data and chromosome changed, process it
                if current_chunk_data and current_chunk_chrom and current_chunk_chrom in chromosomes_to_process:
                    chrom_summaries.append(process_chromosome_data(current_chunk_data, current_chunk_chrom))
                    current_chunk_data = []
                current_chunk_chrom = None
                warned_large_chrom = False  # Reset warning flag
                continue

            # If chromosome changed and we have data, process previous chromosome
            if current_chunk_chrom and chrom != current_chunk_chrom:
                if current_chunk_data:
                    print(f"  Finished chromosome {current_chunk_chrom}: {len(current_chunk_data):,} variants")
                    chrom_summaries.append(process_chromosome_data(current_chunk_data, current_chunk_chrom))
                    current_chunk_data = []

                    # Check memory usage after processing each chromosome
                    if not check_memory_usage(threshold_percent=80):
                        sys.exit("ERROR: Memory usage too high. Exiting to prevent crash.")

                warned_large_chrom = False  # Reset warning flag for new chromosome

            current_chunk_chrom = chrom

            # Safety check: if a single chromosome has too many variants, warn ONCE and check memory
            # This prevents a single huge chromosome from consuming all memory
            MAX_VARIANTS_PER_CHROM = 50_000_000  # 50 million variants per chromosome max
            if len(current_chunk_data) > MAX_VARIANTS_PER_CHROM and not warned_large_chrom:
                print(f"WARNING: Chromosome {current_chunk_chrom} has over {MAX_VARIANTS_PER_CHROM:,} variants!")
                print(f"  This may indicate an issue. Processing anyway but monitoring memory closely...")
                if not check_memory_usage(threshold_percent=75):
                    print(f"ERROR: Memory critical while processing large chromosome {current_chunk_chrom}")
                    sys.exit("Exiting to prevent crash.")
                warned_large_chrom = True  # Set flag so we don't warn again for this chromosome

            # Parse and add to current chromosome data
            try:
                sample = parts[1]
                depth = int(parts[2]) if parts[2] != '.' else None
                genotype = parts[3]

                if depth is not None:
                    current_chunk_data.append({
                        'Sample': sample,
                        'Depth': depth,
                        'Genotype': genotype
                    })
            except (ValueError, IndexError):
                continue

    # Process final chromosome if any data remains
    if current_chunk_data and current_chunk_chrom and current_chunk_chrom in chromosomes_to_process:
        print(f"  Finished chromosome {current_chunk_chrom}: {len(current_chunk_data):,} variants")
        chrom_summaries.append(process_chromosome_data(current_chunk_data, current_chunk_chrom))

    print(f"Combining results from {len(chrom_summaries)} chromosomes...")

    # Combine all chromosome summaries
    if not chrom_summaries:
        print("Warning: No valid data found in file")
        return pd.DataFrame()

    summary = pd.concat(chrom_summaries, ignore_index=True)

    summary['Depth'] = summary['Total_Depth'] / summary['Variant_Count']
    summary['chr_len'] = summary['Chromosome'].map(all_chrom_lengths)
    summary.dropna(subset=['chr_len'], inplace=True)
    summary['Heterozygosity'] = summary['Het_Count'] / summary['chr_len']
    
    # --- MODIFICATION 2: Calculate medians using only the reference chromosomes ---
    print("Calculating normalization medians using reference chromosomes (>= 5Mb)...")
    
    # Create a temporary DataFrame containing only the reference chromosomes for median calculation
    ref_summary = summary[summary['Chromosome'].isin(ref_chrom_keys)]

    median_stats = ref_summary.groupby('Sample', observed=True).agg(
        Depth_SampleMedian=('Depth', 'median'),
        Heterozygosity_SampleMedian=('Heterozygosity', 'median')
    ).reset_index()
    
    # Merge the reference-based medians back into the main summary DataFrame
    summary = summary.merge(median_stats, on='Sample', how='left')

    bad_samples_mask = summary['Depth_SampleMedian'].isna() | summary['Heterozygosity_SampleMedian'].isna()
    if bad_samples_mask.any():
        bad_sample_names = summary.loc[bad_samples_mask, 'Sample'].unique()
        print(f"Warning: Removing {len(bad_sample_names)} sample(s) with failed median calculation: {list(bad_sample_names)}")
        summary = summary[~summary['Sample'].isin(bad_sample_names)]

    epsilon = 1e-9
    summary['Normalized_Depth'] = summary['Depth'] / (summary['Depth_SampleMedian'] + epsilon)
    summary['Normalized_Heterozygosity'] = summary['Heterozygosity'] / (summary['Heterozygosity_SampleMedian'] + epsilon)

    # Replace any remaining inf or nan values that might have been created
    summary.replace([np.inf, -np.inf], np.nan, inplace=True)

    final_df = summary[['Sample', 'Chromosome', 'Depth', 'Heterozygosity', 'Normalized_Depth', 'Normalized_Heterozygosity']]
    final_df.to_csv(output_csv, index=False)
    return final_df

# --- Main Execution ---
ref_chroms, long_chroms, small_chromosomes = get_chromosome_lists(vcf_path, exclude_file, report_csv, output_csv, plot_file)

# --- MODIFICATION 3: Pass the ref_chroms dictionary keys to the processing function ---
df = process_raw_data(raw_data_path, long_chroms, ref_chroms.keys())

if df.empty:
    pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj']).to_csv(report_csv, index=False)
    with open(exclude_file, 'w') as f:
        for chrom in sorted(small_chromosomes):
            f.write(f"{chrom}\n")
    sys.exit("The processed dataframe is empty. Wrote empty reports and exiting.")

df['Cluster'] = np.nan

# --- Identify Potential Sex Chromosomes (Clustering and Stats) ---
print("Identifying sex chromosome candidates...")
results = []
grouped_chroms = df.groupby('Chromosome', observed=True)

for chrom, group in grouped_chroms:
    group_copy = group.copy()
    group_copy.replace([np.inf, -np.inf], np.nan, inplace=True)
    group_copy.dropna(subset=['Normalized_Depth', 'Normalized_Heterozygosity'], inplace=True)

    if len(group_copy['Sample'].unique()) < 20:
        continue

    # Ensure no NaN values remain before KMeans
    if group_copy['Normalized_Depth'].isna().any() or len(group_copy) == 0:
        continue

    kmeans = KMeans(n_clusters=2, random_state=42, n_init='auto')
    group_copy['Cluster'] = kmeans.fit_predict(group_copy[['Normalized_Depth']])
    
    df.loc[group_copy.index, 'Cluster'] = group_copy['Cluster']
    
    cluster_counts = group_copy['Cluster'].value_counts()
    sil_score = silhouette_score(group_copy[['Normalized_Depth']], group_copy['Cluster'])

    if all(count >= 10 for count in cluster_counts) and sil_score > 0.80:
        cluster0_het = group_copy[group_copy['Cluster'] == 0]['Normalized_Heterozygosity']
        cluster1_het = group_copy[group_copy['Cluster'] == 1]['Normalized_Heterozygosity']
        
        if len(cluster0_het.dropna()) > 1 and len(cluster1_het.dropna()) > 1:
            t_stat, p_val = ttest_ind(cluster0_het, cluster1_het, equal_var=False, nan_policy='omit')
            if not np.isnan(p_val):
                results.append([chrom, p_val, sil_score])

# --- Statistical Analysis and Reporting ---
if results:
    results_df = pd.DataFrame(results, columns=['Chromosome', 'p_value', "sil_score"])
    results_df['p_adj'] = multipletests(results_df['p_value'], method='bonferroni')[1]
    significant_results = results_df[results_df['p_adj'] < 0.0001].sort_values(by='p_adj')
    
    # --- ADDED LINE ---
    # Map the chromosome lengths from the long_chroms dictionary to the results
    significant_results['chr_len'] = significant_results['Chromosome'].map(long_chroms)
    
    significant_results.to_csv(report_csv, index=False)
    significant_chromosomes = set(significant_results['Chromosome'])
else:
    # Updated to include the new column in the empty file case
    pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj', 'chr_len']).to_csv(report_csv, index=False)
    significant_chromosomes = set()

print(f"Found {len(significant_chromosomes)} significant chromosome(s): {significant_chromosomes or 'None'}")
df['Cluster'] = df['Cluster'].astype('category')

# --- Write Exclude File ---
chromosomes_to_exclude = set(small_chromosomes).union(significant_chromosomes)
print(f"Writing {len(chromosomes_to_exclude)} chromosomes to exclude file: {exclude_file}")
with open(exclude_file, 'w') as f:
    for chrom in sorted(list(chromosomes_to_exclude)):
        f.write(f"{chrom}\n")

# --- Plotting ---
print("Generating plots...")
for chrom in significant_chromosomes:
    subset = df[df['Chromosome'] == chrom]
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=subset, x='Normalized_Depth', y='Normalized_Heterozygosity',
        hue='Cluster', palette='coolwarm', alpha=0.8, s=50
    )
    plt.xlabel('Normalized Depth (vs. Autosomal Median)')
    plt.ylabel('Normalized Heterozygosity Rate')
    plt.title(f'Sex Chromosome Candidate: {chrom}')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Depth Cluster')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'plot_{chrom}_details.png'), dpi=300)
    plt.close()

if not df.empty:
    plt.figure(figsize=(12, 8))
    autosomes = df[~df['Chromosome'].isin(significant_chromosomes)]
    if not autosomes.empty:
        plt.scatter(autosomes['Normalized_Depth'], autosomes['Heterozygosity'], color='grey', alpha=0.3, label='Autosomes')

    palette = sns.color_palette('bright', n_colors=len(significant_chromosomes))
    for i, chrom in enumerate(sorted(list(significant_chromosomes))):
        subset = df[df['Chromosome'] == chrom]
        if not subset.empty:
            plt.scatter(subset['Normalized_Depth'], subset['Heterozygosity'], color=palette[i], alpha=0.9, label=chrom)

    plt.xlabel('Normalized Depth (vs. Autosomal Median)')
    plt.ylabel('Heterozygosity Rate (Het sites per base)')
    plt.title('Normalized Depth vs. Heterozygosity')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Chromosome', loc='best', markerscale=1.5)
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()

print("Script finished successfully.")