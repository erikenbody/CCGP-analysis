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
def get_chromosome_lists(vcf_path):
    """
    Quickly reads the VCF header, returning a dictionary of chromosomes > 1Mb 
    and a list of chromosomes < 1Mb.
    """
    long_chroms = {}
    short_chroms = []
    with pysam.VariantFile(vcf_path, 'r') as vcf_reader:
        for contig in vcf_reader.header.contigs.values():
            if contig.length >= 5_000_000:
                long_chroms[contig.name] = contig.length
            else:
                short_chroms.append(contig.name)
    
    if not long_chroms:
        sys.exit("Error: No chromosomes > 5Mb found in VCF header. Exiting.")

    print(f"Found {len(long_chroms)} chromosomes > 5Mb and {len(short_chroms)} chromosomes < 5Mb.")
    return long_chroms, short_chroms

def process_raw_data(data_path, chrom_lengths):
    """
    Processes the raw data extracted by bcftools to calculate depth and heterozygosity.
    """
    col_names = ['Chromosome', 'Sample', 'Depth', 'Genotype']
    col_types = {'Chromosome': 'category', 'Sample': 'category', 'Depth': 'Int32', 'Genotype': str}

    df = pd.read_csv(
        data_path, sep='\t', header=None, names=col_names,
        dtype=col_types, na_values=['.']
    ).dropna()
    
    # Filter for chromosomes we are analyzing (>= 1Mb)
    df = df[df['Chromosome'].isin(chrom_lengths.keys())].copy()
    df['is_het'] = df['Genotype'].str.contains('0[/|]1|1[/|]0', regex=True).astype(int)

    print("Aggregating data per sample and chromosome...")
    summary = df.groupby(['Sample', 'Chromosome'], observed=False).agg(
        Total_Depth=('Depth', 'sum'),
        Variant_Count=('Depth', 'size'),
        Het_Count=('is_het', 'sum')
    ).reset_index()

    summary['Depth'] = summary['Total_Depth'] / summary['Variant_Count']
    
    # Map chromosome lengths (in bases)
    summary['chr_len'] = summary['Chromosome'].map(chrom_lengths)
    
    initial_rows = len(summary)
    summary.dropna(subset=['chr_len'], inplace=True)
    final_rows = len(summary)
    if initial_rows > final_rows:
        print(f"Warning: Removed {initial_rows - final_rows} rows corresponding to chromosomes smaller than 1Mb that were not filtered initially.")

    # --- UPDATED: Calculate heterozygosity as a rate per base ---
    summary['Heterozygosity'] = summary['Het_Count'] / summary['chr_len']
    
    # --- UPDATED: Simpler median calculation for normalization ---
    print("Calculating medians for normalization using chromosomes >= 1Mb...")
    
    # The 'summary' DataFrame already only contains long chromosomes.
    # We can directly calculate the medians from it.
    median_stats = summary.groupby('Sample', observed=False).agg(
        Depth_SampleMedian=('Depth', 'median'),
        Heterozygosity_SampleMedian=('Heterozygosity', 'median')
    ).reset_index()
    
    # Merge the new medians back into the main summary
    summary = summary.merge(median_stats, on='Sample', how='left')

    bad_samples_mask = summary['Depth_SampleMedian'].isna() | summary['Heterozygosity_SampleMedian'].isna()
    if bad_samples_mask.any():
        bad_sample_names = summary.loc[bad_samples_mask, 'Sample'].unique()
        print(f"Warning: Removing {len(bad_sample_names)} sample(s) with failed median calculation: {list(bad_sample_names)}")
        summary = summary[~summary['Sample'].isin(bad_sample_names)]

    epsilon = 1e-9
    summary['Normalized_Depth'] = summary['Depth'] / (summary['Depth_SampleMedian'] + epsilon)
    summary['Normalized_Heterozygosity'] = summary['Heterozygosity'] / (summary['Heterozygosity_SampleMedian'] + epsilon)

    final_df = summary[['Sample', 'Chromosome', 'Depth', 'Heterozygosity', 'Normalized_Depth', 'Normalized_Heterozygosity']]
    final_df.to_csv(output_csv, index=False)
    return final_df

# --- Main Execution ---
reference_lengths, small_chromosomes = get_chromosome_lists(vcf_path)
df = process_raw_data(raw_data_path, reference_lengths)
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
grouped_chroms = df.groupby('Chromosome', observed=False)

for chrom, group in grouped_chroms:
    group_copy = group.copy()
    group_copy.replace([np.inf, -np.inf], np.nan, inplace=True)
    group_copy.dropna(subset=['Normalized_Depth'], inplace=True)

    if len(group_copy['Sample'].unique()) < 20:
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
    significant_results.to_csv(report_csv, index=False)
    significant_chromosomes = set(significant_results['Chromosome'])
else:
    pd.DataFrame(columns=['Chromosome', 'p_value', 'sil_score', 'p_adj']).to_csv(report_csv, index=False)
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
    # --- UPDATED: Y-axis label ---
    plt.ylabel('Heterozygosity Rate (Het sites per base)')
    plt.title('Normalized Depth vs. Heterozygosity')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(title='Chromosome', loc='best', markerscale=1.5)
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()

print("Script finished successfully.")
