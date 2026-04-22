#!/usr/bin/env python3
"""
Analyze scaffolds for anomalous depth and heterozygosity patterns.

This script identifies:
1. Scaffolds with average normalized depth < 0.25
2. Scaffolds with heterozygosity close to zero (< 1e-6)
3. Provides diagnostic information about potential causes
4. Creates visualization plots

Can be run standalone or via Snakemake.

Author: Erik (with Claude Code assistance)
Date: 2025-12-19
"""

import pandas as pd
import numpy as np
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def analyze_anomalies(input_csv, vcf_path=None, depth_threshold=0.25, het_threshold=1e-6):
    """
    Analyze depth and heterozygosity data to identify anomalous scaffolds.

    Parameters:
    -----------
    input_csv : str
        Path to the input CSV file
    vcf_path : str, optional
        Path to VCF file to extract scaffold lengths
    depth_threshold : float
        Normalized depth threshold (default: 0.25)
    het_threshold : float
        Heterozygosity threshold (default: 1e-6)

    Returns:
    --------
    df : DataFrame
        Original per-sample data
    scaffold_stats : DataFrame
        Per-scaffold statistics and anomaly classifications
    """

    print(f"Reading data from: {input_csv}")
    df = pd.read_csv(input_csv)

    print(f"Loaded {len(df)} records for {df['Sample'].nunique()} samples "
          f"across {df['Chromosome'].nunique()} chromosomes")

    # Group by chromosome to calculate per-scaffold statistics
    scaffold_stats = df.groupby('Chromosome').agg({
        'Normalized_Depth': ['mean', 'std', 'min', 'max'],
        'Heterozygosity': ['mean', 'std', 'min', 'max'],
        'Sample': 'count'
    }).reset_index()

    # Flatten column names
    scaffold_stats.columns = ['Chromosome',
                              'Mean_Normalized_Depth', 'Std_Normalized_Depth',
                              'Min_Normalized_Depth', 'Max_Normalized_Depth',
                              'Mean_Heterozygosity', 'Std_Heterozygosity',
                              'Min_Heterozygosity', 'Max_Heterozygosity',
                              'Sample_Count']

    # Add scaffold lengths if VCF provided
    if vcf_path:
        print(f"Extracting scaffold lengths from: {vcf_path}")
        try:
            import pysam
            lengths = {}
            with pysam.VariantFile(vcf_path, 'r') as vcf:
                for contig in vcf.header.contigs.values():
                    if contig.length:
                        lengths[contig.name] = contig.length
            scaffold_stats['Scaffold_Length'] = scaffold_stats['Chromosome'].map(lengths)
            print(f"  Added scaffold lengths for {scaffold_stats['Scaffold_Length'].notna().sum()} scaffolds")
        except Exception as e:
            print(f"  Warning: Could not extract scaffold lengths: {e}")
            scaffold_stats['Scaffold_Length'] = np.nan
    else:
        scaffold_stats['Scaffold_Length'] = np.nan

    # Identify anomalous scaffolds
    low_depth = scaffold_stats['Mean_Normalized_Depth'] < depth_threshold
    low_het = scaffold_stats['Mean_Heterozygosity'] < het_threshold
    zero_het = scaffold_stats['Mean_Heterozygosity'] == 0

    # Create flags
    scaffold_stats['Low_Depth_Flag'] = low_depth
    scaffold_stats['Low_Het_Flag'] = low_het
    scaffold_stats['Zero_Het_Flag'] = zero_het

    # Add anomaly category
    scaffold_stats['Anomaly_Type'] = 'Normal'
    scaffold_stats.loc[low_depth & ~low_het, 'Anomaly_Type'] = 'Low_Depth_Only'
    scaffold_stats.loc[~low_depth & low_het, 'Anomaly_Type'] = 'Low_Het_Only'
    scaffold_stats.loc[low_depth & low_het, 'Anomaly_Type'] = 'Low_Depth_And_Het'
    scaffold_stats.loc[zero_het, 'Anomaly_Type'] = 'Zero_Het'

    # Calculate coefficient of variation
    scaffold_stats['CV_Depth'] = (scaffold_stats['Std_Normalized_Depth'] /
                                   (scaffold_stats['Mean_Normalized_Depth'] + 1e-9))

    # Sort
    scaffold_stats = scaffold_stats.sort_values(
        by=['Anomaly_Type', 'Mean_Normalized_Depth'],
        ascending=[False, True]
    )

    return df, scaffold_stats


def print_summary(scaffold_stats, df, depth_threshold, het_threshold, output_file=None):
    """Print summary statistics to console and optionally to file."""

    low_depth = scaffold_stats['Mean_Normalized_Depth'] < depth_threshold
    low_het = scaffold_stats['Mean_Heterozygosity'] < het_threshold

    # Build summary text
    summary_lines = []
    summary_lines.append("="*80)
    summary_lines.append("SUMMARY STATISTICS")
    summary_lines.append("="*80)
    summary_lines.append(f"\nTotal scaffolds: {len(scaffold_stats)}")
    summary_lines.append(f"Total samples: {df['Sample'].nunique()}")
    summary_lines.append(f"\nAnomaly counts:")
    summary_lines.append(scaffold_stats['Anomaly_Type'].value_counts().to_string())
    summary_lines.append(f"\nLow depth (<{depth_threshold}): {low_depth.sum()}")
    summary_lines.append(f"Low het (<{het_threshold}): {low_het.sum()}")

    anomalous = scaffold_stats[scaffold_stats['Anomaly_Type'] != 'Normal']
    if len(anomalous) > 0:
        summary_lines.append(f"\nTop 10 anomalous scaffolds:")
        cols = ['Chromosome', 'Mean_Normalized_Depth', 'Mean_Heterozygosity', 'CV_Depth']
        if 'Scaffold_Length' in anomalous.columns and anomalous['Scaffold_Length'].notna().any():
            cols.insert(1, 'Scaffold_Length')
        summary_lines.append(anomalous[cols].head(10).to_string(index=False))

    # Contamination warning
    potential_contam = scaffold_stats[scaffold_stats['Mean_Normalized_Depth'] < 0.35]
    if len(potential_contam) > 0:
        summary_lines.append(f"\n⚠️  WARNING: {len(potential_contam)} scaffolds with depth <0.35x (potential contamination)")
        summary_lines.append("\nScaffolds to validate:")
        for _, row in potential_contam.head(10).iterrows():
            summary_lines.append(f"  {row['Chromosome']}: depth={row['Mean_Normalized_Depth']:.3f}x, "
                                f"het={row['Mean_Heterozygosity']:.2e}, CV={row['CV_Depth']:.2f}")

    summary_text = "\n".join(summary_lines)

    # Print to console
    print("\n" + summary_text)

    # Write to file if specified
    if output_file:
        with open(output_file, 'w') as f:
            f.write(summary_text + "\n")
        print(f"\nSummary written to: {output_file}")


def plot_anomalies(df, scaffold_stats, output_png):
    """Create visualization."""

    low_het = scaffold_stats[scaffold_stats['Anomaly_Type'].isin(['Low_Het_Only', 'Zero_Het'])]['Chromosome'].tolist()
    print(f"\nCreating plot with {len(low_het)} anomalous scaffolds...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Left: Full data
    autosomes = df[~df['Chromosome'].isin(low_het)]
    ax1.scatter(autosomes['Normalized_Depth'], autosomes['Heterozygosity'],
               color='lightgray', alpha=0.3, s=10, label='Normal')

    if len(low_het) > 0:
        for chrom in low_het[:10]:
            subset = df[df['Chromosome'] == chrom]
            ax1.scatter(subset['Normalized_Depth'], subset['Heterozygosity'],
                       alpha=0.6, s=30, marker='^')
        ax1.scatter([], [], alpha=0.6, s=40, marker='^', color='blue',
                   label=f'Low-het (n={len(low_het)})')

    ax1.set_xlabel('Normalized Depth', fontsize=12)
    ax1.set_ylabel('Heterozygosity', fontsize=12)
    ax1.set_title('All Scaffolds', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Right: Zoom
    if len(low_het) > 0:
        ax2.scatter(autosomes['Normalized_Depth'], autosomes['Heterozygosity'],
                   color='lightgray', alpha=0.3, s=10)
        for chrom in low_het:
            subset = df[df['Chromosome'] == chrom]
            label = chrom if len(low_het) <= 5 else ''
            ax2.scatter(subset['Normalized_Depth'], subset['Heterozygosity'],
                       alpha=0.7, s=40, marker='^', label=label)

        low_het_data = df[df['Chromosome'].isin(low_het)]
        if not low_het_data.empty:
            max_het = low_het_data['Heterozygosity'].max()
            max_depth = low_het_data['Normalized_Depth'].max()
            ax2.set_xlim(-0.1, min(3.0, max_depth * 1.2))
            ax2.set_ylim(-max_het*0.1, max_het * 1.5)
        if len(low_het) <= 5:
            ax2.legend(fontsize=8)

    ax2.set_xlabel('Normalized Depth', fontsize=12)
    ax2.set_ylabel('Heterozygosity', fontsize=12)
    ax2.set_title('Anomalous Scaffolds', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved: {output_png}")


def main_snakemake(snakemake):
    """Run via Snakemake."""

    df, scaffold_stats = analyze_anomalies(
        snakemake.input.depth_het,
        vcf_path=snakemake.input.vcf
    )

    scaffold_stats.to_csv(snakemake.output.csv, index=False)
    print(f"\nWrote CSV: {snakemake.output.csv}")

    print_summary(scaffold_stats, df, 0.25, 1e-6, output_file=snakemake.output.txt)
    plot_anomalies(df, scaffold_stats, snakemake.output.png)


def main_standalone():
    """Run standalone."""

    parser = argparse.ArgumentParser(description='Analyze scaffold anomalies')
    parser.add_argument('input_csv', help='Input depth/het CSV')
    parser.add_argument('output_csv', help='Output anomaly CSV')
    parser.add_argument('output_png', help='Output plot PNG')
    parser.add_argument('--vcf', help='VCF for scaffold lengths')
    parser.add_argument('--summary-txt', help='Optional text file for summary output')
    args = parser.parse_args()

    df, scaffold_stats = analyze_anomalies(args.input_csv, vcf_path=args.vcf)
    scaffold_stats.to_csv(args.output_csv, index=False)
    print(f"\nWrote CSV: {args.output_csv}")

    print_summary(scaffold_stats, df, 0.25, 1e-6, output_file=args.summary_txt)
    plot_anomalies(df, scaffold_stats, args.output_png)


if __name__ == '__main__':
    if 'snakemake' in globals():
        main_snakemake(snakemake)
    else:
        main_standalone()
