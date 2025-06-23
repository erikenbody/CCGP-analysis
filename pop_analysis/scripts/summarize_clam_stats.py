import pandas as pd

pi_out, het_out, pi_out2, het_out2, summary, k = snakemake.input.pi_out, snakemake.input.het_out, snakemake.input.pi_out2, snakemake.input.het_out2, snakemake.output.summary, snakemake.wildcards.k

def avg_pi(path):
    df = pd.read_csv(path, sep="\t")
    return df['pi'].mean()

def avg_het(path):
    df = pd.read_csv(path, sep="\t")
    # Compute per-row heterozygosity
    df['het'] = df['count_het_sites'] / (df['callable_bases'])
    # Average per sample, then across samples
    return df.groupby('sample_name')['het'].mean().mean()

avg_pi_outside_roh = avg_pi(pi_out)
avg_pi_genome = avg_pi(pi_out2)
avg_het_outside_roh = avg_het(het_out)
avg_het_genome = avg_het(het_out2)

with open(summary, "w") as f:
    f.write("pop\tavg_pi_genome\tavg_het_genome\tavg_pi_outside_roh\tavg_het_outside_roh\n")
    f.write(f"{k}\t{avg_pi_genome}\t{avg_het_genome}\t{avg_pi_outside_roh}\t{avg_het_outside_roh}\n")