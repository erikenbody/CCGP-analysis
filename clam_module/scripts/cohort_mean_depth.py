"""Write a bash-sourceable .env file with MIN_MEAN_DEPTH and MAX_MEAN_DEPTH for
`clam loci`.

Two modes (selected by config; resolved by clamg_cohort_mean_depth in the
parent Snakefile):

  flat: snakemake.params.flat_min_mean_depth and flat_max_mean_depth are both
        set => write those values through verbatim. No QC inputs required.

  estimate: estimate cohort per-base mean depth from QC bam_sumstats + .fai:
        per_sample_depth ~ (Total_Reads * Percent_mapped/100 - Num_duplicates)
                            * 150 / sum(scaffold lengths)
        cohort_mean = mean across samples in samples.tsv
        bounds      = mean +/- sd_scale * sqrt(mean)
        150 bp = CCGP PE150 nominal read length. Mosdepth's default flag mask
        excludes DUP, so dups are subtracted to approximate what snpArcher's
        cov_filter actually counted.
"""

import math
from pathlib import Path


def _read_samples(path):
    wanted = set()
    with open(path) as fh:
        for line in fh:
            s = line.split("\t")[0].strip()
            if s:
                wanted.add(s)
    return wanted


def _genome_length(fai_path):
    n = 0
    with open(fai_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                n += int(parts[1])
    return n


def _per_sample_depths(sumstats_path, wanted, read_length, genome_len):
    depths = []
    seen = set()
    with open(sumstats_path) as fh:
        fh.readline()  # Sample Total_Reads Percent_mapped Num_duplicates ...
        for line in fh:
            fields = [c.strip() for c in line.rstrip("\n").split("\t")]
            if len(fields) < 4:
                continue
            name = fields[0]
            if name not in wanted:
                continue
            seen.add(name)
            total_reads = float(fields[1])
            pct_mapped = float(fields[2])
            num_dups = float(fields[3])
            mapped_unique = max(0.0, total_reads * pct_mapped / 100.0 - num_dups)
            depths.append(mapped_unique * read_length / genome_len)
    return depths, seen


def _write_flat(out_path, lo, hi):
    Path(out_path).write_text(
        f"# Source: flat config (clam_min_mean_depth / clam_max_mean_depth)\n"
        f"SOURCE=flat_config\n"
        f"MIN_MEAN_DEPTH={lo:.6f}\n"
        f"MAX_MEAN_DEPTH={hi:.6f}\n"
    )


def _write_estimate(out_path, *, cohort_mean, sd, sd_scale, lo, hi,
                    n_samples, genome_len, read_length, sumstats_path):
    Path(out_path).write_text(
        f"# Source: estimated from {sumstats_path}\n"
        f"# n_samples={n_samples} genome_length={genome_len} read_length={read_length}\n"
        f"# cohort_mean_depth={cohort_mean:.6f} poisson_sd={sd:.6f} sd_scale={sd_scale}\n"
        f"SOURCE=estimated_bam_sumstats\n"
        f"COHORT_MEAN_DEPTH={cohort_mean:.6f}\n"
        f"POISSON_SD={sd:.6f}\n"
        f"SD_SCALE={sd_scale}\n"
        f"MIN_MEAN_DEPTH={lo:.6f}\n"
        f"MAX_MEAN_DEPTH={hi:.6f}\n"
        f"N_SAMPLES={n_samples}\n"
        f"GENOME_LENGTH={genome_len}\n"
        f"READ_LENGTH={read_length}\n"
    )


def main(snakemake):
    out_path = snakemake.output.depth_file
    flat_min = snakemake.params.flat_min_mean_depth
    flat_max = snakemake.params.flat_max_mean_depth

    if flat_min is not None and flat_max is not None:
        _write_flat(out_path, float(flat_min), float(flat_max))
        return

    samples_tsv = snakemake.input.samples_tsv
    bam_sumstats = snakemake.input.bam_sumstats
    fai = snakemake.input.fai
    sd_scale = float(snakemake.params.sd_scale)
    read_length = int(snakemake.params.read_length)

    wanted = _read_samples(samples_tsv)
    genome_len = _genome_length(fai)
    if genome_len == 0:
        raise ValueError(f"Empty genome length from {fai}")

    depths, seen = _per_sample_depths(bam_sumstats, wanted, read_length, genome_len)

    missing = wanted - seen
    if missing:
        raise ValueError(
            f"{len(missing)} samples in samples.tsv not found in "
            f"{bam_sumstats}: {sorted(missing)[:5]}..."
        )
    if not depths:
        raise ValueError("No samples matched between samples.tsv and bam_sumstats")

    cohort_mean = sum(depths) / len(depths)
    sd = math.sqrt(cohort_mean)
    lo = max(0.0, cohort_mean - sd_scale * sd)
    hi = cohort_mean + sd_scale * sd

    _write_estimate(
        out_path,
        cohort_mean=cohort_mean, sd=sd, sd_scale=sd_scale, lo=lo, hi=hi,
        n_samples=len(depths), genome_len=genome_len, read_length=read_length,
        sumstats_path=bam_sumstats,
    )


main(snakemake)  # noqa: F821
