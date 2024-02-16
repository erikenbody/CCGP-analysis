#!/usr/bin/env bash
# use snakemake conda
bcftools index -t /scratch2/erik/CCGP-reruns/projects/3-Setophaga/results/GCA_024362935.1/3-Setophaga_raw.vcf.gz
bcftools index -t /scratch2/erik/CCGP-reruns/projects/25-Enhydra/results/GCA_002288905.2/25-Enhydra_raw.vcf.gz
bcftools index -t /scratch2/erik/CCGP-reruns/projects/46-Haliotis/results/GCA_022045235.1/46-Haliotis_raw.vcf.gz
bcftools index -t /scratch2/erik/CCGP-reruns/projects/61-Agelaius/results/GCA_023055355.1/61-Agelaius_raw.vcf.gz
bcftools index -t /scratch2/erik/CCGP-reruns/projects/78-Hetaerina/results/GCA_022747635.1/78-Hetaerina_raw.vcf.gz
bcftools index -t /scratch2/erik/CCGP-reruns/projects/97-Lynx/results/GCA_022079265.1/97-Lynx_raw.vcf.gz
#go rename haiotis after download