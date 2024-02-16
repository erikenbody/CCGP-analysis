from snakemake.utils import min_version
min_version("7.0")
configfile: config["config_file"]
from pathlib import Path

def get_scaffolds(fai):
    scaffolds = []
    with open(fai, 'r') as fai_file:
        for line in fai_file:
            scaffold = line.split('\t')[0]
            scaffolds.append(scaffold)
    return scaffolds


fai = Path("results", config["refgenome"], "data", "genome", config["refgenome"] + ".fna.fai")
#output_path = Path("results", config["refgenome"], "data", "genome", "scaffolds")

rule all:
    input:
        expand("results/{refGenome}/CCGP/{prefix}_roh_pi_ridges.pdf",refGenome=config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_filtered.froh",refGenome=config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/{prefix}_clean_snps.vcf.gz", refGenome=config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/{prefix}_clean_indels.vcf.gz", refGenome=config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_annotated_pruned_0.3_dosage.txt.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_annotated_pruned_0.6_dosage.txt.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_annotated_no_pruning_dosage.txt.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        #expand("results/{refGenome}/trackhub/trackDb.txt", refGenome=config['refgenome']), #original trackhub output, no longer relavent in 2024
        expand("results/{refGenome}/trackhub/index.html", refGenome=config['refgenome']),
        expand("results/{refGenome}/CCGP/{prefix}_dist_done.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}_TESS_qmatrix.csv", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/algatr/{prefix}_RDA_outliers_best_Zscores.csv",refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_pruned_mil.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        #expand("results/{refGenome}/CCGP/{prefix}_annotated_pruned_0.6.dist.id", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/local_pca/{prefix}_done.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/local_pca/bcf_subsets/{scaff}/{prefix}_{scaff}.bcf", scaff=extract_scaffolds("results/{refGenome}/data/genome/{refGenome}.fna.fai"), refGenome = config['refgenome'], prefix=config['final_prefix'])
        # expand("results/{refGenome}/local_pca/bcf_subsets/{scaff}/{prefix}_{scaff}.bcf", scaff=extract_scaffolds("/scratch2/erik/CCGP-reruns/projects/63A-Quercus/results/GCA_001633185.5/data/genome/GCA_001633185.5.fna.fai"), refGenome = config['refgenome'], prefix=config['final_prefix']),
        #expand("results/{refGenome}/feems/{prefix}_pruned_0.6_cv.png", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}.coords.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_annotated.vcf.gz.csi", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}.samps4coords.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        #expand("results/{refGenome}/data/genome/scaffolds/{prefix}_{scaff}.txt", refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        expand("results/{refGenome}/algatr/subsets/{prefix}_{scaff}_annotated.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        expand("results/{refGenome}/algatr/subsets/RDA/{scaff}/{prefix}_RDA_anova_best.csv",refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai))



include: "ccgp/Snakefile"
include: "reference/Snakefile"
include: "local_pca/Snakefile"
include: "alagtr/Snakefile"
#include: "../ccgp_feems/Snakefile" #need to get envs set up for this to work

if config['bed']:
   include: "postprocess/Snakefile"
   include: "trackhub/Snakefile" #need to add the output
else:
   include: "postprocess_variant/Snakefile"

