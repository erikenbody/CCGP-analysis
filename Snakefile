from snakemake.utils import min_version
min_version("7.0")
configfile: config["config_file"]
from pathlib import Path

def get_scaffolds(fai):
    scaffolds = []
    with open(fai, 'r') as fai_file:
        for line in fai_file:
            scaffold = line.split('\t')[0]
            length = line.split('\t')[1]
            if int(length) >= 500000:
                scaffolds.append(scaffold)
    return scaffolds

fai = Path("results", config["refgenome"], "data", "genome", config["refgenome"] + ".fna.fai")
pops = Path("results", config["refgenome"], "algatr", config["final_prefix"] + "_populations")
contigs = Path("results", config["refgenome"], "algatr", config["final_prefix"] + "_contigs.tsv")
output_contigs = Path("results", config["refgenome"], "algatr", config["final_prefix"] + "_chromosome_names.txt")
#output_path = Path("results", config["refgenome"], "data", "genome", "scaffolds")

output = [
        expand("results/{refGenome}/CCGP/{prefix}_roh_pi_ridges.pdf",refGenome=config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/trackhub/index.html", refGenome=config['refgenome']),
        expand("results/{refGenome}/CCGP/{prefix}_dist_done.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_pruned_mil.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        
        # Local PCA outputs, in progress:
        # expand("results/{refGenome}/local_pca/{prefix}_done.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/local_pca/bcf_subsets/{scaff}/{prefix}_{scaff}.bcf", scaff=extract_scaffolds("results/{refGenome}/data/genome/{refGenome}.fna.fai"), refGenome = config['refgenome'], prefix=config['final_prefix'])
        # expand("results/{refGenome}/local_pca/bcf_subsets/{scaff}/{prefix}_{scaff}.bcf", scaff=extract_scaffolds("/scratch2/erik/CCGP-reruns/projects/63A-Quercus/results/GCA_001633185.5/data/genome/GCA_001633185.5.fna.fai"), refGenome = config['refgenome'], prefix=config['final_prefix']),
        
        # feems output, in progress
        #expand("results/{refGenome}/feems/{prefix}_pruned_0.6_cv.png", refGenome = config['refgenome'], prefix=config['final_prefix']),
        
        # Algatr
        expand("results/{refGenome}/algatr/{prefix}.coords.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}_TESS_qmatrix.csv", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}_TESS_pieplot.pdf", refGenome = config['refgenome'], prefix=config['final_prefix']),

        # phasing output, in progress:
        #expand("results/{refGenome}/algatr/{prefix}_shapeit5_phased.bcf", refGenome = config['refgenome'], prefix=config['final_prefix']),
        #expand("results/{refGenome}/algatr/subsets/{prefix}_{scaff}_annotated.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        #expand("results/{refGenome}/algatr/haplotypes/{prefix}_{scaff}_phased.bcf", refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        
        # RDA output:
        #expand("results/{refGenome}/algatr/subsets/RDA/{scaff}/{prefix}_RDA_cortest_full.csv",refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        
        # RDA LD pruned output:
         
        # # GONE stuff.
       
        # expand("results/{refGenome}/GONE/Linux/{prefix}_plink.map", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/GONE/Linux/{prefix}_plink.ped", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/GONE/Linux/Output_d2_{prefix}_plink", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/GONE/Linux/Output_Ne_{prefix}_plink", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/GONE/Linux/OUTPUT_{prefix}_plink", refGenome = config['refgenome'], prefix=config['final_prefix']),
        
        expand("results/{refGenome}/pop_analysis/{prefix}_k.done", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/pop_analysis/{prefix}_admixture_final.pdf", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # admixture
        expand("results/{refGenome}/algatr/admixture/logs/{prefix}_best_K.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}_admixture_composite.pdf", refGenome = config['refgenome'], prefix=config['final_prefix']),

        #not using eval admix, although it is implemented
        #expand("results/{refGenome}/algatr/admixture/evaladmix/{prefix}_complete.txt", refGenome = config['refgenome'], prefix=config['final_prefix'])


        
]
if config["rename_contigs"]:
    output.append(
        expand("results/{refGenome}/algatr/{prefix}_renamed_clean_snps.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/pop_analysis/{prefix}_chromosome_names.txt", refGenome = config['refgenome'], prefix=config['final_prefix'])
    )


rule all:
    input: output

include: "ccgp/Snakefile"
#include: "reference/Snakefile"
include: "local_pca/Snakefile"
include: "alagtr/Snakefile"
include: "pop_analysis/Snakefile"
#include: "../ccgp_feems/Snakefile" #need to get envs set up for this to work

if config['bed']:
   include: "postprocess/Snakefile"
   include: "trackhub/Snakefile" 
else:
   include: "postprocess_variant/Snakefile"
