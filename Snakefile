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

# rename_contigs = config["rename_chromosome"]
# def get_chrom_names(contigs, output_contigs):
#     if rename_contigs:
#         if contigs.is_dir():
#             df = pd.read_csv(contigs, sep='\t')
#             with open(output_contigs, 'w') as output_file:
#                 for i, row in df.iterrows():
#                     sequence_name = row["sequence_name"]
#                     accession = row["accession"]
#                     line = f"{accession}\t{sequence_name}\n"
#                     output_file.write(line)
#             return output_contigs
#     else:
#         with open(output_contigs, 'w') as output_file:
#             pass

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

        # phasing output, in progress:
        #expand("results/{refGenome}/algatr/{prefix}_shapeit5_phased.bcf", refGenome = config['refgenome'], prefix=config['final_prefix']),
        #expand("results/{refGenome}/algatr/subsets/{prefix}_{scaff}_annotated.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        #expand("results/{refGenome}/algatr/haplotypes/{prefix}_{scaff}_phased.bcf", refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        
        # RDA output:
        expand("results/{refGenome}/algatr/subsets/RDA/{scaff}/{prefix}_RDA_cortest_full.csv",refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        
        # RDA LD pruned output:
         
        # # GONE stuff.
       
        # expand("results/{refGenome}/GONE/Linux/{prefix}_plink.map", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/GONE/Linux/{prefix}_plink.ped", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/GONE/Linux/Output_d2_{prefix}_plink", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/GONE/Linux/Output_Ne_{prefix}_plink", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/GONE/Linux/OUTPUT_{prefix}_plink", refGenome = config['refgenome'], prefix=config['final_prefix']),
        
        #expand("results/{refGenome}/pop_analysis/{prefix}/k.done", refGenome = config['refgenome'], prefix=config['final_prefix']),
        #expand("results/{refGenome}/{prefix}_GONE_downlaoded.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
    
        
]
if config["rename_contigs"]:
    output.append(
        expand("results/{refGenome}/algatr/{prefix}_renamed_clean_snps.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/pop_analysis/{prefix}_chromosome_names.txt", refGenome = config['refgenome'], prefix=config['final_prefix'])
    )


rule all:
    input: output

include: "ccgp/Snakefile"
include: "reference/Snakefile"
include: "local_pca/Snakefile"
include: "alagtr/Snakefile"
include: "pop_analysis/Snakefile"
#include: "../ccgp_feems/Snakefile" #need to get envs set up for this to work

if config['bed']:
   include: "postprocess/Snakefile"
   include: "trackhub/Snakefile" 
else:
   include: "postprocess_variant/Snakefile"
