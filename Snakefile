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

def get_k(pops):
    k_values = []
    if pops.is_dir():
        files = [file.name for file in pops.iterdir() if file.is_file()]
        k_values = [file.split("_")[1].split(".")[0] for file in files]

    return k_values

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
        expand("results/{refGenome}/CCGP/{prefix}_filtered.froh",refGenome=config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/{prefix}_clean_snps.vcf.gz", refGenome=config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/{prefix}_clean_indels.vcf.gz", refGenome=config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_annotated_pruned_0.3_dosage.txt.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_annotated_pruned_0.6_dosage.txt.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_annotated_no_pruning_dosage.txt.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/trackhub/index.html", refGenome=config['refgenome']),
        expand("results/{refGenome}/CCGP/{prefix}_dist_done.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_pruned_mil.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
        
        # Local PCA outputs, in progress:
        # expand("results/{refGenome}/local_pca/{prefix}_done.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        # expand("results/{refGenome}/local_pca/bcf_subsets/{scaff}/{prefix}_{scaff}.bcf", scaff=extract_scaffolds("results/{refGenome}/data/genome/{refGenome}.fna.fai"), refGenome = config['refgenome'], prefix=config['final_prefix'])
        # expand("results/{refGenome}/local_pca/bcf_subsets/{scaff}/{prefix}_{scaff}.bcf", scaff=extract_scaffolds("/scratch2/erik/CCGP-reruns/projects/63A-Quercus/results/GCA_001633185.5/data/genome/GCA_001633185.5.fna.fai"), refGenome = config['refgenome'], prefix=config['final_prefix']),
        
        # feems output, in progress
        #expand("results/{refGenome}/feems/{prefix}_pruned_0.6_cv.png", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}.coords.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/CCGP/{prefix}_annotated.vcf.gz.csi", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}.samps4coords.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        
        # phasing output, in progress:
        #expand("results/{refGenome}/algatr/{prefix}_shapeit5_phased.bcf", refGenome = config['refgenome'], prefix=config['final_prefix']),
        #expand("results/{refGenome}/algatr/subsets/{prefix}_{scaff}_annotated.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        #expand("results/{refGenome}/algatr/haplotypes/{prefix}_{scaff}_phased.bcf", refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        
        # RDA output:
        #expand("results/{refGenome}/algatr/subsets/RDA/{scaff}/{prefix}_imputed_simple.txt",refGenome = config['refgenome'], prefix=config['final_prefix'], scaff=get_scaffolds(fai)),
        
        # RDA LD pruned output:
        #expand("results/{refGenome}/algatr/{prefix}_imputed_simple.txt",refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}_TESS_qmatrix.csv", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}_populations-done.txt", refGenome = config['refgenome'], prefix=config['final_prefix']),
        expand("results/{refGenome}/algatr/{prefix}_population_vcf/population_{k}.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix'], k=get_k(pops)),
        expand("results/{refGenome}/algatr/{prefix}_population_pi/pi_population_{k}.csv", refGenome = config['refgenome'], prefix=config['final_prefix'], k=get_k(pops)),
        expand("results/{refGenome}/algatr/{prefix}_roh/population_{k}.rg.roh", refGenome = config['refgenome'], prefix=config['final_prefix'], k=get_k(pops)),
        expand("results/{refGenome}/algatr/{prefix}_roh/population_{k}.roh.gz", refGenome = config['refgenome'], prefix=config['final_prefix'], k=get_k(pops))
        
]
if config["rename_contigs"]:
    output.append(
        expand("results/{refGenome}/algatr/{prefix}_renamed_clean_snps.vcf.gz", refGenome = config['refgenome'], prefix=config['final_prefix']),
    )


rule all:
    input: output
        


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

