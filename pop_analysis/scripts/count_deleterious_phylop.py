from cyvcf2 import VCF
import gzip

# Input and output file paths from Snakemake
vcf_file = snakemake.input["vcf"]
bedgraph_file = snakemake.input["bedgraph"]
sites = snakemake.output["sites"]
counts_output = snakemake.output["counts"]
af_output = snakemake.output["af"]
#af_output = "test.af"
"""
Reads a gzipped bedGraph file line by line, filtering positions where 
column 4 (PhyloP score) > 0 and the position is present in the VCF.
Writes (CHROM, POS, PhyloP, "DELETERIOUS") to the output file.
"""
# Step 1: Open VCF and extract positions for filtering
output_file=sites
vcf = VCF(vcf_file)
vcf_positions = set((variant.CHROM, variant.POS) for variant in vcf)
#vcf_chroms = set(variant.CHROM for variant in vcf)

with gzip.open(bedgraph_file, 'rt') as bg, open(output_file, 'w') as out:
    for line in bg:
        chrom, start, _, phylop = line.strip().split()
        phylop_score = float(phylop)
        pos = int(start) + 1 #assuming this column is 0 based 
        #print(chrom)
        #if chrom in vcf_chroms:
        #            print(f"Found matching CHROM: {chrom}")
        #print(vcf_positions)
        if phylop_score > 1 and (chrom, pos) in vcf_positions:
                #print(pos)
        #if (chrom, pos) in vcf_positions:
                out.write(f"{chrom}\t{pos}\t{phylop_score}\tDELETERIOUS\n")

file=sites
deleterious_sites = set()
with open(file, 'r') as f:
    for line in f:
        chrom, pos, _, _ = line.strip().split('\t')
        deleterious_sites.add((chrom, int(pos)))

vcf = VCF(vcf_file)
sample_names = vcf.samples
homo_alt_counts = {sample: 0 for sample in sample_names}
# Open AF output file for writing
with open(af_output, "w") as af_out:

    for variant in vcf:

        if (variant.CHROM, variant.POS) in deleterious_sites:
            # Initialize counters for allele frequencies
            non_ref_count = 0
            valid_samples = 0

            for genotype in variant.genotypes:
                if genotype[0] != -1 and genotype[1] != -1:  # Not missing
                    valid_samples += 1
                    non_ref_count += (genotype[0]) + (genotype[1])

            if valid_samples > 0:
                # Calculate non-reference allele frequency (AF)
                af = non_ref_count / (2 * valid_samples)
                af_out.write(f"{variant.CHROM}\t{variant.POS}\t{af:.5f}\n")

            # Count homozygous alternate (1/1) genotypes for individuals
            for i, genotype in enumerate(variant.genotypes):
                if genotype[0] == genotype[1] == 1:
                    sample_name = sample_names[i]
                    homo_alt_counts[sample_name] += 1

# Write homozygous alternate counts to the output file
with open(counts_output, "w") as f:
    f.write("Sample\tHomozygous_Alternate_Count\n")
    for sample, count in homo_alt_counts.items():
        f.write(f"{sample}\t{count}\n")

# if __name__ == "__main__":

#     # Step 2: Filter bedGraph to keep only deleterious sites present in the VCF
#     filter_bedgraph_to_deleterious(bedgraph_file, sites)

#     # Step 3: Load deleterious sites into a set
#     deleterious_sites = load_deleterious_sites(sites)

#     # Step 4: Count homozygous alternate sites and calculate non-reference AF
#     count_homozygous_alternate_and_af(deleterious_sites, counts_output, af_output)