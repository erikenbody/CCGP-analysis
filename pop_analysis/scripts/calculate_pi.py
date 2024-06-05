import allel
import csv

def calculate_mean_pairwise_diff(vcf_path):

    callset = allel.read_vcf(vcf_path)
    genotypes = allel.GenotypeArray(callset['calldata/GT'])
    ac = genotypes.count_alleles()

    # Calculate mean pairwise difference for all variant sites
    pairwise_diff = allel.mean_pairwise_difference(ac, an=None, fill=0)
    sum_pairwise_diff = sum(pairwise_diff)
    return sum_pairwise_diff

def load_bed_file(bed_path):
    # calculate the total length of callable sequence based on a bed file of callable regions
    
    total_length = 0

    with open(bed_path, 'r') as bed_file:
        for line in bed_file:
            fields = line.strip().split()
            chrom, start, end = fields[:3]
            total_length += int(end) - int(start)

    return total_length

def main(vcf_path, bed_path, k_value, ccgp_project_id, ref_genome):

    sum_pairwise_diff = calculate_mean_pairwise_diff(vcf_path)

    total_length = load_bed_file(bed_path)

    pi = sum_pairwise_diff / total_length

    output_path = f"results/{ref_genome}/pop_analysis/{ccgp_project_id}_population_pi/pi_population_{k_value}.csv"
    with open(output_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([ccgp_project_id, k_value, pi])



    print(sum_pairwise_diff)
    print(total_length)
    print(pi)

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Calculate ratio of mean pairwise difference to total length of BED file")
    # parser.add_argument("vcf_path", help="Path to the VCF file")
    # parser.add_argument("bed_path", help="Path to the BED file")
    vcf_input = snakemake.input["pop_vcf"]
    bed_input = snakemake.input["bed"]
    k_value = snakemake.params["k_value"]
    ccgp_project_id = snakemake.params["project_id"]
    ref_genome = snakemake.params["ref_genome"]
    #args = parser.parse_args()
    main(vcf_input, bed_input, k_value, ccgp_project_id, ref_genome)
