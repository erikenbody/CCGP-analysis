import pandas as pd
import argparse
import subprocess




def grab_chrom_names(file, ccgp_project_id, contig_names, rename_bool):
    if rename_bool:
        df = pd.read_csv(file, sep='\t')
        with open(contig_names, 'w') as f:
            for i, row in df.iterrows():
                sequence_name = row["sequence_name"]
                accession = row["accession"]
                line = f"{accession}\t{sequence_name}\n"
                f.write(line)

    else:
         with open(contig_names, 'w') as f:
            pass


def main():

    ccgp_project_id = snakemake.params["project_id"]
    #rename_bool = snakemake.params["rename_contigs"]
    rename_bool = True
    contig_file = snakemake.input["contigs"]
    contig_names = snakemake.output["renamed"]

    grab_chrom_names(contig_file, ccgp_project_id, contig_names, rename_bool)
    

if __name__ == "__main__":
    main()