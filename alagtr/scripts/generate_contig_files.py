import argparse
import os


def create_files(contig, refgenome, prefix):
    filename = f"results/{refgenome}/data/genome/scaffolds/{prefix}_{contig}.txt"
    # cwd = os.getcwd()
    # print(cwd)
    with open(filename, 'w') as scaffold_file:
        scaffold_file.write(contig)
        print(f'made {contig}')
    

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--contig-id", type=str, help="Specify contig ID")
    parser.add_argument("--refgenome", type=str, help="Specify ref genome")
    parser.add_argument("--prefix", type=str, help="Specify project ID")
    args = parser.parse_args()

    if args.contig_id:
        create_files(args.contig_id, args.refgenome, args.prefix)


if __name__ == "__main__":
    main()