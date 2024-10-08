import sys
import shutil

'''
This script converts scaffold names in a `.map` file to numeric chromosome identifiers based on their order in an input file. It creates a conversion dictionary from the input file, updates the `.map` file using this dictionary, and writes the result to an output file.
'''


def generate_mapping(fai_file, map_file, output_file):

    # conversion_dict = {}
    # chrom_num = 1
    # with open(fai_file, 'r') as f:
    #     for line in f:
    #         line = line.strip().split()
    #         conversion_dict[line[0]] = str(chrom_num)
    #         chrom_num += 1

    entries = []
    with open(fai_file, 'r') as f:
        for line in f:
            line = line.strip().split()
            entries.append((line[0], int(line[1])))
    entries.sort(key=lambda x: x[1], reverse=True)
    conversion_dict = {entry[0]: str(idx + 1) for idx, entry in enumerate(entries)}
    #for contig, number in conversion_dict.items():
        #print(f"{contig}: {number}")
    # read bim file and replace the scaffold names with numbering 1:n (n = number of scaffolds)
    updated_lines = []
    with open(map_file, 'r') as f:
        for line in f:
            elements = line.strip().split('\t')
            scaffold = elements[0]
            if scaffold in conversion_dict:
                elements[0] = conversion_dict[scaffold]
            updated_lines.append('\t'.join(elements))

    with open(output_file, 'w') as f:
        for line in updated_lines:
            f.write(line + '\n')

fai_file = snakemake.params["fai"]
map_file = snakemake.input["map_temp"]
output_file = snakemake.output["mapp"]
generate_mapping(fai_file, map_file, output_file)