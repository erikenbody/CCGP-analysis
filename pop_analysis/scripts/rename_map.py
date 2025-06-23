import sys
import shutil

'''
This script converts scaffold names in a `.map` file to numeric chromosome identifiers based on their order in an input file. It creates a conversion dictionary from the input file, updates the `.map` file using this dictionary, and writes the result to an output file.
'''


def generate_mapping_from_map(map_file, output_file):
    scaffold_to_chr = {}
    current_chr = 1
    updated_lines = []

    with open(map_file, 'r') as f:
        for line in f:
            elements = line.strip().split('\t')
            scaffold = elements[0]
            if scaffold not in scaffold_to_chr:
                scaffold_to_chr[scaffold] = str(current_chr)
                current_chr += 1
            elements[0] = scaffold_to_chr[scaffold]
            updated_lines.append('\t'.join(elements))

    with open(output_file, 'w') as f:
        for line in updated_lines:
            f.write(line + '\n')

map_file = snakemake.input["map_temp"]
output_file = snakemake.output["mapp"]
generate_mapping_from_map(map_file, output_file)