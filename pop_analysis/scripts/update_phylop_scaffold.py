import pandas as pd
import gzip
import argparse

def load_mapping(mapping_file):
    """Load the mapping of old CHROM names to new CHROM names."""
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            old_chrom, new_chrom = line.strip().split()
            mapping[new_chrom] = old_chrom  # Corrected order
    return mapping

def update_chrom(data_file, mapping, output_file):
    """Update the CHROM column in a gzipped bedGraph file based on the mapping."""
    with gzip.open(data_file, 'rt') as f:
        df = pd.read_csv(f, sep='\t', header=None, names=["CHROM", "START", "END", "SCORE"])

    # Replace CHROM values using the mapping dictionary
    df['CHROM'] = df['CHROM'].map(mapping)

    # Save the updated DataFrame to a gzipped output file
    with gzip.open(output_file, 'wt') as out:
        df.to_csv(out, sep='\t', index=False, header=False)

def main():
    parser = argparse.ArgumentParser(description="Update CHROM values in a bedGraph file based on a mapping file.")
    parser.add_argument("mapping_file", help="Path to the file with old-to-new CHROM mappings.")
    parser.add_argument("data_file", help="Path to the gzipped bedGraph file with original CHROM column.")
    parser.add_argument("output_file", help="Path to save the updated gzipped bedGraph file.")

    args = parser.parse_args()

    # Load mapping and update the CHROM values
    mapping = load_mapping(args.mapping_file)
    update_chrom(args.data_file, mapping, args.output_file)

if __name__ == "__main__":
    main()