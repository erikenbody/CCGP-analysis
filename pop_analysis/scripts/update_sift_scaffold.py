import pandas as pd
import argparse

def load_mapping(mapping_file):
    """Load the mapping of new CHROM names."""
    mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            old_chrom, new_chrom = line.strip().split()
            mapping[new_chrom] = old_chrom
    return mapping

def update_chrom(data_file, mapping, output_file):
    """Update the CHROM column based on the mapping."""
    # Load the second file into a DataFrame
    df = pd.read_csv(data_file, sep='\t')
    
    # Create a backup column for original CHROM
    df['CHROM.bak'] = df['CHROM']
    
    # Replace CHROM values using the mapping dictionary
    df['CHROM'] = df['CHROM'].map(mapping)
    
    # Save the updated DataFrame to a new file
    df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Update CHROM values based on a mapping file.")
    parser.add_argument("mapping_file", help="Path to the file with old-to-new CHROM mappings.")
    parser.add_argument("data_file", help="Path to the data file with the original CHROM column.")
    parser.add_argument("output_file", help="Path to save the updated data file.")

    args = parser.parse_args()

    # Load mapping and update the CHROM values
    mapping = load_mapping(args.mapping_file)
    update_chrom(args.data_file, mapping, args.output_file)

if __name__ == "__main__":
    main()