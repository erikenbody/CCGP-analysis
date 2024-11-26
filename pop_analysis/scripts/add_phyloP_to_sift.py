import sys
import gzip

def bedgraph_generator(bedgraph_path):
    """Generator that yields (chrom, pos) -> value from the bedGraph file."""
    with gzip.open(bedgraph_path, 'rt') as f:
        for line in f:
            chrom, start, _, value = line.strip().split()
            pos = int(start) + 1  # Convert to 1-based
            yield (chrom, pos), float(value)

def create_bedgraph_cache(bedgraph_path):
    """Create a cache for quick lookup from the bedGraph file."""
    cache = {}
    for key, value in bedgraph_generator(bedgraph_path):
        cache[key] = value
    return cache

def process_files(sift_file, bedgraph_file, output_file="SIFT_merged.tsv"):
    # Create the bedGraph cache for lookup
    bedgraph_cache = create_bedgraph_cache(bedgraph_file)

    with open(sift_file, 'r') as sf, open(output_file, 'w') as out:
        # Write the header with the new 'PhyloP' column
        header = sf.readline().strip() + "\tPhyloP\n"
        out.write(header)

        # Process each line of the SIFT file
        for line in sf:
            fields = line.strip().split('\t')
            chrom = fields[-1]  # CHROM.bak
            pos = int(fields[1])  # POS

            # Lookup the PhyloP score, defaulting to None if not found
            phyloP = bedgraph_cache.get((chrom, pos), None)

            # Write the updated line with the PhyloP score
            out.write(line.strip() + f"\t{phyloP}\n")

    print(f"Merged file saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python merge_sift_bedgraph.py <SIFT_file> <bedGraph_file>")
        sys.exit(1)

    sift_file = sys.argv[1]
    bedgraph_file = sys.argv[2]
    process_files(sift_file, bedgraph_file)