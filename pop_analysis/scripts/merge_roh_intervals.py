import pandas as pd
import pybedtools

# Define input and output files
bed_file = snakemake.input["rg"]
output_bed = snakemake.output["merged_rg"]
output_filtered = snakemake.output["filtered_roh"]

columns = ["chrom", "start", "end", "individual"]
df = pd.read_csv(
    bed_file,
    sep='\t', 
    header=0, 
    names=columns, 
    dtype={"start": int, "end": int}
)

# Step 2: Filter intervals larger than 500,000 bp
df = df[df["end"] - df["start"] > 500000]
df = df.sort_values(by=["chrom", "start", "end"])

df[["chrom", "start", "end"]].to_csv(output_filtered, sep='\t', header=False, index=False)

bed = pybedtools.BedTool(output_filtered)
merged_intervals = bed.merge()

merged_intervals.saveas(output_bed)