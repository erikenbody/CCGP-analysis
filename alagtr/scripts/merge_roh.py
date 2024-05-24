import pandas as pd
import pybedtools
import tempfile
import os

# input_file = snakemake.input["rg"]
# output_file = snakemake.output["merged_rg"]
# gap_size = snakemake.params["gap_size"]

def merge_rohs(input_file, output_file, gap_size):
    df = pd.read_csv(input_file, sep="\t", header=None, skiprows=1, dtype={0: str, 1: str, 2: str, 3: int, 4: int, 5: float, 6: float, 7: float})

    df_reordered = df[[2, 3, 4, 1]]

    sample_list = []
    for i, row in df_reordered.iterrows():
        sample = row[1]
        if sample not in sample_list:
            sample_list.append(sample)

    with open(output_file, "w") as outfile: 
        for sample in sample_list:
            # print(df_reordered)
            filtered_df = df_reordered[df_reordered[1] == sample] 

            temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False)  # Create a temporary file for filtered dataframe.
            filtered_df.to_csv(temp_file.name, sep="\t", header=False, index=False)  # Export filtered DataFrame to temporary bed file
            temp_file.close()

            bed_file = pybedtools.BedTool(temp_file.name)
            #sorted_bed = bed_file.sort()
            #merged_bed = sorted_bed.merge(c=[4], o='distinct', d=5000)
            merged_bed = bed_file.merge(c=[4], o='distinct', d=gap_size)

            for interval in merged_bed:
                outfile.write(str(interval))

            os.unlink(temp_file.name)
            df_reordered = df_reordered[df_reordered[1] != sample]



def main():
    input_file = snakemake.input["rg"]
    output_file = snakemake.output["merged_rg"]
    gap_size = snakemake.params["gap_size"]


    merge_rohs(input_file, output_file, gap_size)


if __name__ == "__main__":
    main()

    