import pandas as pd
import statistics
import gzip
import tempfile
import psutil
import os

def calc_roh(roh, fai, output, population):
    dffai = pd.read_table(fai, sep='\t', header=None)
    dffai = dffai[dffai[1] > 10000000]  # only include chromosomes gr 1Mbp
    chroms = dffai[0].values

    glength = dffai[1].sum()

    dfroh = pd.read_table(roh, sep='\t', header=0)
    dfroh = dfroh[dfroh['[6]Length (bp)'] > 500000]  # only look at 500kbp ROH and longer
    dfroh = dfroh[dfroh['[3]Chromosome'].isin(chroms)]
    print('hi')

    dffroh = dfroh.groupby(
        ['[2]Sample']
    ).agg(
        {
            '[6]Length (bp)': "sum"
        }
    ).div(glength)
    dffroh = dffroh.reset_index()

    if dffroh.empty:
        print("No individuals found with ROH greater than 500kbp.")
        with open(output, 'w'):
            pass
        with open(output.replace(".froh", "_top.froh"), "w"):
            pass
        return False

    dffroh.to_csv(output, sep='\t', index=False, header=None)
    print(dffroh)

    if len(dffroh) < 8:
        combined_df = dffroh
    else:
        # get the largest 2 ROH
        df_lg = dffroh.nlargest(2, '[6]Length (bp)')
        # get 8 random rows so we have 10 individuals total
        random_rows = dffroh.sample(n=8)
        combined_df = pd.concat([df_lg, random_rows], axis=0)
    combined_df.to_csv(output.replace(".froh", "_top.froh"), sep='\t', index=False, header=False, columns=[combined_df.columns[0]])
    #combined_df.to_csv(output, sep='\t', index=False, header=False, columns=[combined_df.columns[0]])
    return True


def main():
    roh_input = snakemake.input["rg"]
    fai_input = snakemake.input["fai"]
    froh_output = snakemake.output["froh"]
    population = snakemake.params["population"]
    interval_width = snakemake.params["interval_width"]

    calc_roh(roh_input, fai_input, froh_output, population)

if __name__ == "__main__":
    main()