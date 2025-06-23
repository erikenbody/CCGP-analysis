import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot(gone_files, output_png, output_avg_png, ccgp_id, population_id):
    plt.figure(figsize=(10, 6))

    all_dfs = []
    
    for seed, gone_file in gone_files.items():
        df = pd.read_csv(gone_file, sep='\t', skiprows=2, header=None, names=["Generation", "Ne"])
        df["Ne"] = pd.to_numeric(df["Ne"])
        df["Seed"] = seed
        all_dfs.append(df)
        
        plt.plot(df["Generation"], df["Ne"], linestyle='-', label=f'Seed {seed}')
    
    plt.xlabel("Generations")
    plt.ylabel("Ne")
    plt.title(f"Effective Population Size (Ne) over Generations - {ccgp_id} - {population_id}")
    plt.legend()
    plt.grid(True)
    
    plt.savefig(output_png, format='png')
    plt.close()

    # Combine all data
    combined_df = pd.concat(all_dfs, ignore_index=True)

    # Compute mean and confidence intervals
    summary_df = combined_df.groupby("Generation")["Ne"].agg(['mean', 'std']).reset_index()
    summary_df["lower_CI"] = summary_df["mean"] - summary_df["std"]
    summary_df["upper_CI"] = summary_df["mean"] + summary_df["std"]

    # Plot with average Ne and confidence intervals
    plt.figure(figsize=(10, 6))
    plt.fill_between(summary_df["Generation"], summary_df["lower_CI"], summary_df["upper_CI"], color="lightgray", alpha=0.5)
    plt.plot(summary_df["Generation"], summary_df["mean"], color="blue", label="Mean Ne")

    plt.xlabel("Generations")
    plt.ylabel("Average Ne")
    plt.title(f"Average Effective Population Size (Ne) over Generations - {ccgp_id} - {population_id}")
    plt.legend()
    plt.grid(True)

    plt.savefig(output_avg_png, format='png')
    plt.close()

def main():
    gone_files = {seed: path for seed, path in zip(snakemake.params["seeds"], snakemake.input["GONE_out_Ne"])}
    output_png = snakemake.output["GONE_plot"]
    output_avg_png = snakemake.output["GONE_avg_plot"]
    ccgp_id = snakemake.params["project_id"]
    population_id = snakemake.params["population_id"]
    
    plot(gone_files, output_png, output_avg_png, ccgp_id, population_id)

if __name__ == "__main__":
    main()