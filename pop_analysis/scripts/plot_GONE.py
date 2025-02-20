import pandas as pd
import matplotlib.pyplot as plt

def plot(gone_file, output_png, ccgp_id, population_id):
    df = pd.read_csv(gone_file, sep='\t', skiprows=2, header=None)
    df.columns = ["Generation", "Ne"]
    df['Ne'] = pd.to_numeric(df['Ne'])
    
    plt.figure(figsize=(10, 6))
    plt.plot(df["Generation"], df["Ne"], linestyle='-')
    plt.xlabel("Generations")
    plt.ylabel("Ne")
    plt.title(f"Effective Population Size (Ne) over Generations - {ccgp_id} - {population_id}")
    plt.grid(True)
    
    # Save the plot as a PNG file
    plt.savefig(output_png, format='png')
    plt.close()






def main():
    gone_output = snakemake.input["GONE_out_Ne"]
    output_png = snakemake.output["GONE_plot"]
    ccgp_id = snakemake.params["project_id"]
    population_id = snakemake.params["population_id"]


    plot(gone_output, output_png, ccgp_id, population_id)

if __name__ == "__main__":
    main()