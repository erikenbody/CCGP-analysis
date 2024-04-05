import os
import pandas as pd

def generate_population_output(sample_dict, total_populations, ccgp_project_id, ref_genome, trashed_samples):
    for key, value in sample_dict.items():
        output = f"results/{ref_genome}/algatr/{ccgp_project_id}_populations/population_{key}.txt"
        os.makedirs(os.path.dirname(output), exist_ok=True)
        with open(output, 'w') as file:
            for sample in value:
                file.write(f'{sample}\n')

    output_done = f"results/{ref_genome}/algatr/{ccgp_project_id}_populations-done.txt"
    
    # output_trash = f"results/{ref_genome}/algatr/{ccgp_project_id}_populations/trashed_samples.txt"
    # os.makedirs(os.path.dirname(output_trash), exist_ok=True)
    # with open(output_trash, 'w') as file:
    #     for sample in trashed_samples:
    #         file.write(f'{sample}\n')

    with open(output_done, 'w') as file:
        file.write(f"k_total = {total_populations}\n")
        file.write('\n')
        for key, value in sample_dict.items():
            file.write(f'@{key}\n')
            for sample in value:
                file.write(f'{sample}\n')
            file.write('\n')

    return output_done

    

def process_df(df_matrix, total_populations):

    for i, row in df_matrix.iterrows():
        total_K = row["total_K"]
        if total_K != total_populations:
            df_matrix.drop(i, axis=0, inplace=True)

    return df_matrix, total_populations

def assign_pop(df_matrix, total_populations, ccgp_project_id, ref_genome):
   
    sample_dict = {}
    trashed_samples = []

    for i, row in df_matrix.iterrows():
        pop_num = row["pop_assignment"]
        sample = row["individual"]
        qval = row["qvalue"]

        # if qval <= 0.80:
        #     sample.append(trashed_samples)
        #     break

        if pop_num in sample_dict.keys():
            if sample not in sample_dict[pop_num]:
                sample_dict[pop_num].append(sample)
            # if sample not in sample_dict[pop_num]:
            #     if qval >= 0.80:
            #         sample_dict[pop_num].append(sample)
            #     else:
            #         trashed_samples.append(sample)
        else:
            sample_dict[pop_num] = [sample]

    output_file = generate_population_output(sample_dict, total_populations, ccgp_project_id, ref_genome, trashed_samples)
    
    return output_file

def main():
    #snakemake stuff
    ccgp_project_id = snakemake.params["project_id"]
    ref_genome = snakemake.params["ref_genome"]
    manual_k = snakemake.params["manual_k_assignment"]

    matrix_file = snakemake.input["matrix"]

    df = pd.read_csv(matrix_file)
    # if type(manual_k) == str and manual_k.upper() == "NA":
    #     total_pops = df.at[1, 'best_k'] # change this to a snakemake param later
    if type(manual_k) == int:
        total_pops = manual_k
    else:
        total_pops = df.at[1, 'best_k']

    df_matrix, total_populations = process_df(df, total_pops)
    output_file = assign_pop(df_matrix, total_populations, ccgp_project_id, ref_genome)

if __name__ == "__main__":
    main()