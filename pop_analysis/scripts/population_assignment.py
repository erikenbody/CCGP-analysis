import os
import pandas as pd

def generate_population_output(sample_dict, total_populations, ccgp_project_id, ref_genome, trashed_samples):
    for key, value in sample_dict.items():
        output = f"results/{ref_genome}/pop_analysis/{ccgp_project_id}_populations/population_{key}.txt"
        os.makedirs(os.path.dirname(output), exist_ok=True)
        with open(output, 'w') as file:
            for sample in value:
                file.write(f'{sample}\n')

    output_done = f"results/{ref_genome}/pop_analysis/{ccgp_project_id}_populations-done.txt"
    
    # output_trash = f"results/{ref_genome}/pop_analysis/{ccgp_project_id}_trashed_samples/trashed_samples.txt"
    # os.makedirs(os.path.dirname(output_trash), exist_ok=True)
    # with open(output_trash, 'w') as file:
    #     for sample in trashed_samples:
    #         file.write(f'{sample}\n')

    with open(output_done, 'w') as file:
        file.write(f"k_total = {total_populations}\n")
        file.write('\n')
        file.write('Excluded samples:\n')
        for samp in trashed_samples:
            file.write(f'{samp}\n')
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
        pop_num = row["pop_assignment"]
        kval = row["kval"].replace("K", "")
        if total_K != total_populations:
            df_matrix.drop(i, axis=0, inplace=True)
        else:
            if int(pop_num) != int(kval):
                df_matrix.drop(i, axis=0, inplace=True)


    return df_matrix, total_populations

def assign_pop(df_matrix, total_populations, ccgp_project_id, ref_genome, threshold):
   
    sample_dict = {}
    trashed_samples = []

    for i, row in df_matrix.iterrows():
        pop_num = row["pop_assignment"]
        kval = row["kval"].replace("K", "")
        sample = row["individual"]
        qval = row["qvalue"]
        print(row)
        


        if pop_num in sample_dict.keys():
            # if sample not in sample_dict[pop_num]:
            #     sample_dict[pop_num].append(sample)
            if sample not in sample_dict[pop_num]:
                if qval >= threshold:
                    sample_dict[pop_num].append(sample)
                else:
                    if sample not in trashed_samples:
                        trashed_samples.append(sample)
        else:
            sample_dict[pop_num] = [sample]
    # for key, value in sample_dict.items():
    #     print(f'{key} {value}')
    output_file = generate_population_output(sample_dict, total_populations, ccgp_project_id, ref_genome, trashed_samples)
    
    return output_file

def main():
    #snakemake stuff
    ccgp_project_id = snakemake.params["project_id"]
    ref_genome = snakemake.params["ref_genome"]
    manual_k = snakemake.params["manual_k_assignment"]
    threshold = snakemake.params["threshold_qval"]
    matrix_file = snakemake.input["matrix"]

    df = pd.read_csv(matrix_file)
    # if type(manual_k) == str and manual_k.upper() == "NA":
    #     total_pops = df.at[1, 'best_k'] # change this to a snakemake param later
    if type(manual_k) == int:
        total_pops = manual_k
    else:
        total_pops = df.at[1, 'best_k']

    df_matrix, total_populations = process_df(df, total_pops)
    output_file = assign_pop(df_matrix, total_populations, ccgp_project_id, ref_genome, threshold)

if __name__ == "__main__":
    main()