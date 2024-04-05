import sys
sys.path.append(
    ".."
)
import pandas as pd
import os
from db import get_mongo_client
#these should be native
#import argparse
#mport re

client = get_mongo_client()
db = client["ccgp_dev"]
collection = db["sample_metadata"]
wd_scripts = os.getcwd()
print(wd_scripts)
os.chdir("..")
wd_parent = os.getcwd() #this should be ccgp-reruns/

def mongo_get_coords(project, ref_genome, sample_names, vcf_samples):

    data_list = []
    query_criteria = {"ccgp-project-id": project}

    sample_entries_for_project = collection.find(query_criteria)
    #print(len(vcf_samples))
    for sample in sample_entries_for_project:
        if sample["*sample_name"] not in vcf_samples:
            #print(f"Sample {sample["*sample_name"]} not in VCF... Skipping.")
            continue

        sample_name = sample.get("*sample_name", "")
        lat = sample.get("lat", "")
        long = sample.get("long", "")

        data_list.append({"Sample Name": sample_name, "Long": long, "Lat": lat})
        print(sample_name)
    
    df = pd.DataFrame(data_list)

    num_rows = df.shape[0]
    if num_rows == 0:
        print(f"No results found for project {project}. Skipping...")


    print(f"Number of samples in the {project}: {num_rows}")

    output_file_path = os.path.join(wd_scripts, "results", ref_genome, "algatr", f"{project}.coords.txt")
    #print(output_file_path)

    df.to_csv(output_file_path, sep="\t", index=False, header=False)
    print(f'Got {project} coords.')
    print('')

if __name__ == "__main__":
    ccgp_project_id = snakemake.params["project_id"]
    ref_genome = snakemake.params["ref_genome"]
    samples = snakemake.params["sample_id"]

    vcf_samples = []
    coord_samples_file = f"/scratch2/erik/CCGP-reruns/projects/{ccgp_project_id}/results/{ref_genome}/algatr/{ccgp_project_id}.samps4coords.txt"
    with open(coord_samples_file) as file:
        for sample in file:
            vcf_samples.append(sample.strip())

    mongo_get_coords(ccgp_project_id, ref_genome, samples, vcf_samples)

    

