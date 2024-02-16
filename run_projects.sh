#!/bin/bash

export MAMBA_LOCK_TIMEOUT=600  # Set timeout to 10 minutes

# Assuming you are in the parent directory containing "projects"
parent_directory="/scratch2/erik/CCGP-reruns"

# Define the function to run Snakemake for a given directory
run_snakemake() {
    directory="$1"
    dir_name=$(basename "$directory")
    snakemake --dir "$directory" --config config_file="config.yaml" -c 10 --use-conda --rerun-triggers mtime
    # Add any additional logic or flags as needed
}

# Export the function to make it available to GNU Parallel
export -f run_snakemake

# Use GNU Parallel to run the function in parallel for each directory
find "$parent_directory"/projects/*/ -maxdepth 0 -type d | parallel -j 10 run_snakemake

#this means 10 jobs running at one time, i.e. 10 snakemake runs
