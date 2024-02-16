#!/bin/bash

echo "" > run_projects.sh

# Loop through project directories
for project_dir in projects/*/; do
    # Get project and genome names from the directory structure
    PROJECTNAME=$(basename "$project_dir")
    GENOMENAME=$(basename "$(find "$project_dir/results" -mindepth 1 -maxdepth 1 -type d )")

    # Create the new project config file
    new_config_file="projects/$PROJECTNAME/config.yaml"
    cp config/config.yaml "$new_config_file"

    # Update the final_prefix and refgenome values in the new config file
    sed -i "s/final_prefix: .*/final_prefix: $PROJECTNAME/" "$new_config_file"
    sed -i "s/refgenome: .*/refgenome: $GENOMENAME/" "$new_config_file"

    bed="projects/${PROJECTNAME}/results/${GENOMENAME}/${PROJECTNAME}_callable_sites.bed"
    if [[ -e "$bed" ]]; then
        sed -i "s/bed: .*/bed: true/" "$new_config_file"
    fi

    #create empty file if no samps file yet

    samps_file="projects/${PROJECTNAME}/results/${GENOMENAME}/${PROJECTNAME}_samps.txt"

    if [[ ! -e "$samps_file" ]]; then
        echo "" > "$samps_file"
    fi

    echo "snakemake -c 20 --use-conda --config config_file="config.yaml" --dir "projects/${PROJECTNAME}" " >> run_projects.sh
done
