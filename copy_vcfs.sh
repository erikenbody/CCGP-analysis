#!/usr/bin/env bash
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/17-Semicossyphus/17-Semicossyphus_final.vcf.gz 17-Semicossyphus_raw.vcf.gz
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/17-Semicossyphus/17-Semicossyphus_final.vcf.gz.tbi 17-Semicossyphus_raw.vcf.gz.tbi
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/26-Cebidichthys/26-Cebidichthys_final.vcf.gz 26-Cebidichthys_raw.vcf.gz
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/26-Cebidichthys/26-Cebidichthys_final.vcf.gz.tbi 26-Cebidichthys_raw.vcf.gz.tbi
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/41-Callipepla/41-Callipepla_final.vcf.gz 41-Callipepla_raw.vcf.gz
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/41-Callipepla/41-Callipepla_final.vcf.gz.tbi 41-Callipepla_raw.vcf.gz.tbi
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/41-Melospiza/41-Melospiza_final.vcf.gz 41-Melospiza_raw.vcf.gz
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/41-Melospiza/41-Melospiza_final.vcf.gz.tbi 41-Melospiza_raw.vcf.gz.tbi
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/eenbody/41-Cyanocitta/41-Cyanocitta_final.vcf.gz 41-Cyanocitta_raw.vcf.gz
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/eenbody/41-Cyanocitta/41-Cyanocitta_final.vcf.gz.tbi 41-Cyanocitta_raw.vcf.gz.tbi
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/eenbody/41-Passerculus/41-Passerculus_final.vcf.gz 41-Passerculus_raw.vcf.gz
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/eenbody/41-Passerculus/41-Passerculus_final.vcf.gz.tbi 41-Passerculus_raw.vcf.gz.tbi
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/eenbody/8-Embiotoca/8-Embiotoca_final.vcf.gz 8-Embiotoca_raw.vcf.gz
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/eenbody/8-Embiotoca/8-Embiotoca_final.vcf.gz.tbi 8-Embiotoca_raw.vcf.gz.tbi

# mkdir -p projects/77-Linanthus/results/GCA_023055425.1
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/77-Linanthus/77-Linanthus_final.vcf.gz projects/77-Linanthus/results/GCA_023055425.1/77-Linanthus_raw.vcf.gz
# scp eenbody@plaza.gi.ucsc.edu:/public/groups/corbettlab/mara/77-Linanthus/77-Linanthus_final.vcf.gz.tbi projects/77-Linanthus/results/GCA_023055425.1/77-Linanthus_raw.vcf.gz.tbi

# gdrive files download 1QivaSUb5fE5H8EcDMEr2ABOWCLKbTB5h --destination projects/73-Haliotis/results/GCA_023055435.1/
# gdrive files download 1bC8nFnW3wfWSD0Y0RdvF9dYuYridr2Gn --destination projects/7-Clinocottus/results/GCA_023055335.1/
# gdrive files download 1nwjbpNZBJg477imbu_kqWb3r16XbYFFo --destination projects/7-Clinocottus/results/GCA_023055335.1/

#get QC directories

project="41-Callipepla"

for project in "41-Cyanocitta" "41-Melospiza" "41-Passerculus" "41-Perognathus" "77-Linanthus"
do
    echo $project
    rsync -av --exclude='*filtered*.vcf.gz'  hgdownload.soe.ucsc.edu::ccgp/${project}/QC/ ./projects/${project}/results/
done

gdrive files download 126mg2geM4BhlOdh8FEMG4vMLQ0bneDQc --destination Rmen --recursive
gdrive files download 1n-nrQMxj0qSf7I6lI4ll5b3502uTeJmM --destination . --recursive