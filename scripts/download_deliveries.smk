import pdfkit
from pathlib import Path

rule all:
    input:
        expand("projects/{project_id}/results/{refGenome}/done_downloading.txt", project_id=config['pid'], refGenome=config["ref"]),
        expand("projects/{project_id}/results/{refGenome}/CCGP/CCGP_downloading.txt", project_id=config['pid'], refGenome=config["ref"]),
        expand("projects/{project_id}/results/{refGenome}/trackhub/trackhub_downloading.txt", project_id=config['pid'], refGenome=config["ref"]),
        expand("projects/{project_id}/config_edited.txt", project_id=config['pid'])
        
rule download_from_google:
    output:
        touch("projects/{project_id}/results/{refGenome}/done_downloading.txt")
    shell:
        """
        mkdir -p ./projects/{config[pid]}/results/{config[ref]}/
        
        gsutil -m cp -r \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_callable_sites.bed" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_raw.vcf.gz" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_raw.vcf.gz.tbi" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_clean_indels.vcf.gz" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_clean_indels.vcf.gz.tbi" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_clean_snps.vcf.gz" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_clean_snps.vcf.gz.tbi" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_filtered.vcf.gz" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/{config[pid]}_filtered.vcf.gz.csi" \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/QC" \
            "./projects/{config[pid]}/results/{config[ref]}/"        
        
        gsutil -m cp -r \
            "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/data" \
            "./projects/{config[pid]}/results/{config[ref]}/"
        """

rule ccgp_from_google:
    output:
        touch("projects/{project_id}/results/{refGenome}/CCGP/CCGP_downloading.txt")
    shell:
        """
        gsutil -m rsync -r "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/CCGP" \
            "./projects/{config[pid]}/results/{config[ref]}/CCGP"
        """
rule trackhub_from_google:
    output:
        touch("projects/{project_id}/results/{refGenome}/trackhub_downloading.txt")
    shell:
        """
        mkdir -p ./projects/{config[pid]}/results/{config[ref]}/trackhub
        gsutil -m rsync -r "gs://ccgp-workflow-results/{config[pid]}/results/{config[ref]}/trackhub" \
            "./projects/{config[pid]}/results/{config[ref]}/trackhub"
        """

rule copy_config_yaml:
    input:
        "config/config.yaml"
    output:
        "projects/{project_id}/config.yaml"
    shell:
        """
        cp config/config.yaml projects/{config[pid]}/config.yaml
        """

rule edit_config_yaml:
    input:
        config = "projects/{project_id}/config.yaml"
    output:
        touch(temp("projects/{project_id}/config_edited.txt"))
    shell:
        """
        sed -i 's/GCA_023055505.1/{config[ref]}/g' {input.config}
        sed -i 's/cali_small/{config[pid]}/g' {input.config}
        """