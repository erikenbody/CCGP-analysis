# CCGP-analysis

This repo contains code associated with re-running certain steps of snpArcher (for example, if we decide to remove an additional sample) and running downstream analyses. In brief this means:

This repo contains the following snakemake modules from [snpArcher](https://github.com/harvardinformatics/snpArcher):
* postprocess
* ccgp
* trackhubs

On top of this, some code has been added to the ccgp module. The distinction between ccgp and [algatr](https://github.com/TheWangLab/algatr/) at the moment is nebulous and warrants revisiting sometime later, probably they can be made the same module or possibly algatr modded down to be more simplified.

New modules here include:
* algatr
* local_pca
* postprocess_variant

# How to

1) Install mamba and snakemake 
2) Clone this repo
3) Run the following to test out a small subset of data (within the repo directory):
```
snakemake --dir test_data/ --use-conda -c 10 --profile profiles/silverbullet_reruns/
```

# Notes on implementation

If you want to run on a different project, you must match the file naming and structure of the `test_data` directory. This contains a handful of outputs from the snpArcher workflow:

1) Raw vcf
2) Callable sites bed file
3) List of samples to remove
4) The fai index

A config file is required named `config.yaml` which defines the parameters for the run. This is similar to the snpArcher config, but has a few custom additions. 

# Notes on filters

* clean SNPs has MAF filter set of < 0.01
* setting MAF filter to < 0.05 for the mil sub sampled SNP file. But this hasnt been re run yet
* All VCFs should be filtered for F_MISSING < 0.25 (must be 75% of individuals with data)
