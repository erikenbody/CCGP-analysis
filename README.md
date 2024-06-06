# CCGP-analysis

This repo contains code associated with re-running certain steps of snpArcher (for example, if we decide to remove an additional sample) and running downstream analyses. In brief this means:

This repo contains the following snakemake modules from snpArcher:
* postprocess
* ccgp
* trackhubs

On top of this, some code has been added to ccgp. The distinction between ccgp and algatr at the moment is nebulus and warrants revisiting sometime later, probably theyc an be made the same module or possibly algatr modded down to very simple just algatr and tesss.

New modules here include:
* algatr
* local_pca
* postprocess_variant

# How to

1) Install mamba and snakemake 
2) clone this repo
3) run the following to test out a small subset of data (within the repo directory):
```
snakemake --dir test_data/ --use-conda -c 10 --profile profiles/silverbullet_reruns/
```

# Notes on implementation

The way this directory 


If you want to run on a different project, you must match the file naming and structure of the test_data directory. This contains a handful of outputs from the snparcher workflow:

1) raw vcf
2) callable sites bed file
3) list of samples to remove
4) the fai index

A config file is required `config.yaml` which defines the parameters for the run. This is similar to the snpArcher config, but has a few custom additions. 