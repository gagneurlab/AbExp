# AbExp variant effect prediction pipeline

AbExp is a tool to predict aberrant gene underexpression in 48 human tissue based on DNA sequence variants.
It was trained on aberrant gene expression calls from the GTEx dataset.

This repository contains a bioinformatics software pipeline for calculating **AbExp variant effect predictions**, taking vcf files as input.
The preprint to this method can be found on [BioRxiv](https://doi.org/10.1101/2023.12.04.569414).

## Minimum resource requirements

- disk space: 1.0-1.4TB for cache files (hg19 + hg38)
  - LOFTEE: 25GB
  - VEP v108: 113GB
  - CADD v1.6: 854GB
  - SpliceAI-RocksDB (optional): 349GB
- RAM: 64GB
- GPU supporting CUDA for SpliceAI annotation

## Setup

1) install conda and mamba on your system
2) To download the VEP cache (if not existing yet):
   ```bash
   VEP_CACHE_PATH="<your cache path here>"
   VEP_VERSION=108

   mamba env create -f scripts/veff/vep_env.v108.yaml --name vep_v108
   conda activate vep_v108
   
   bash misc/install_vep_cache/install_cache_for_version.sh $VEP_VERSION $VEP_CACHE_PATH
   
   conda deactivate
   ```
3) To download the CADD cache (if not existing yet):
   ```bash
   CADD_CACHE_PATH="<your cache path here>"

   bash misc/download_CADD_v1.6.sh $CADD_CACHE_PATH
   ```
4) To download LOFTEE data and scripts:
   ```bash
   LOFTEE_DIR="<your path here>"

   bash misc/install_vep_cache/download_loftee.sh $LOFTEE_DIR
   ```
4) configure `system_config.yaml`:
   - specify paths to the VEP cache, CADD cache, LOFTEE data and LOFTEE source code as defined in steps 2-4
   - (optional) Disable downloading the SpliceAI-RocksDB cache for pre-computed SpliceAI annotations
   - (optional) Change file paths of automatically downloaded annotations to shared location
   - (optional) Any options defined in `defaults.yaml` can be overwritten in this file if necessary
2) run `mamba env create -f envs/abexp-veff-py.yaml`
3) activate the created environment: `conda activate abexp-veff-py`

## Usage

1) Edit `config.yaml` to your needs:
   - Specify `vcf_input_dir`. All `.vcf|.vcf.gz|.vcf.bgz|.bcf` files in this folder will be annotated.
   - Specify `vcf_is_normalized: True` if all variants are left-normalized and biallelic (`bcftools norm -cs -m`).
     Otherwise, the pipeline will normalize the variants before annotation.
   - Specify the `output_dir`
   - Specify `fasta_file` and `gtf_file` corresponding to the `human_genome_version`.
     The `gtf_file` needs to contain Ensembl gene and transcript identifiers.
     Therefore, it is highly recommended to use the [Gencode genome annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/).
   
   An example is pre-configured and can be used to test the pipeline.

2) Run `snakemake --use-conda -c all`.
   All rules are annotated with resource requirements s.t. snakemake can submit jobs to HPC clusters or cloud environments.
   It is highly recommended to use snakemake with some batch submission system, e.g. SLURM.
   For further information, please visit the Snakemake documentation.
3) The resulting variant effect predictions will be stored in `<output_dir>/predict/abexp_v1.0/<input_vcf_file>.parquet`. It should contain the following columns:
   - 'chrom'
   - 'start'
   - 'end'
   - 'ref'
   - 'alt'
   - 'gene'
   - 'tissue': GTEx tissue, e.g. "Artery - Tibial"
   - 'tissue_type', GTEx tissue type, e.g. "Blood Vessel"
   - 'abexp_v1.0': The predicted AbExp score
   - a set of features used to predict the AbExp score

## Development setup
Advanced users who want to edit this pipeline can use the following steps to convert the python scripts back to Jupyter notebooks:
1) Make sure that the `jupytext` command is available, e.g. via `mamba install jupytext`
2) run `find scripts/ -iname "*[.py.py|.R.R]" -exec jupytext --sync {} \;` to convert all percent scripts to jupyter notebooks
Jupyter will then automatically synchronize the percent scripts with the corresponding notebook files.

