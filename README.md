# AbExp variant effect prediction pipeline

AbExp is a tool to predict aberrant gene expression in 49 human tissue based on DNA sequence variants.
It was trained on aberrant gene expression calls from the GTEx dataset.

This repository contains a bioinformatics software pipeline for calculating **AbExp variant effect predictions**, taking vcf files as input.
The publication to this method can be found in [Nature Communications](https://www.nature.com/articles/s41467-025-58210-w).

## Minimum resource requirements

- Linux
- Disk space: 1.0-1.4TB for cache files (hg19 + hg38)
  - LOFTEE: 25GB
  - VEP v108: 113GB
  - CADD v1.6: 854GB
  - SpliceAI-RocksDB (optional): 349GB
- RAM: 64GB
- GPU supporting CUDA for SpliceAI annotation

## Setup

1) Install conda and mamba on your system.

   Tip: Use the [conda-libmamba-solver](https://conda.github.io/conda-libmamba-solver/user-guide/) for an improved experience when installing environments with `conda`
2) Download the VEP cache (if not existing yet):
   ```bash
   VEP_CACHE_PATH="<your cache path here>"
   VEP_VERSION=108

   mamba env create -f scripts/veff/vep_env.v108.yaml --name vep_v108
   conda activate vep_v108
   
   bash misc/install_vep_cache/install_cache_for_version.sh $VEP_VERSION $VEP_CACHE_PATH
   
   conda deactivate
   ```
3) Download the CADD cache (if not existing yet):
   ```bash
   CADD_CACHE_PATH="<your cache path here>"

   bash misc/download_CADD_v1.6.sh $CADD_CACHE_PATH
   ```
4) Download LOFTEE data and scripts:
   ```bash
   LOFTEE_DIR="<your path here>"

   bash misc/install_vep_cache/download_loftee.sh $LOFTEE_DIR
   ```
4) Configure `system_config.yaml`:
   - specify paths to the VEP cache, CADD cache, LOFTEE data and LOFTEE source code as defined in steps 2-4
   - (optional) Disable downloading the SpliceAI-RocksDB cache for pre-computed SpliceAI annotations by setting absplice.use\_spliceai\_rocksdb to False
   - (optional) Change file paths of automatically downloaded annotations to shared location
   - (optional) Any options defined in `defaults.yaml` can be overwritten in this file if necessary
2) Run `mamba env create -f envs/abexp-veff-py.yaml`
3) Activate the created environment: `conda activate abexp-veff-py`

## Usage

1) Edit the `config.yaml` and specify the following parameters:
   - `vcf_input_dir`. All `.vcf|.vcf.gz|.vcf.bgz|.bcf` files in this folder will be annotated. Genotypes are not required.
   - `vcf_is_normalized: True` if all variants are left-normalized and biallelic (`bcftools norm -cs -m`).
     Otherwise, the pipeline will normalize the variants before annotation.
   - `output_dir`
   - `fasta_file` and `gtf_file` corresponding to the `human_genome_version`.
     The `gtf_file` needs to contain Ensembl gene and transcript identifiers.
     Therefore, it is highly recommended to use the [Gencode genome annotations](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/).

   An example is pre-configured and can be used to test the pipeline.

2) Run `snakemake --use-conda -c all --configfile=config.yaml`.
   All rules are annotated with resource requirements s.t. snakemake can submit jobs to HPC clusters or cloud environments.
   It is highly recommended to use snakemake with some batch submission system, e.g. SLURM.
   For further information, please visit the [Snakemake documentation](https://snakemake.readthedocs.io/).
3) The resulting variant effect predictions will be stored in `<output_dir>/predict/abexp_v1.1/<input_vcf_file>.parquet`. It will contain the following columns:
   - 'chrom': chromosome of the variant
   - 'start': start position of the variant (0-based)
   - 'end': end position of the variant (1-based)
   - 'ref': reference allele
   - 'alt': alternate allele
   - 'gene': the gene affected by the variant
   - 'tissue': GTEx tissue, e.g. "Artery - Tibial"
   - 'tissue_type', GTEx tissue type, e.g. "Blood Vessel"
   - 'abexp_v1.1': The predicted AbExp score
   - a set of features used to predict the AbExp score

## License
All source code and model weights in this repository are licensed under the [MIT license](./LICENSE).

**Please note:** AbExp relies on [CADD](https://cadd.gs.washington.edu/) and [SpliceAI](https://github.com/Illumina/SpliceAI/), both of which are free to use only in non-commercial settings.
If you plan to use AbExp in a commercial context, please ensure that you have the appropriate permissions or licenses to use both tools.

## Development setup
Advanced users who want to edit this pipeline can use the following steps to convert the python scripts back to Jupyter notebooks:
1) Make sure that the `jupytext` command is available, e.g. via `mamba install jupytext`
2) run `find scripts/ -iname "*[.py.py|.R.R]" -exec jupytext --sync {} \;` to convert all percent scripts to jupyter notebooks
Jupyter will then automatically synchronize the percent scripts with the corresponding notebook files.

