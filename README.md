# AbExp variant effect prediction pipeline

AbExp is a tool to predict aberrant gene expression in 49 human tissue based on DNA sequence variants.
It was trained on aberrant gene expression calls from the GTEx dataset.

This repository contains a bioinformatics software pipeline for calculating **AbExp variant effect predictions**, taking vcf files as input.
The publication to this method can be found in [Nature Communications](https://www.nature.com/articles/s41467-025-58210-w). We also offer a [web interface](https://abexp.cmm.cit.tum.de/) for querying AbExp scores on any SNP.

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

   mamba env create -f workflow/modules/veff/vep/envs/vep_env.v108.yaml --name vep_v108
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
5) Configure the `system` section of `config/config.yaml`:
   - specify paths to the VEP cache, CADD cache, LOFTEE data and LOFTEE source code as defined in steps 2-4
   - (optional) Disable downloading the SpliceAI-RocksDB cache for pre-computed SpliceAI annotations by setting absplice.use\_spliceai\_rocksdb to False
   - (optional) Change file paths of automatically downloaded annotations to shared location
   - (optional) Any option in `workflow/schemas/config.schema.yaml` can be set in this section.
     The schema also lists the defaults. The options of `vep`, `mehari`, `absplice` and `enformer` are in
     `workflow/modules/veff/<module>/config.schema.yaml`.
6) Run `mamba env create -f workflow/envs/abexp-veff-py.yaml`. The environment contains Snakemake 9.
7) Activate the created environment: `conda activate abexp-veff-py`

## Usage

1) Edit the `config/config.yaml` and specify the following parameters:
   - `vcf_input_dir`. All `.vcf|.vcf.gz|.vcf.bgz|.bcf` files in this folder will be annotated. Genotypes are not required.
   - `vcf_is_normalized: True` if all variants are left-normalized and biallelic (`bcftools norm -cs -m`).
     Otherwise, the pipeline will normalize the variants before annotation.
   - `output_dir`
   - `fasta_file` and `gff3_file` corresponding to the `human_genome_version`.
     `gff3_file` is a [GENCODE genome annotation](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/) in GFF3 format,
     e.g. `gencode.v42.annotation.gff3.gz`.

   - (optional) `veff.annotator: "mehari"` to annotate transcript consequences with
     [mehari](https://github.com/varfish-org/mehari) instead of VEP, together with
     `veff.mehari_gencode_transcripts_fasta`, see below.

   An example is pre-configured and can be used to test the pipeline.
   A config passed with `--configfile` replaces `config/config.yaml`, see step 2. A key that it leaves out
   takes the default of `workflow/schemas/config.schema.yaml`, not the value of `config/config.yaml`.

2) Run `snakemake --sdm conda -c all`. Snakemake reads `config/config.yaml` by default;
   use `--configfile my_config.yaml` for another config file.
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

## Using mehari instead of VEP

[mehari](https://github.com/varfish-org/mehari) can replace VEP for the transcript consequence annotation.
It needs no VEP cache and no LOFTEE, and it runs much faster.
The module `workflow/modules/veff/mehari` calls the mehari Python package and keeps its output.
The per-transcript table has mehari's own consequence terms, and there are no LoF or NMD calls and no
CADD, SIFT, PolyPhen or Condel scores.
The shipped AbExp models were trained on VEP features and do not run on this route.
It is meant for training new models.

Setup:
1) Download the GENCODE transcript FASTA of the release of your `gff3_file`, e.g.
   `bash misc/mehari/download_gencode.sh 42 GRCh38 data/gencode/release_42` (GENCODE 42 = Ensembl 108).
   The script also fetches the GFF3 of that release.
2) In `config/config.yaml`, set `veff.annotator: "mehari"` and `veff.mehari_gencode_transcripts_fasta`.
   The pipeline builds the mehari transcript database from `gff3_file` and the FASTA (about 6 minutes and 4 GB RAM for
   a full GRCh38 release) and stores it at `system.mehari.transcripts_db` (see `workflow/modules/veff/mehari/config.schema.yaml`).
3) Run snakemake with `--sdm conda`. The mehari rules use the environment `workflow/modules/veff/mehari/envs/mehari_env.yaml`,
   which installs the mehari Python package from bioconda.
   If you already have a conda environment with the mehari Python package, set `system.mehari.conda_env`
   in the config to its name instead.

## Using AbExp as Snakemake modules

Other workflows can import all of AbExp, only the variant annotation, or single steps like VEP,
with the Snakemake [`module`](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html) directive.

```
workflow/Snakefile                          # AbExp: vcf_prep, gtf_transcripts, veff, feature sets, prediction
workflow/modules/vcf_prep/                  # normalizes and strips the VCFs, sets the variant IDs
workflow/modules/gtf_transcripts/           # transcripts of the GFF3 file as parquet
workflow/modules/veff/Snakefile             # variant annotation: prepared VCFs in, per-variant features out
workflow/modules/veff/vep/                  # VEP with LOFTEE and CADD
workflow/modules/veff/mehari/               # mehari
workflow/modules/veff/loftee/               # reloftee: LOFTEE loss-of-function calls
workflow/modules/veff/tissue_specific_vep/  # consequences per GTEx tissue
workflow/modules/veff/absplice/             # AbSplice-DNA
workflow/modules/veff/enformer/             # Enformer
workflow/modules/veff/nmd_scanner/          # NMD-Scanner
```

`workflow/modules/veff/loftee` runs reloftee, a VEP-free reimplementation of LOFTEE that is not
published yet. It calls its own loss-of-function consequences from
`gff3_file`, or from another Ensembl or GENCODE annotation in `system.loftee.genome_annotation`.
The module is off by default. Set `system.loftee.enabled: true` and `system.loftee.conda_env`, and request
`<output_dir>/veff/loftee/veff.parquet/<vcf_file>.parquet` as target.
See its config.schema.yaml for the required and optional inputs.

`workflow/modules/veff/nmd_scanner` adds [NMD-Scanner](https://github.com/gagneurlab/NMD-Scanner),
which scans variants for premature termination codons and evaluates the NMD escape rules on the
transcripts of a GTF file, and keeps its own per-transcript, per-variant table. It is off by
default; set `system.nmd_scanner.enabled: true` to run it. NMD-Scanner 0.2.0 reads only GTF, so
also set `system.nmd_scanner.gtf_file` to the GENCODE GTF of the same release as `gff3_file`.

Both modules keep the output of their tool as it is, and neither is read by tissue_specific_vep,
fset or predict.

Each module has its own `config.schema.yaml` with its options and defaults, its own scripts and
conda environments (`envs/`), and the files it ships. Every rule that needs more than a shell has a
`conda:` environment, so the importing workflow needs only Snakemake 9 and `--sdm conda`.

Example: VEP only. The VEP module reads the variant IDs that `vcf_prep` sets, so it takes the
stripped VCFs of `vcf_prep` as input:

```python
ABEXP = "gagneurlab/AbExp"
TAG = "<release tag>"
# mapping of chromosome names, e.g. a copy of workflow/modules/vcf_prep/resources/chromAlias.tsv
CHROM_ALIAS = "chromAlias.tsv"

VCF_PREP_CONFIG = {
    "vcf": "input/{vcf_file}",
    "output_dir": "results/abexp",
    "fasta_file": "genome.fa",
    "chrom_alias_tsv": CHROM_ALIAS,
}

module vcf_prep:
    snakefile: github(ABEXP, path="workflow/modules/vcf_prep/Snakefile", tag=TAG)
    config: VCF_PREP_CONFIG

use rule * from vcf_prep as vcf_prep_*

VEP_CONFIG = {
    "vcf": str(rules.vcf_prep_extract_vcf_variants.output.vcf_file),
    "output_dir": "results/abexp/veff",
    "human_genome_version": "hg38",
    "fasta_file": "genome.fa",
    "chrom_alias_tsv": CHROM_ALIAS,
    "vep_cache_dir": "<VEP cache>/{vep_version}",
    "cadd_dir": "<CADD v1.6>/{human_genome_assembly}",
    "loftee_data_dir": "<LOFTEE data>/{human_genome_assembly}",
    "loftee_src_path": "<LOFTEE source>/{human_genome_assembly}_src",
}

module vep:
    snakefile: github(ABEXP, path="workflow/modules/veff/vep/Snakefile", tag=TAG)
    config: VEP_CONFIG

use rule * from vep as vep_*
```

The VEP table is then at `results/abexp/veff/vep/veff.parquet/<vcf_file>.parquet`.

- Build the config dicts before the `module` statement. Snakemake does not accept a dict that
  spans several lines after `config:`.
- Connect modules with `rules.<rule>.output`, as in the example. `workflow/modules/veff/Snakefile`
  does the same for all steps.
- For several configurations, e.g. hg19 and hg38, import a module once per configuration, each with
  its own output directory. The instances must not share `resources_dir`, or their download rules
  produce the same files.

## License
All source code and model weights in this repository are licensed under the [MIT license](./LICENSE).

**Please note:** AbExp relies on [CADD](https://cadd.gs.washington.edu/) and [SpliceAI](https://github.com/Illumina/SpliceAI/), both of which are free to use only in non-commercial settings.
If you plan to use AbExp in a commercial context, please ensure that you have the appropriate permissions or licenses to use both tools.

## Development setup
Advanced users who want to edit this pipeline can use the following steps to convert the python scripts back to Jupyter notebooks:
1) Make sure that the `jupytext` command is available, e.g. via `mamba install jupytext`
2) run `find workflow/ -iname "*[.py.py|.R.R]" -exec jupytext --sync {} \;` to convert all percent scripts to jupyter notebooks
Jupyter will then automatically synchronize the percent scripts with the corresponding notebook files.

