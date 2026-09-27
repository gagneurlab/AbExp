# AbExp variant effect prediction pipeline

AbExp is a tool to predict aberrant gene expression in 49 human tissue based on DNA sequence variants.
It was trained on aberrant gene expression calls from the GTEx dataset.

This repository contains a bioinformatics software pipeline for calculating **AbExp variant effect predictions**, taking vcf files as input.
The publication to this method can be found in [Nature Communications](https://www.nature.com/articles/s41467-025-58210-w). We also offer a [web interface](https://abexp.cmm.cit.tum.de/) for querying AbExp scores on any SNP.

## Minimum resource requirements

- Linux with tar and gzip
- Disk space: about 130GB for the downloaded resources of one genome assembly (hg38), plus the optional SpliceAI-RocksDB
  - VEP v108 cache: 26GB (hg38), 16GB (hg19); twice as much while the download is extracted
  - CADD v1.6: 88GB (hg38), 84GB (hg19)
  - LOFTEE: 14GB (hg38), 1.3GB (hg19)
  - GTEx and SpliceMap tables: 2GB
  - SpliceAI-RocksDB (optional): 349GB for hg19 + hg38
- RAM: 64GB
- (optional) a GPU with CUDA for Enformer and SpliceAI, see [GPU](#gpu)

## Setup

1) Install conda and mamba on your system, and create the environment with Snakemake 9 and aria2c:
   ```bash
   mamba env create -f workflow/envs/abexp-veff-py.yaml
   conda activate abexp-veff-py
   ```
   Tip: Use the [conda-libmamba-solver](https://conda.github.io/conda-libmamba-solver/user-guide/) for an improved experience when installing environments with `conda`
2) Configure `config/config.yaml`, see [Usage](#usage).
   The pipeline downloads the resources it needs, e.g. the VEP cache, CADD and LOFTEE, for the configured
   `human_genome_version` only: more than 100 GB. By default, they go to `<output_dir>/resources`, so each output folder gets its own copy.
   Set `dirs.resources_dir` in the `system` section to a shared folder, so that all output folders use one copy. Other options in the `system` section:
   - `vep.vep_cache_dir`, `vep.cadd_dir`, `vep.loftee_data_dir` and `vep.loftee_src_path`: existing caches to reuse instead of downloading
   - `absplice.use_spliceai_rocksdb: False` to skip the SpliceAI-RocksDB download
   - Any option in `workflow/schemas/config.schema.yaml` can be set in this section.
     The schema also lists the defaults. The options of `vep`, `mehari`, `absplice` and `enformer` are in
     `workflow/modules/veff/<module>/config.schema.yaml`.
3) (optional) Download the resources before the first run: `snakemake setup -c 4`.
   Otherwise, the first run downloads them.
   With a cluster executor, the downloads run as cluster jobs. If the compute nodes have no internet access,
   run `snakemake setup -c 4` without the executor on a host with internet access.
   An interrupted download of the VEP cache, CADD, LOFTEE or SpliceAI-RocksDB resumes on the next run, and
   Snakemake write-protects these files.
   If a later version of AbExp changes a download rule, Snakemake stops with `ProtectedOutputException`.
   Run `snakemake --cleanup-metadata <files>` to keep the downloaded files.
   Earlier versions downloaded SpliceAI-RocksDB in a conda environment. Its download rule has changed, and the
   databases are not write-protected, so Snakemake would delete them and download them again. To keep them, run
   `snakemake --cleanup-metadata <resources_dir>/spliceai_rocksdb/spliceAI_hg38_chr*.db` once (for hg19:
   `spliceAI_hg19_chr*.db`; or the paths of `absplice.spliceai_rocksdb_path`). The `*.db_backup.tar.gz` files
   and `*.db_backup` folders next to the databases are not needed; deleting them frees about 300GB per genome
   assembly.
4) (optional) Create the conda environments before the first run: `snakemake --conda-create-envs-only`.

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
   - (optional) `use_gpu: True` to run Enformer and SpliceAI on a GPU, see [GPU](#gpu).

   An example is pre-configured and can be used to test the pipeline.

2) Run `snakemake -c all`. Snakemake reads `config/config.yaml` by default;
   use `--configfile my_config.yaml` for another config file.
   The workflow profile `workflow/profiles/default` sets `--sdm conda`; `--sdm` on the command line overrides it.
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

## GPU

`use_gpu: True` in the config makes Enformer and SpliceAI use the CUDA variant of their TensorFlow
environment (`workflow/modules/veff/envs/abexp-tensorflow-cuda.yaml`) instead of the CPU variant. The CUDA
variant also runs on hosts without a GPU, on the CPU. Only with `use_gpu: True`, these rules request one GPU
(resource `gpu`), e.g. from SLURM.

Conda creates the CUDA variant only on a host with a CUDA driver. To create it on a host without one, e.g. a
login node, set `CONDA_OVERRIDE_CUDA` to a CUDA version that the driver of the GPU nodes supports, e.g. 12.9:
```bash
CONDA_OVERRIDE_CUDA=12.9 snakemake --conda-create-envs-only
```

## Using mehari instead of VEP

[mehari](https://github.com/varfish-org/mehari) can replace VEP for the transcript consequence annotation.
It needs no VEP cache, CADD or LOFTEE, so the pipeline does not download them, and it runs much faster.
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
3) Run snakemake. The mehari rules use the environment `workflow/modules/veff/mehari/envs/mehari_env.yaml`,
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
workflow/modules/veff/envs/                 # conda environments that several veff modules use
```

`workflow/modules/veff/loftee` runs [reloftee](https://github.com/gagneurlab/reloftee), a
VEP-free reimplementation of LOFTEE. It calls its own loss-of-function consequences from
`gff3_file`, or from another Ensembl or GENCODE annotation in `system.loftee.genome_annotation`.
It is off by default; set `system.loftee.enabled: true` and `system.loftee.conda_env` to run it.
See its config.schema.yaml for the required and optional inputs.

`workflow/modules/veff/nmd_scanner` adds [NMD-Scanner](https://github.com/gagneurlab/NMD-Scanner),
which scans variants for premature termination codons and evaluates the NMD escape rules on the
transcripts of a GTF file, and keeps its own per-transcript, per-variant table. It is off by
default; set `system.nmd_scanner.enabled: true` to run it. NMD-Scanner 0.2.0 reads only GTF, so
also set `system.nmd_scanner.gtf_file` to the GENCODE GTF of the same release as `gff3_file`.

Both modules keep the output of their tool as it is, and neither is read by tissue_specific_vep,
fset or predict.

Each module has its own `config.schema.yaml` with its options and defaults, its own scripts and
conda environments (`envs/`), and the files it ships. The conda environments that several veff modules
use are in `workflow/modules/veff/envs/`. Every rule that needs more than a shell has a
`conda:` environment. The download rules have no conda environment and use aria2c, tar and gzip from the
PATH. So the importing workflow needs Snakemake 9 and aria2c (conda-forge package `aria2`) in one
environment, and `--sdm conda`.

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
    # downloads of the VEP cache, CADD and LOFTEE; default is <output_dir>/resources
    "resources_dir": "resources",
    # optional: existing caches to reuse instead of downloading
    # "vep_cache_dir": "<VEP cache>/{vep_version}",
    # "cadd_dir": "<CADD v1.6>/{human_genome_assembly}",
    # "loftee_data_dir": "<LOFTEE data>/{human_genome_assembly}",
    # "loftee_src_path": "<LOFTEE source>/{human_genome_assembly}_src",
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
- Modules with downloads have a target rule that lists them, e.g. `veff__vep_setup` (`vep_veff__vep_setup`
  in the example). Run it to download the resources before the first run.
- The modules vep and loftee download the same LOFTEE data. To use both, point the loftee module at the
  files of the vep module and import it without its download rules, as `workflow/modules/veff/Snakefile` does.

## License
All source code and model weights in this repository are licensed under the [MIT license](./LICENSE).

**Please note:** AbExp relies on [CADD](https://cadd.gs.washington.edu/) and [SpliceAI](https://github.com/Illumina/SpliceAI/), both of which are free to use only in non-commercial settings.
If you plan to use AbExp in a commercial context, please ensure that you have the appropriate permissions or licenses to use both tools.

## Development setup
Advanced users who want to edit this pipeline can use the following steps to convert the python scripts back to Jupyter notebooks:
1) Make sure that the `jupytext` command is available, e.g. via `mamba install jupytext`
2) run `find workflow/ -iname "*[.py.py|.R.R]" -exec jupytext --sync {} \;` to convert all percent scripts to jupyter notebooks
Jupyter will then automatically synchronize the percent scripts with the corresponding notebook files.

