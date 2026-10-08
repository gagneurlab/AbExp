# AbExp variant effect prediction pipeline

AbExp is a tool to predict aberrant gene expression in 49 human tissue based on DNA sequence variants.
It was trained on aberrant gene expression calls from the GTEx dataset.

This repository contains a bioinformatics software pipeline for calculating **AbExp variant effect predictions**, taking vcf files as input.
The publication to this method can be found in [Nature Communications](https://www.nature.com/articles/s41467-025-58210-w). We also offer a [web interface](https://abexp.cmm.cit.tum.de/) for querying AbExp scores on any SNP.

## Getting started

1) Create the environment with Snakemake, see step 1 of [Setup](#setup).
2) Run the pre-configured example from the repository root: `snakemake -c 4`. It annotates the ClinVar variants on
   chr22 in `example/clinvar_chr22` (hg38) and predicts with `abexp_v1.1`. The results go to
   `example/output_hg38/predict/abexp_v1.1/`. The first run downloads the resources, see
   [Minimum resource requirements](#minimum-resource-requirements). `snakemake -n` lists the jobs without running them.
3) For your own VCFs, edit `config/config.yaml`, see [Usage](#usage). To run only the annotation without the CADD and
   LOFTEE downloads, see [Turning off CADD or LOFTEE](#turning-off-cadd-or-loftee).

## Minimum resource requirements

- Linux with tar and gzip
- Disk space for the downloaded resources of one genome assembly (hg38): about 300GB with the default config, or
  about 130GB with `absplice.use_spliceai_rocksdb: False` in the `system` section
  - VEP v108 cache: 26GB (hg38), 16GB (hg19); twice as much while the download is extracted
  - CADD v1.6: 88GB (hg38), 84GB (hg19)
  - LOFTEE: 14GB (hg38), 1.3GB (hg19)
  - GTEx and SpliceMap tables: 2GB
  - SpliceAI-RocksDB: about 175GB per genome assembly, 349GB for hg19 + hg38
  - Enformer model, for models with Enformer features: 1GB
- RAM: 64GB
- (optional) a GPU with CUDA for Enformer, SpliceAI and Pangolin, see [GPU](#gpu)

## Setup

1) Install conda and mamba on your system, and create the environment with Snakemake 9, aria2c and kagglehub:
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
   - `vep.cadd_plugin: False` and `vep.loftee_plugin: False` to skip CADD and LOFTEE; the shipped models then do not
     run, see [Turning off CADD or LOFTEE](#turning-off-cadd-or-loftee)
   - `absplice.use_spliceai_rocksdb: False` to skip the SpliceAI-RocksDB download
   - Any option in `workflow/schemas/config.schema.yaml` can be set in this section.
     The schema also lists the defaults. The options of each veff module, e.g. `vep` or `absplice`, are in
     `workflow/modules/veff/<module>/config.schema.yaml`.
3) (optional) Download the resources before the first run: `snakemake setup -c 4`.
   Otherwise, the first run downloads them. For models with Enformer features, the resources include the
   Enformer model in `<resources_dir>/kagglehub`, and the Enformer jobs load it from there.
   With a cluster executor, the downloads run as cluster jobs. If the compute nodes have no internet access,
   run `snakemake setup -c 4` without the executor on a host with internet access.
   If several runs share `resources_dir`, run `setup` once before them. Snakemake locks files only within one
   working directory, so two runs from different directories could download the same file at the same time.
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
   - (optional) `use_gpu: True` to run Enformer, SpliceAI and Pangolin on a GPU, see [GPU](#gpu).

   An example is pre-configured and can be used to test the pipeline.
   A config passed with `--configfile` replaces `config/config.yaml`, see step 2. A key that it leaves out
   takes the default of `workflow/schemas/config.schema.yaml`, not the value of `config/config.yaml`.

2) Run `snakemake -c all`. Snakemake reads `config/config.yaml` by default;
   use `--configfile my_config.yaml` for another config file.
   The workflow profile `workflow/profiles/default` sets `--sdm conda` and the rerun triggers, see
   [When jobs rerun](#when-jobs-rerun). Options on the command line override it.
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

`use_gpu: True` in the config makes Enformer, MMSplice, SpliceAI and AbSplice-DNA use the CUDA variant of their
TensorFlow environment (`workflow/modules/veff/envs/abexp-tensorflow-cuda.yaml`) instead of the CPU variant.
Pangolin of the absplice2 module likewise uses the CUDA variant of its PyTorch environment
(`workflow/modules/veff/absplice2/envs/absplice2-pangolin-cuda.yaml`). The CUDA variants also run on hosts
without a GPU, on the CPU. Only with `use_gpu: True`, the rules of Enformer, MMSplice, SpliceAI and Pangolin request
one GPU (resource `gpu`), e.g. from SLURM. AbSplice-DNA requests none.

Conda creates the CUDA variants only on a host with a CUDA driver. To create them on a host without one, e.g.
a login node, set `CONDA_OVERRIDE_CUDA` to a CUDA version that the driver of the GPU nodes supports, e.g. 12.9:
```bash
CONDA_OVERRIDE_CUDA=12.9 snakemake --conda-create-envs-only
```

## Turning off CADD or LOFTEE

By default, VEP runs the plugins CADD and LOFTEE. Their data take 88GB (CADD) and 14GB (LOFTEE) for hg38.
`system.vep.cadd_plugin: False` or `system.vep.loftee_plugin: False` turns a plugin off. VEP then runs without it,
and neither the run nor `snakemake setup` downloads its data. The LOFTEE source is still downloaded (5MB),
because the VEP plugin MaxEntScan reads it. reloftee (`system.loftee.enabled: True`) still downloads the LOFTEE
data, because it reads them. Changing either option reruns VEP and the steps after it.

The shipped models (abexp_v1.0, abexp_v1.1, abexp_v1.1_Enformer and abexp_v1.1_nobcv) need both plugins: they read
the features `cadd_raw.max` and `LoF_HC.proportion`. So with a plugin off, only the annotation runs, without
prediction:
- Set `predict_abexp_models: []`. Otherwise the workflow stops at the start with an error that names the model and
  the option.
- Request the annotation as target, e.g.
  `snakemake -c 4 <output_dir>/veff/tissue_specific_vep.py/veff.parquet/<input_vcf_file>.parquet`.

Turn a plugin off only if you do not predict with the shipped models, e.g. to train a new model without CADD or
LOFTEE features, or to try the annotation without the large downloads.

## Using mehari instead of VEP

[mehari](https://github.com/varfish-org/mehari) can replace VEP for the transcript consequence annotation.
It needs no VEP cache, CADD or LOFTEE, so the pipeline does not download them, and it runs much faster.
The module `workflow/modules/veff/mehari` calls the mehari Python package and keeps its output.
The per-transcript table has mehari's own consequence terms, and there are no LoF or NMD calls and no
CADD, SIFT, PolyPhen or Condel scores.
The shipped AbExp models were trained on VEP features and do not run on this route.
It is meant for training new models. Set `predict_abexp_models: []`, and request the annotation as target as in
[Turning off CADD or LOFTEE](#turning-off-cadd-or-loftee). Otherwise the workflow stops at the start.

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

## Model format

The AbExp models are LightGBM text models, e.g. `resources/models/abexp_v1.1/model.txt`.
Any LightGBM 4.x reads them with `lightgbm.Booster(model_file=...)`.
The feature names of a model are in the `feature_names=` line of its model.txt.
The predict rule selects the feature columns by these names, in the order of training.
It fails if the feature set lacks one of them.
To add a model, set `model` of an entry in `system.models`, see `workflow/schemas/config.schema.yaml`.

Earlier versions stored the models as joblib pickles (`model.joblib`), which need LightGBM 3.3.
The predict rule no longer reads them. To keep using a custom joblib model, convert it once to a text model.
The conversion needs LightGBM 3.3, scikit-learn, joblib, and `packages/aberrant_expression` for the AbExp model
wrapper. The pickles name the module of the wrapper by its old name, `abexp_utils.models.wrappers` or
`rep.models.wrappers`, so the conversion maps both to `abexp.utils.models.wrappers`.
Run it in the root of this repository:
```bash
mamba create -n abexp-convert -c conda-forge python=3.11 "lightgbm~=3.3" "scikit-learn<1.8" joblib
PYTHONPATH=packages/aberrant_expression/src mamba run -n abexp-convert python -c '
import joblib
import joblib.numpy_pickle

class Unpickler(joblib.numpy_pickle.NumpyUnpickler):
    def find_class(self, module, name):
        if module in ("abexp_utils.models.wrappers", "rep.models.wrappers"):
            module = "abexp.utils.models.wrappers"
        return super().find_class(module, name)

joblib.numpy_pickle.NumpyUnpickler = Unpickler
model = joblib.load("my_model/model.joblib")
# the shipped models wrap an LGBMRegressor; for a bare LGBMRegressor, use model.booster_
model.model.booster_.save_model("my_model/model.txt")
'
```
The `feature_names=` line of the new model.txt must list the column names of the feature set.
With scikit-learn below 1.8, `predict()` of the joblib model also runs in this environment, e.g. to compare it with
the text model. For an `AbExpZscoreRegressor` like the shipped models, `Booster.predict()` of the text model returns
exactly the values of `predict()` of the joblib model.

## NMD efficiency model

The rule `veff__nmd_scanner_score` of the nmd_scanner module predicts `nmd_pred_score` with the NMD efficiency random forest of
[NMD-Scanner](https://github.com/gagneurlab/NMD-Scanner).
The model is `workflow/modules/veff/nmd_scanner/resources/nmd_efficiency_rf.onnx`
(sha256 `e19f3c2b6c55c70b6452f657fea7e06a98ac33295b6128e0de0f2ef8fbd700b8`).
The rule runs it with onnxruntime and needs neither scikit-learn nor a pickle.
The model has one input `input` of type double and shape [N, 19].
It has one output `variable` of type float and shape [N, 1].
The metadata key `feature_names` holds the names of the 19 input columns in input order, as a JSON list.
The rule takes the input names from `nmd_scanner.schema.MODEL_INPUTS`, which `veff__nmd_scanner_annotation`
stores in the parquet metadata of its output, and stops with an error if they differ from `feature_names`.
It passes the boolean columns as 0.0 and 1.0.

The ONNX file is [`nmd_efficiency_rf.onnx`](https://github.com/gagneurlab/NMD-Scanner/blob/v0.5.0/nmd_efficiency_rf.onnx) of NMD-Scanner 0.5.0, byte for byte.
NMD-Scanner's `scripts/train_model.py` trains this random forest with scikit-learn 1.9.1 on the NMD-Scanner 0.4.0
features of the NMDEff TCGA benchmark, and writes it as ONNX.
The script writes the `TreeEnsembleRegressor` node (`ai.onnx.ml` opset 3) itself, so the split thresholds and leaf
values stay in double precision.
The `RandomForestRegressor` converter of skl2onnx would store them as float32.
By the ONNX specification, this node outputs float32.
The script checks that the output equals `predict()` of the random forest, rounded to float32, on all training rows.
On the 14342 scored rows of the ClinVar chr22 example, the output also equaled it in every row.

To redo the model, run in the root of a clone of NMD-Scanner at tag v0.5.0:
```bash
uv run scripts/train_model.py --gff3 gencode.v42.annotation.gff3.gz --fasta GRCh38.fa --out-dir out/
sha256sum out/models/nmd_efficiency_rf.onnx
```
Then copy `out/models/nmd_efficiency_rf.onnx` to `workflow/modules/veff/nmd_scanner/resources/`.
The script pins its dependencies, and two runs gave the same file byte for byte.

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
workflow/modules/veff/absplice2/            # AbSplice2-DNA
workflow/modules/veff/enformer/             # Enformer
workflow/modules/veff/nmd_scanner/          # NMD-Scanner, NMD efficiency, NMD features
workflow/modules/veff/envs/                 # conda environments that several veff modules use
```

`workflow/modules/veff/loftee` runs reloftee, a VEP-free reimplementation of LOFTEE that is not
published yet. It calls its own loss-of-function consequences from
`gff3_file`, or from another Ensembl or GENCODE annotation in `system.loftee.genome_annotation`.
The module is off by default. Set `system.loftee.enabled: true` and `system.loftee.conda_env`, and request
`<output_dir>/veff/loftee/veff.parquet/<vcf_file>.parquet` as target.
See its config.schema.yaml for the required and optional inputs.

`workflow/modules/veff/nmd_scanner` adds [NMD-Scanner](https://github.com/gagneurlab/NMD-Scanner),
which scans variants for premature termination codons (PTCs) and evaluates the NMD escape rules on the
transcripts of `gff3_file`. It keeps NMD-Scanner's own per-transcript, per-variant table without the
sequence columns. Two more rules build the NMD features from this table:
- `veff__nmd_scanner_score` predicts the NMD efficiency `nmd_pred_score` of each transcript with a PTC,
  with the random forest of NMD-Scanner (see [NMD efficiency model](#nmd-efficiency-model)). It scores the
  rows whose `nmd_model_status` is "ok", and its log gives the number of rows per status.
- `veff__nmd_scanner_features` aggregates the scores per variant, gene and GTEx tissue. The GTEx isoform
  proportions of the tissue_specific_vep module weight them: the weight of a transcript is its median
  transcript proportion in the tissue, divided by the sum of these medians over all transcripts of its gene.
  So the weights of a gene add up to 1 in each tissue. The output has the key columns and a struct
  column `features` with 28 fields: 12 from the scores and the PTC counts, and a weighted proportion and
  a maximum for each of 8 NMD-Scanner flags (start and stop loss, the 5 NMD escape rules and
  `ptc_less_than_150nt_to_start`). `alt_has_ptc.proportion` and `num_ptc` count every transcript in which
  the variant creates a PTC, also those that the model cannot score. As in tissue_specific_vep, each variant
  and gene gets a row for every tissue of the isoform proportions. A transcript that GTEx does not list under
  its gene has no weight, so it counts only in the features without weights, such as `num_ptc` and
  `nmd_pred_score`.

The module is off by default. Set `system.nmd_scanner.enabled: true` and request
`<output_dir>/veff/nmd_scanner/features.parquet/<vcf_file>.parquet` as target.

The loftee and nmd_scanner modules keep the output of their tool as it is. Neither module is read by
tissue_specific_vep, fset or predict.

`workflow/modules/veff/absplice2` adds [AbSplice2-DNA](https://github.com/gagneurlab/absplice2).
It runs [Pangolin](https://github.com/tkzeng/Pangolin) and scores each variant, gene and GTEx tissue
with the AbSplice2 model, from Pangolin and from the MMSplice and SpliceMap results of the absplice
module. It is off by default; set `system.absplice2.enabled: true` and request
`<output_dir>/veff/absplice2/veff.parquet/<vcf_file>.parquet` as target. Its output is not read by
tissue_specific_vep, fset or predict. The module downloads the AbSplice2 model (0.7 MB).

Pangolin needs a GPU for large VCFs (`use_gpu: True`, see [GPU](#gpu)). On a CPU with 4 threads, it takes
about 3 s per variant, i.e. about 9 hours for 10,000 variants. A whole-genome VCF with millions of variants
would take months.

The module builds the Pangolin annotation database from `gff3_file`: all genes, and the transcripts and
exons with the transcript tags of the databases published with Pangolin (Ensembl_canonical for hg38).
AbSplice2 was trained with Pangolin's GENCODE v38 database. Built from the GENCODE v38 GFF3, the database
has the same genes, transcripts and exons. Other GENCODE releases differ in their genes and canonical
transcripts. For every 4th variant of the chr22 ClinVar example VCFs, GENCODE v40 instead of v38 added 14
variant/gene pairs, all with AbSplice_DNA below 0.001, and did not change the other 2,603 pairs. Set
`system.absplice2.pangolin_annotation_db` to use an existing database instead, e.g. the GENCODE v38
database of Pangolin.

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
TAG = "v2.0.0"  # x-release-please-version
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
    # optional: run VEP without CADD or LOFTEE, and skip their downloads
    # "cadd_plugin": False,
    # "loftee_plugin": False,
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
- The rerun triggers of AbExp's workflow profile do not apply in the importing workflow, see
  [When jobs rerun](#when-jobs-rerun).

## When jobs rerun

The workflow profile `workflow/profiles/default` sets `rerun-triggers: [mtime, params, input]`.
With it, Snakemake reruns a job only if an input file is newer than its outputs, if its set of input files
changed, or if its params changed. A changed conda environment or changed code does not rerun jobs.
Snakemake's default triggers would rerun them after any change of a conda environment file, even of a comment,
and whenever a script is newer than the outputs, e.g. after a `git pull` or a fresh clone.

Instead, each rule that computes its outputs has the param `output_version`. Its value is one entry of the dict
`OUTPUT_VERSION` in the Snakefile of the module, or in `workflow/Snakefile` for the feature sets and predictions.
The dict has one key per step. A step is one rule or several rules that must change together, e.g. the Enformer
predictions of the reference and the alternative sequences. A new value reruns the rules of the step, and the
rules downstream of them rerun because their input changed.
- In a commit that changes outputs, bump the key of the earliest step whose outputs change, e.g. after a new tool
  version in a conda environment or a changed script. For example, if a new aberrant-expression release changes
  the scores of the Enformer tissue mapper, bump `"tissue"` in `workflow/modules/veff/enformer/Snakefile`. The
  tissue rules and all rules downstream of them rerun, but not the Enformer predictions.
- A conda environment can serve several steps, e.g. the TensorFlow environment of Enformer and SpliceAI. After a
  change of such an environment, bump each step whose outputs change.
- Do not bump a key for changes that keep the outputs, e.g. a comment.
- Snakemake deletes temporary outputs once no job needs them. So a bump of a step that reads a temporary output
  also reruns the step that wrote it. Therefore the MMSplice, SpliceAI and Pangolin scores and the aggregated
  and tissue Enformer predictions are not temporary: storing them costs less than recomputing them. The raw
  Enformer predictions and the VEP output stay temporary. So the Enformer aggregation shares the key
  `"predict"`, and the VEP annotation and its parsing share one key.
- The download rules have no version; their URL is the version. Rules that only index a file or convert the
  format of a download have none either.

Snakemake uses this profile when it runs `workflow/Snakefile`, also with `--snakefile` from another directory.
`--rerun-triggers` on the command line overrides the setting. The setting does not apply:
- with `--workflow-profile none` or another workflow profile
- in a workflow that imports AbExp as a module, which uses its own profile. Pass
  `--rerun-triggers mtime params input` there to get the same behavior.

The recorded params of outputs from earlier AbExp versions lack `output_version`. So the first run after the
update reruns all rules with an `output_version` and the rules downstream of them, once. To keep the existing
outputs instead, run Snakemake once with `--touch` and the same config and targets. It marks the outputs as up
to date and records the new params. Do this only if the last run finished and nothing else changed since,
because `--touch` also marks outputs as up to date that need a rerun. Alternatively,
`--cleanup-metadata <files>` deletes the recorded params of the given outputs.

## Python package

`packages/aberrant_expression` contains the Python distribution aberrant-expression, which the workflow imports as
`abexp`. Its subpackages are `abexp.utils`, `abexp.enformer`, and `abexp.mmsplice`, `abexp.absplice` and
`abexp.spliceai_rocksdb`: the parts of MMSplice, AbSplice and SpliceAI-RocksDB that the absplice module runs, ported
to kipoiseq2. The conda environments install the package from this repository, pinned to its release tag. So a
change of the package takes effect in the workflow only with its next release.
[packages/README.md](packages/README.md) describes how to test, release and publish it.

## Releases

release-please releases the workflow from the conventional commits on main, see
`.github/workflows/release-please.yml`. The release PR updates `CHANGELOG.md` and `TAG` in the module example
above. Merging it creates the tag `vX.Y.Z` and the GitHub release. Commits that change only files in `packages/`
do not count for the workflow, because the package has its own releases.

## License
All source code and model weights in this repository are licensed under the [MIT license](./LICENSE).
The Python package in `packages/` carries its own MIT license file.

**Please note:** AbExp relies on [CADD](https://cadd.gs.washington.edu/) and [SpliceAI](https://github.com/Illumina/SpliceAI/), both of which are free to use only in non-commercial settings.
If you plan to use AbExp in a commercial context, please ensure that you have the appropriate permissions or licenses to use both tools.

The optional module absplice2 installs [Pangolin](https://github.com/tkzeng/Pangolin) and downloads the model of
[AbSplice2](https://github.com/gagneurlab/absplice2) when it runs. Both are licensed under the GPL-3.0 and are not part of this repository.

## Development setup
Advanced users who want to edit this pipeline can use the following steps to convert the python scripts back to Jupyter notebooks:
1) Make sure that the `jupytext` command is available, e.g. via `mamba install jupytext`
2) run `find workflow/ -iname "*[.py.py|.R.R]" -exec jupytext --sync {} \;` to convert all percent scripts to jupyter notebooks
Jupyter will then automatically synchronize the percent scripts with the corresponding notebook files.

