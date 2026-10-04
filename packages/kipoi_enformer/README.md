# Kipoi-enformer
Variant effect prediction of promoter variants using the [Enformer](https://github.com/google-deepmind/deepmind-research/tree/master/enformer) model.

## Installation

Install a release from the AbExp repository with pip:

<!-- x-release-please-start-version -->
```bash
pip install "kipoi-enformer @ git+https://github.com/gagneurlab/AbExp.git@kipoi-enformer-v0.0.1#subdirectory=packages/kipoi_enformer"
```
<!-- x-release-please-end -->

For development, install the package in editable mode from a checkout of AbExp:

```bash
pip install -e "packages/kipoi_enformer[dev]"
```

For GPU support on Linux, also install `tensorflow[and-cuda]`.

## Enformer model

`Enformer()` loads the Enformer model from Kaggle Models with kagglehub, under the handle
`deepmind/enformer/tensorFlow2/enformer/1`. These are the same files that TF Hub served at
`https://tfhub.dev/deepmind/enformer/1`. On first use, kagglehub downloads the model, about 1 GB, and caches it in
`~/.cache/kagglehub`. Set `KAGGLEHUB_CACHE` to use another directory. Once the model is cached, kagglehub loads it
without contacting Kaggle. Version 0.0.1 downloaded the model with tensorflow-hub instead and cached it in
`TFHUB_CACHE_DIR`.

## Tests

Most tests read example files: the files in `tests/data/` of this package and the chr22 sequence in
`example/chr22_hg19.fa` of the AbExp repository. git-lfs stores the files in `tests/data/`. `.lfsconfig` excludes
them from the default fetch, so that pip installs from a git URL do not download them. With git-lfs installed, fetch
them in a checkout of AbExp and run the tests:

```bash
git lfs pull --include="packages/kipoi_enformer/tests/data/**" --exclude=""
pytest packages/kipoi_enformer/tests
```

Pytest skips a test if one of its files is missing or not fetched, e.g. in an sdist.

The tests run across all cores by default, through pytest-xdist. Pass `-n0` to run them in one process, which a
debugger needs and which restores per-test output order.
The tests that run Enformer share their outputs, so one worker runs them, one after the other.

With `-n0`, the tests need about 3 GB of memory.

## Genome annotation

The dataloaders and `EnformerVeff` take the genome annotation as `genome_annotation`: the path to a GFF3 file, such
as a GENCODE annotation, or a polars or pandas DataFrame with pyranges-style columns, see
`kipoi_enformer.utils.genome_annotation_to_polars`. Version 0.0.1 read a GTF file instead, and the parameter was
`gtf`.

## Tissue mapper

`EnformerTissueMapper` maps the Enformer CAGE tracks to an expression score per GTEx tissue, with one linear model per
tissue. `train` fits a scikit-learn pipeline of a `StandardScaler` and a linear model for each tissue. The model must
be linear, e.g. from `sklearn.linear_model`; other models raise a `TypeError`. `train` writes the parameters of the
pipelines to a parquet file with one row per tissue: `tissue`, the `mean` and `scale` of the `StandardScaler`, and
the `coef` and `intercept` of the linear model. The features are the tracks in the order of the tracks yaml file.
`predict` computes the scores from this file with numpy and gives the same scores as the pipelines with scikit-learn
1.5 to 1.7.

## Usage

The example reads the chr22 sequence and genome annotation in `example/` of the AbExp repository and the test
data of this package. Run it in the root of a checkout of AbExp, after you have fetched the test data, see
[Tests](#tests).

```python
from kipoi_enformer.dataloader import RefTSSDataloader, VCFTSSDataloader
from kipoi_enformer.enformer import Enformer, EnformerAggregator, EnformerTissueMapper, EnformerVeff
from pathlib import Path
from sklearn import linear_model

fasta_file = 'example/chr22_hg19.fa'
genome_annotation = 'example/chr22.gencode.v40lift37.annotation.gff3.gz'
data_dir = Path('packages/kipoi_enformer/tests/data')

# define output dirs
output_dir = Path('output')
(output_dir / 'raw/ref.parquet/chrom=chr22').mkdir(exist_ok=True, parents=True)
(output_dir / 'raw/alt.parquet').mkdir(exist_ok=True, parents=True)
(output_dir / 'aggregated/ref.parquet/chrom=chr22').mkdir(exist_ok=True, parents=True)
(output_dir / 'aggregated/alt.parquet').mkdir(exist_ok=True, parents=True)
(output_dir / 'tissue/ref.parquet/chrom=chr22').mkdir(exist_ok=True, parents=True)
(output_dir / 'tissue/alt.parquet').mkdir(exist_ok=True, parents=True)

# define enformer objects
enformer = Enformer()
enformer_aggregator = EnformerAggregator()
enformer_tissue_mapper = EnformerTissueMapper(
    tracks_path=data_dir / 'human_cage_nonuniversal_enformer_tracks.yaml')
enformer_veff = EnformerVeff(genome_annotation=genome_annotation)

# Reference sequences
# define reference dataloader
ref_dl = RefTSSDataloader(fasta_file=fasta_file, genome_annotation=genome_annotation, chromosome='chr22',
                          canonical_only=False, protein_coding_only=True)
# run enformer on reference genome
enformer.predict(ref_dl, batch_size=2, filepath=output_dir / 'raw/ref.parquet/chrom=chr22/data.parquet',
                 num_output_bins=11)
# aggregate reference enformer scores
enformer_aggregator.aggregate(output_dir / 'raw/ref.parquet/chrom=chr22/data.parquet',
                              output_dir / 'aggregated/ref.parquet/chrom=chr22/data.parquet')
# train tissue mapper using reference genome
enformer_tissue_mapper.train([output_dir / 'aggregated/ref.parquet/chrom=chr22/data.parquet'],
                             output_path=output_dir / 'tissue_mapper.parquet',
                             expression_path=data_dir / 'transcripts_tpms.zarr',
                             model=linear_model.ElasticNetCV(cv=2))
# map reference to tissues
enformer_tissue_mapper.predict(output_dir / 'aggregated/ref.parquet/chrom=chr22/data.parquet',
                               output_dir / 'tissue/ref.parquet/chrom=chr22/data.parquet')

# Alternative sequences
# define alternative dataloader
alt_dl = VCFTSSDataloader(fasta_file=fasta_file, genome_annotation=genome_annotation,
                          vcf_file=data_dir / 'chr22_var.vcf.gz', variant_upstream_tss=50,
                          variant_downstream_tss=200, canonical_only=False, protein_coding_only=True)
# run enformer on alternative genome
enformer.predict(alt_dl, batch_size=2, filepath=output_dir / 'raw/alt.parquet/chr22_var.vcf.gz.parquet',
                 num_output_bins=11)
# aggregate alternative enformer scores
enformer_aggregator.aggregate(output_dir / 'raw/alt.parquet/chr22_var.vcf.gz.parquet',
                              output_dir / 'aggregated/alt.parquet/chr22_var.vcf.gz.parquet')
# map alternative to tissues
enformer_tissue_mapper.predict(output_dir / 'aggregated/alt.parquet/chr22_var.vcf.gz.parquet',
                               output_dir / 'tissue/alt.parquet/chr22_var.vcf.gz.parquet')

# Variant effect prediction
enformer_veff.run(ref_paths=[output_dir / 'tissue/ref.parquet/chrom=chr22/data.parquet'],
                  alt_path=output_dir / 'tissue/alt.parquet/chr22_var.vcf.gz.parquet',
                  output_path=output_dir / 'veff.parquet', aggregation_mode='canonical')
```