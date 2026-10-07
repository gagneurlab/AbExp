# aberrant-expression

tl;dr: The Python packages of the [AbExp](https://github.com/gagneurlab/AbExp) pipeline, in one distribution. pip
installs it as `aberrant-expression`, and Python imports it as `abexp`. The base dependencies cover only
`abexp.utils.common` and `abexp.utils.models`; each extra adds the dependencies of other subpackages.

| import                   | contents                                                                                  | extra      |
| ------------------------ | ----------------------------------------------------------------------------------------- | ---------- |
| `abexp.utils`            | functions on polars and Spark DataFrames, and the scikit-learn model wrappers of AbExp     | `polars`, `spark` |
| `abexp.enformer`         | variant effect prediction of promoter variants with Enformer                              | `enformer` |
| `abexp.mmsplice`         | MMSplice and the junction VCF dataloaders                                                 | `splicing` |
| `abexp.absplice`         | AbSplice-DNA from MMSplice with SpliceMaps and from SpliceAI                              | `splicing` |
| `abexp.spliceai_rocksdb` | SpliceAI scores from SpliceAI-RocksDB                                                     | `rocksdb`  |

`abexp.utils` and `abexp.enformer` were the distributions abexp-utils (`abexp_utils`) and kipoi-enformer
(`kipoi_enformer`). `abexp.mmsplice`, `abexp.absplice` and `abexp.spliceai_rocksdb` were the distribution
abexp-splicing. The modules keep their paths below the new top level, e.g. `kipoi_enformer.enformer` is now
`abexp.enformer.enformer`.

The name `abexp` on PyPI belongs to an unrelated library, which also installs into `site-packages/abexp`. So the two
cannot share an environment.

## Installation

Install a release from the AbExp repository with pip. Choose the extras you need:

<!-- x-release-please-start-version -->
```bash
pip install "aberrant-expression[polars] @ git+https://github.com/gagneurlab/AbExp.git@aberrant-expression-v0.0.1#subdirectory=packages/aberrant_expression"
```
<!-- x-release-please-end -->

| extra      | for                                                                       |
| ---------- | ------------------------------------------------------------------------- |
| `polars`   | `abexp.utils.polars_functions`                                            |
| `spark`    | `abexp.utils.spark_functions`. The Spark functions need a Java runtime; pyspark 4 requires Java 17 or newer. |
| `splicing` | `abexp.mmsplice` and `abexp.absplice`                                     |
| `rocksdb`  | `abexp.spliceai_rocksdb`, with `splicing`                                 |
| `enformer` | `abexp.enformer`. For GPU support on Linux, also install `tensorflow[and-cuda]`. |
| `all`      | all extras but `rocksdb`                                                  |
| `dev`      | pytest, pytest-cov and pytest-xdist for the tests                         |

The extra `rocksdb` installs python-rocksdb, which pip builds only with the RocksDB library installed; conda-forge
has a build of it. SpliceAI predictions of the variants that are not in SpliceAI-RocksDB also need SpliceAI.
Install it separately. The AbExp pipeline uses the SpliceAI fork [hoeze/SpliceAI](https://github.com/hoeze/SpliceAI),
which also predicts multi-nucleotide variants:

```bash
pip install "spliceai @ git+https://github.com/hoeze/SpliceAI.git@e3470ae185d1418fe40cf35e8c4e33431cbc5df2"
```

For development, install the package in editable mode from a checkout of AbExp, and add `rocksdb` to the extras if
you work on `abexp.spliceai_rocksdb`:

```bash
pip install -e "packages/aberrant_expression[all,dev]"
```

## abexp.utils

`abexp.utils.polars_functions` and `abexp.utils.spark_functions` reshape polars and Spark DataFrames and join the
feature sets of AbExp. `abexp.utils.models` holds the scikit-learn model wrappers of the AbExp models.

## abexp.enformer

### Enformer model

`Enformer()` loads the Enformer model from Kaggle Models with kagglehub, under the handle
`deepmind/enformer/tensorFlow2/enformer/1`. These are the same files that TF Hub served at
`https://tfhub.dev/deepmind/enformer/1`. On first use, kagglehub downloads the model, about 1 GB, and caches it in
`~/.cache/kagglehub`. Set `KAGGLEHUB_CACHE` to use another directory. Once the model is cached, kagglehub loads it
without contacting Kaggle. kipoi-enformer 0.0.1 downloaded the model with tensorflow-hub instead and cached it in
`TFHUB_CACHE_DIR`.

### Genome annotation

The dataloaders and `EnformerVeff` take the genome annotation as `genome_annotation`: the path to a GFF3 file, such
as a GENCODE annotation, or a polars or pandas DataFrame with pyranges-style columns, see
`abexp.enformer.utils.genome_annotation_to_polars`. kipoi-enformer 0.0.1 read a GTF file instead, and the parameter
was `gtf`. `gtf` still works as a deprecated alias and emits a `DeprecationWarning`, but a path must now point to a
GFF3 file.

### Tissue mapper

`EnformerTissueMapper` maps the Enformer CAGE tracks to an expression score per GTEx tissue, with one linear model per
tissue. `train` fits a scikit-learn pipeline of a `StandardScaler` and a linear model for each tissue. The model must
be linear, e.g. from `sklearn.linear_model`; other models raise a `TypeError`. `train` writes the parameters of the
pipelines to a parquet file with one row per tissue: `tissue`, the `mean` and `scale` of the `StandardScaler`, and
the `coef` and `intercept` of the linear model. The features are the tracks in the order of the tracks yaml file.
`predict` computes the scores from this file with numpy and gives the same scores as the pipelines with scikit-learn
1.5 to 1.7.

### Usage

The example reads the chr22 sequence and genome annotation in `example/` of the AbExp repository and the test
data of this package. Run it in the root of a checkout of AbExp, after you have fetched the test data, see
[Tests](#tests).

```python
from abexp.enformer.dataloader import RefTSSDataloader, VCFTSSDataloader
from abexp.enformer.enformer import Enformer, EnformerAggregator, EnformerTissueMapper, EnformerVeff
from pathlib import Path
from sklearn import linear_model

fasta_file = 'example/chr22_hg19.fa'
genome_annotation = 'example/chr22.gencode.v40lift37.annotation.gff3.gz'
data_dir = Path('packages/aberrant_expression/tests/enformer/data')

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

## abexp.mmsplice, abexp.absplice and abexp.spliceai_rocksdb

The parts of MMSplice, AbSplice and SpliceAI-RocksDB that the AbExp pipeline runs, on kipoiseq2. kipoi, kipoiseq 0.7,
cyvcf2 and pyranges are no longer needed. The outputs match the upstream packages, except for the row order,
rounding in the last digits and the choice among tied junctions, see
[Differences](#differences-from-the-upstream-packages).

| import                    | ported from                                                                                     | what AbExp uses                                                     |
| ------------------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| `abexp.mmsplice`          | [mmsplice](https://github.com/gagneurlab/MMSplice_MTSplice) 2.4.0 (31513da)                     | `MMSplice` and the junction VCF dataloaders                         |
| `abexp.absplice`          | [absplice](https://github.com/gagneurlab/absplice) daad7b6, [splicemap](https://github.com/gagneurlab/splicemap) cf922eb | `SpliceOutlierDataloader`, `SpliceOutlier`, `SplicingOutlierResult`, `read_spliceai_vcf` |
| `abexp.spliceai_rocksdb`  | [spliceai_rocksdb](https://github.com/gagneurlab/spliceai_rocksdb) 3c40d6e                      | `SpliceAI`                                                          |

### Usage

The example runs AbSplice-DNA on the test data of this package. Run it in the root of a checkout of AbExp.

```python
import polars as pl
from abexp.absplice import SpliceOutlier, SpliceOutlierDataloader, SplicingOutlierResult

data_dir = 'packages/aberrant_expression/tests/splicing/data'

# MMSplice with SpliceMaps: delta PSI per variant, junction and tissue
dl = SpliceOutlierDataloader(
    'example/chr22_hg38.fa', f'{data_dir}/clinvar_chr22.vcf',
    splicemap5=[f'{data_dir}/Whole_Blood_splicemap_psi5.csv.gz'],
    splicemap3=[f'{data_dir}/Whole_Blood_splicemap_psi3.csv.gz'],
)
SpliceOutlier().predict_save(dl, 'mmsplice_splicemap.csv')

# AbSplice-DNA per variant, gene and tissue, from MMSplice and SpliceAI
result = SplicingOutlierResult(
    df_mmsplice=pl.read_csv('mmsplice_splicemap.csv'),
    df_spliceai=pl.read_csv(f'{data_dir}/expected/spliceai_vcf.csv'),
)
df = result.predict_absplice_dna()
```

`abexp.spliceai_rocksdb.SpliceAI(fasta, annotation='grch38', db_path={'22': <path>}).predict_save(vcf, csv)`
looks up the SpliceAI scores in SpliceAI-RocksDB and runs SpliceAI for the variants that are not in it.
`abexp.absplice.read_spliceai_vcf` reads a VCF file that the SpliceAI command line tool annotated.

### Differences from the upstream packages

- kipoiseq2 yields the variant-junction pairs in the order of the VCF file, and kipoiseq 0.7 in the order of
  pyranges. So the rows of the MMSplice table come in another order.
- MMSplice therefore predicts batches with other samples. TensorFlow rounds them differently, so `delta_logit_psi`
  can differ in the last digits.
- A variant can have the same delta PSI at several junctions of a gene. absplice daad7b6 picked any one of them,
  and the row order and the numpy version decided which one. abexp.absplice picks the first in the order of
  junction, event_type, splice_site and the other columns. So the reported junction can differ from absplice, and
  with its `median_n` also `splice_site_is_expressed` and `AbSplice_DNA`.
- `read_spliceai_vcf` names the column `acceptor_loss_position`, not `acceptor_loss_positiin`.
- `read_spliceai_vcf` gives each ALT allele only the SpliceAI entries of that allele. absplice gave it the entries
  of all ALT alleles of its record. The workflow splits multi-allelic records before SpliceAI, so its tables do
  not change.
- abexp.mmsplice, abexp.absplice and abexp.spliceai_rocksdb use polars and numpy instead of pandas. The tables in
  and out are polars DataFrames. `predict_absplice_dna` returns variant, gene_id and tissue as columns, sorted by
  them, and not as an index.
- The CSV files of `predict_save` format floats as polars does, e.g. `6.4e-6` instead of `6.4e-06`. The values do
  not change.
- A VCF file without variants yields no samples. mmsplice compared the contigs of the VCF header with the FASTA
  file, and abexp.mmsplice compares the chromosomes of the variants.
- Only the parts that AbExp uses are kept. Not included: MTSplice, the VEP plugin, the dataloaders of GTF exons
  and exon tables, the VCF writers, the variant filters of absplice, AbSplice-RNA, the CADD-Splice and
  sample-based features, and the SpliceMap count tables. `predict_absplice_dna` reads only ONNX models.

## Tests

The tests are in one directory per area, and each area needs the extras of its subpackages:

| directory        | extras            | reads                                                                                   |
| ---------------- | ----------------- | --------------------------------------------------------------------------------------- |
| `tests/utils`    | `polars`, `spark` | nothing; the Spark tests need a Java runtime                                            |
| `tests/enformer` | `enformer`        | the files in `tests/enformer/data/` and the hg19 chr22 sequence `example/chr22_hg19.fa` |
| `tests/splicing` | `splicing`        | the files in `tests/splicing/data/` and the hg38 chr22 sequence `example/chr22_hg38.fa` |

The example sequences are in the AbExp repository, outside of the package. Pytest skips a test if one of its files
is missing, e.g. in an sdist.

git-lfs stores the files in `tests/enformer/data/`. `.lfsconfig` excludes them from the default fetch, so that pip
installs from a git URL do not download them. With git-lfs installed, fetch them in a checkout of AbExp and run the
tests:

```bash
git lfs pull --include="packages/aberrant_expression/tests/enformer/data/**" --exclude=""
cd packages/aberrant_expression
pytest tests/utils
pytest tests/enformer
pytest tests/splicing
```

The tests run across all cores by default, through pytest-xdist. Pass `-n0` to run them in one process, which a
debugger needs and which restores per-test output order. The tests of one xdist group run on one worker, one after
the other: the Spark tests, and the tests that run Enformer, which share their outputs. With `-n0`, the enformer
tests need about 3 GB of memory.

The expected outputs in `tests/splicing/data/expected/` come from the upstream packages, see
`tests/splicing/make_expected.py`. The SpliceAI-RocksDB tests also need the extra `rocksdb` and the hg38 chr22
database (2.6 GB). Set `ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22` to the path of `spliceAI_hg38_chr22.db`, otherwise pytest
skips them. The test of the SpliceAI predictions also needs SpliceAI.

## License

The code is under the MIT license, see [LICENSE](LICENSE). It lists the copyright notices of the upstream
projects of `abexp.mmsplice`, `abexp.absplice` and `abexp.spliceai_rocksdb`.
