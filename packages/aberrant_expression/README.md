# aberrant-expression

tl;dr: The Python packages of the [AbExp](https://github.com/gagneurlab/AbExp) pipeline, in one distribution. pip
installs it as `aberrant-expression`, and Python imports it as `abexp`. The base has no dependencies and covers
`abexp.utils.common`; each extra adds the dependencies of other subpackages.

| import                   | contents                                                                                  | extra      |
| ------------------------ | ----------------------------------------------------------------------------------------- | ---------- |
| `abexp.utils`            | functions on polars and Spark DataFrames, a GFF3 reader, and the scikit-learn model wrappers of AbExp | `polars`, `spark`, `gff3`, `models` |
| `abexp.enformer`         | variant effect prediction of promoter variants with Enformer                              | `enformer` |
| `abexp.mmsplice`         | MMSplice and the junction VCF dataloaders                                                 | `splicing` |
| `abexp.absplice`         | AbSplice-DNA from MMSplice with SpliceMaps and from SpliceAI                              | `splicing` |
| `abexp.spliceai_rocksdb` | SpliceAI scores from SpliceAI-RocksDB                                                     | `rocksdb`  |
| `abexp.pangolin`         | Pangolin's splice scores per variant and gene, with the published weights of Pangolin     | `pangolin` |
| `abexp.absplice2`        | AbSplice2-DNA from Pangolin and from MMSplice with SpliceMaps                             | `absplice2` |

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
pip install "aberrant-expression[polars] @ git+https://github.com/gagneurlab/AbExp.git@aberrant-expression-v0.1.0#subdirectory=packages/aberrant_expression"
```
<!-- x-release-please-end -->

| extra      | for                                                                       |
| ---------- | ------------------------------------------------------------------------- |
| `models`   | `abexp.utils.models`                                                      |
| `polars`   | `abexp.utils.polars_functions`                                            |
| `spark`    | `abexp.utils.spark_functions`. The Spark functions need a Java runtime; pyspark 4 requires Java 17 or newer. |
| `gff3`     | `abexp.utils.gff3`                                                        |
| `splicing` | `abexp.mmsplice` and `abexp.absplice`                                     |
| `rocksdb`  | `abexp.spliceai_rocksdb`, with `splicing`                                 |
| `pangolin` | `abexp.pangolin`. conda-forge has PyTorch as `pytorch-cpu` and `pytorch-gpu`. |
| `absplice2` | `abexp.absplice2`                                                        |
| `enformer` | `abexp.enformer`. For GPU support on Linux, also install `tensorflow[and-cuda]`. |
| `all`      | all extras but `rocksdb`                                                  |
| `numpy`, `kipoiseq2`, `tensorflow`, `tqdm` | one requirement each, which several extras share |

The extra `rocksdb` installs python-rocksdb, which pip builds only with the RocksDB library installed; conda-forge
has a build of it. SpliceAI predictions of the variants that are not in SpliceAI-RocksDB also need SpliceAI.
Install it separately. The AbExp pipeline uses the SpliceAI fork [hoeze/SpliceAI](https://github.com/hoeze/SpliceAI),
which also predicts multi-nucleotide variants:

```bash
pip install "spliceai @ git+https://github.com/hoeze/SpliceAI.git@e3470ae185d1418fe40cf35e8c4e33431cbc5df2"
```

For development, use the uv workspace of a checkout of AbExp. In the root of the checkout, this installs the
package in editable mode, with all extras but `rocksdb` and with the dependency group `test` for the tests:

```bash
uv sync --extra all --group test
```

Add `--extra rocksdb` if you work on `abexp.spliceai_rocksdb`.
[packages/README.md](https://github.com/gagneurlab/AbExp/blob/main/packages/README.md) in the AbExp repository
describes the workspace and the tests.

## abexp.utils

`abexp.utils.polars_functions` and `abexp.utils.spark_functions` reshape polars and Spark DataFrames and join the
feature sets of AbExp. `abexp.utils.models` holds the scikit-learn model wrappers of the AbExp models, and needs the
extra `models`. `abexp.utils.gff3` reads GFF3 files with polars-bio for `abexp.enformer` and `abexp.pangolin`, and
needs the extra `gff3`. With `par_y_suffix=True`, it adds the suffix `_PAR_Y` to `gene_id` and `transcript_id` of the
chrY PAR copies, as `abexp.enformer` does.

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
was `gtf`.

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

### Known issues

`EnformerAggregator` averages the bins around the TSS bin of each shift. As the TSS, it takes the middle of the bins
that `Enformer.predict` saved. That is right for an even `num_output_bins`, such as all 896. With 21 output bins, the
workflow default, Enformer saves the bins 438 to 458. Their middle lies 64 bp downstream of the TSS. So shift +43
averages the bins 447 to 449, although the TSS lies in bin 447. With 896 output bins, it averages 446 to 448. Shifts
0 and -43 average 447 to 449 in both cases.

The training pipeline [gtsitsiridis/kipoi_veff_analysis](https://github.com/gtsitsiridis/kipoi_veff_analysis)
trained the published tissue mapper (`elasticnet_cage_gtexv8`) on 21 output bins, so the tissue mapper learned this
behavior. The reference scores that the workflow downloads come from the same pipeline. So the code keeps the
behavior. A fix would replace `pred_seq_length // 2` with `pred.shape[2] // 2 * bin_size` in
`EnformerAggregator._aggregate_batch`. It needs a retrained tissue mapper and recomputed reference scores. The offset
came in with the crop to 21 output bins, in commit e2bbf5ec (2024-05-16) of the training pipeline.

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
- If no variant lies on a chromosome of the junctions, e.g. in a VCF file with only chrM variants, abexp.mmsplice
  logs a warning and yields no samples. mmsplice raised a ValueError if no contig of the VCF header had junctions.
- Only the parts that AbExp uses are kept. Not included: MTSplice, the VEP plugin, the dataloaders of GTF exons
  and exon tables, the VCF writers, the variant filters of absplice, AbSplice-RNA, the CADD-Splice and
  sample-based features, and the SpliceMap count tables. `predict_absplice_dna` reads only ONNX models.

## abexp.pangolin

Pangolin's splice scores of variants (Zeng and Li, Genome Biology 2022, [Pangolin](https://github.com/tkzeng/Pangolin)),
on kipoiseq2 and PyTorch. Per variant and gene, it gives the largest gain and loss of splice site usage within
`distance` bases, like `pangolin -d 50 -m True`. The network and the data path are our own code; only the weights
come from Pangolin. The network uses unpadded convolutions, which give the same outputs with about a third less
computation, and it scores ref and alt sequences of equal length in batches.

The weights are not part of this package. They are GPL-3, as the Pangolin repository. `download_models` downloads
the 12 files of `abexp.pangolin.MODEL_FILES` (35 MB) from `pangolin/models` of the Pangolin repository at the
pinned commit `PANGOLIN_COMMIT` (5cf94b8), and checks their SHA-256 sums. It keeps the files that exist with the right
sum. `download_model(file, models_dir)` does the same for one file. `PangolinModels.from_dir` checks the SHA-256
sums too, and stops with a ValueError on files of another commit.

```python
from abexp.pangolin import Pangolin, PangolinModels, download_models, read_gff3_genes

download_models('pangolin_models')
genes = read_gff3_genes('gencode.v42.annotation.gff3.gz', transcript_tags=['Ensembl_canonical'])
pangolin = Pangolin('genome.fa', genes, PangolinModels.from_dir('pangolin_models'), distance=50, mask=True)
pangolin.predict_df('variants.vcf').write_parquet('pangolin.parquet')
```

### Differences from upstream Pangolin

The scores match Pangolin (fork neverov-am/Pangolin 232cba0) up to the last float32 digits, except here:

- Every ALT allele of a record is scored. Pangolin scores only the first one.
- A gene counts if the ref allele overlaps it. Pangolin misses a gene that starts at the variant position, and a gene
  that a deletion reaches only after its first base.
- Beyond the chromosome ends, the sequence is N. Pangolin skips variants within 5050 bases of the chromosome start,
  and stops with an error for most variants within 5050 bases of its end.
- The ref check ignores the case of the FASTA file, and the sequence is upper case. Pangolin skips a variant whose
  FASTA ref is in lower case.
- Variants with letters other than A, C, G, T and N in ref or alt are skipped. Pangolin stops with an error on them.
  kipoiseq2 drops ALT alleles with N; Pangolin scores them with N as 0.
- The output is a polars DataFrame with unrounded float32 scores, not a VCF file with scores rounded to 2 decimals.
- The genes and splice sites come from a GENCODE GFF3 file, not from a gffutils database. An exon belongs to the
  gene with its `gene_id` on its chromosome.
- With `mask`, each gene is masked on its own, as in the fork. tkzeng/Pangolin 5cf94b8 masks the loss and gain
  arrays of a strand in place. So there, each gene starts from the arrays that the genes before it on the same
  strand have masked (`pangolin/pangolin.py`, lines 133 to 153). The fork copies the arrays for each gene (lines 144
  and 145 of 232cba0).

## abexp.absplice2

AbSplice2-DNA per variant, gene and GTEx tissue, from Pangolin, MMSplice with SpliceMaps, and the SpliceMaps. The
code is our own, in polars, for the steps of the [AbSplice2](https://github.com/gagneurlab/absplice2) example
workflow after MMSplice and Pangolin, at commit a30120f. `absplice2_dna` matches Pangolin's gain and loss sites
with the SpliceMap sites within 2 bp, joins them with MMSplice, scores each row with the model, and keeps the row
with the largest score per variant, gene and tissue.

The AbSplice2 model is not part of this package. It is an ExplainableBoostingClassifier of interpret 0.2.7, which
the caller loads and passes as a function:

```python
import pickle

import polars as pl
from abexp.absplice2 import absplice2_dna, read_mmsplice_splicemap

with open('AbSplice_2_DNA.pkl', 'rb') as fd:
    model = pickle.load(fd)
df = absplice2_dna(
    pangolin=pl.read_parquet('pangolin.parquet'),  # the output of abexp.pangolin
    splicemap5=['Whole_Blood_splicemap_psi5.csv.gz'],
    splicemap3=['Whole_Blood_splicemap_psi3.csv.gz'],
    mmsplice_splicemap=read_mmsplice_splicemap('mmsplice_splicemap.csv'),
    predict=lambda features: model.predict_proba(features.to_pandas())[:, 1],
)
```

The output matches the AbSplice2 scripts, except here:

- The variant is split into `chrom`, `start` (0-based), `end`, `ref` and `alt`, and `gene_id` is named `gene`.
- Pangolin's scores are rounded to 2 decimals in float64, Pangolin rounds them in float32. The two differ by 0.01
  only for scores within float32 precision of a rounding boundary.
- Of the rows with the same largest score, AbSplice2 keeps any one. `absplice2_dna` keeps the first in the order
  of the output columns, so its MMSplice columns can differ from AbSplice2.

## Tests

The tests are in one directory per area, and each area needs the extras of its subpackages:

| directory        | extras            | reads                                                                                   |
| ---------------- | ----------------- | --------------------------------------------------------------------------------------- |
| `tests/utils`    | `polars`, `spark`, `gff3` | nothing; the Spark tests need a Java runtime                                    |
| `tests/enformer` | `enformer`        | the files in `tests/enformer/data/` and the hg19 chr22 sequence `example/chr22_hg19.fa` |
| `tests/splicing` | `splicing`        | the files in `tests/splicing/data/` and the hg38 chr22 sequence `example/chr22_hg38.fa` |
| `tests/pangolin` | `pangolin`        | the files in `tests/pangolin/data/` and the published weights of Pangolin (35 MB)       |
| `tests/absplice2` | `absplice2`      | nothing                                                                                 |

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
pytest tests/pangolin
pytest tests/absplice2
```

The tests run across all cores by default, through pytest-xdist. Pass `-n0` to run them in one process, which a
debugger needs and which restores per-test output order. The tests of one xdist group run on one worker, one after
the other: the Spark tests, and the tests that run Enformer, which share their outputs. With `-n0`, the enformer
tests need about 3 GB of memory.

The pangolin tests run small networks with fixed random weights on a synthetic genome, and the published weights of
Pangolin on an excerpt of GRCh38 chr22 in `tests/pangolin/data/`. They download the published weights once into the
pytest cache, or into the folder in `ABEXP_PANGOLIN_MODELS_DIR` if set. Without network access, the download fails,
and so do the tests with the published weights. Their expected rows come from upstream Pangolin, see
`tests/pangolin/make_expected.py` and `tests/pangolin/make_expected_published.py`.

The absplice2 tests run `absplice2_dna` with a stand-in for the AbSplice2 model on the made-up inputs of
`tests/absplice2/constellations.py`, one constellation per test. Their expected rows come from the AbSplice2 scripts
with the same stand-in, see `tests/absplice2/make_expected.py`.

The expected outputs in `tests/splicing/data/expected/` come from the upstream packages, see
`tests/splicing/make_expected.py`. The SpliceAI-RocksDB tests also need the extra `rocksdb` and the hg38 chr22
database (2.6 GB). Set `ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22` to the path of `spliceAI_hg38_chr22.db`, otherwise pytest
skips them. The test of the SpliceAI predictions also needs SpliceAI.

## License

The code is under the MIT license, see [LICENSE](LICENSE). It lists the copyright notices of the upstream
projects of `abexp.mmsplice`, `abexp.absplice` and `abexp.spliceai_rocksdb`. `abexp.pangolin` loads the weights of Pangolin, which are
GPL-3 and not part of the package, see [abexp.pangolin](#abexppangolin).
