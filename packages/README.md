# Python package

tl;dr: `aberrant_expression/` holds the Python distribution aberrant-expression, which the workflow imports as
`abexp`. The workflow installs it from this repository, pinned to its release tag. release-please releases it. A
change of the package reaches the workflow only with its next release.

| directory             | distribution          | import                                                                                                          | tag                              |
| --------------------- | --------------------- | --------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| `aberrant_expression` | `aberrant-expression` | `abexp.utils`, `abexp.enformer`, `abexp.mmsplice`, `abexp.absplice`, `abexp.spliceai_rocksdb`, `abexp.pangolin` | `aberrant-expression-v<version>` |

The base of aberrant-expression has no dependencies. Each extra adds the requirements of some subpackages:

| extra      | subpackages                        | requirements, and the extras it refers to                                                     |
| ---------- | ---------------------------------- | --------------------------------------------------------------------------------------------- |
| `models`   | `abexp.utils.models`               | scikit-learn, `numpy`                                                                         |
| `polars`   | `abexp.utils.polars_functions`     | polars                                                                                        |
| `spark`    | `abexp.utils.spark_functions`      | pyspark, pandas                                                                               |
| `gff3`     | `abexp.utils.gff3`                 | polars-bio, `polars`                                                                          |
| `splicing` | `abexp.mmsplice`, `abexp.absplice` | onnxruntime, `polars`, `numpy`, `kipoiseq2`, `tensorflow`, `tqdm`                             |
| `rocksdb`  | `abexp.spliceai_rocksdb`           | python-rocksdb, `splicing`                                                                    |
| `enformer` | `abexp.enformer`                   | kagglehub, pyarrow, xarray, zarr, pyyaml, `gff3`, `models`, `kipoiseq2`, `tensorflow`, `tqdm` |
| `pangolin` | `abexp.pangolin`                   | torch, `gff3`, `numpy`, `kipoiseq2`                                                           |
| `all`      | all but `abexp.spliceai_rocksdb`   | all extras but `rocksdb`                                                                      |

`pyproject.toml` names each third-party requirement once, with its version range. An extra that needs a requirement
refers to the extra that holds it, e.g. `pangolin` lists `aberrant-expression[gff3]` for polars-bio. The extras
`numpy`, `kipoiseq2`, `tensorflow` and `tqdm` hold one requirement each, which several extras share.

aberrant-expression replaces three distributions: abexp-utils (`abexp_utils`, now `abexp.utils`), kipoi-enformer
(`kipoi_enformer`, now `abexp.enformer`) and abexp-splicing (`abexp.mmsplice`, `abexp.absplice` and
`abexp.spliceai_rocksdb`). The modules keep their paths below the new top level, e.g. `kipoi_enformer.enformer` is
now `abexp.enformer.enformer`. The name `abexp` on PyPI belongs to an unrelated library, which also installs into
`site-packages/abexp`. So the two cannot share an environment.

## How the workflow installs the package

The conda environments install aberrant-expression with pip from a git URL of this repository, pinned to a release
tag, e.g. `...@aberrant-expression-v0.0.1#subdirectory=packages/aberrant_expression`:
- the extra `polars` in `workflow/envs/abexp-veff-py.yaml`
- the extra `pangolin` in `workflow/modules/veff/absplice2/envs/absplice2-pangolin-cpu.yaml` and
  `absplice2-pangolin-cuda.yaml`
- the extras `enformer`, `splicing` and `rocksdb` in `workflow/modules/veff/envs/abexp-tensorflow-cpu.yaml` and
  `abexp-tensorflow-cuda.yaml`

A relative path to `packages/` does not work. Snakemake copies the environment files to `.snakemake/conda/`
before conda runs pip, so the path would resolve from there. Also, a workflow that imports AbExp as a module
has no copy of `packages/`.

Between releases, main installs the package of the last release.

## Running the tests

The tests are in one directory per area: `tests/utils`, `tests/enformer`, `tests/splicing` and `tests/pangolin`.
Each area needs the extras of its subpackages: `polars`, `spark` and `gff3` for `tests/utils`, and the extra of the
same name for the other areas. The tests run across all cores by default, through pytest-xdist. Pass `-n0` to run
them in one process, which a debugger needs and which restores per-test output order.

```bash
pip install -e "./packages/aberrant_expression[all,dev]"
cd packages/aberrant_expression && pytest tests/utils
```
The spark tests need Java.

The enformer tests read example files from `packages/aberrant_expression/tests/enformer/data/` and the chr22
sequence from `example/chr22_hg19.fa`. git-lfs stores the files in `tests/enformer/data/`, and a clone fetches them
only on request. Install git-lfs and fetch them once:
```bash
git lfs pull --include="packages/aberrant_expression/tests/enformer/data/**" --exclude=""
cd packages/aberrant_expression && pytest tests/enformer
```
Without the files, pytest skips the tests that need them.

The fetch is opt-in, because `.lfsconfig` excludes `tests/enformer/data/` from the default fetch. The conda
environments pip-install the package from a git URL of this repository. Without the exclude, each environment build
on a machine with git-lfs would download the test data and use up the LFS quota of the repository. Once the quota
runs out, the download fails, and so does the pip install. `--exclude=""` in the command above lifts the exclude.

The splicing tests read the files in `packages/aberrant_expression/tests/splicing/data/`, which git stores directly,
and the chr22 sequence from `example/chr22_hg38.fa`:
```bash
cd packages/aberrant_expression && pytest tests/splicing
```
The SpliceAI-RocksDB test also needs the extra `rocksdb` and the hg38 chr22 database. Set
`ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22` to the path of `spliceAI_hg38_chr22.db`, otherwise pytest skips it. pip
builds python-rocksdb only with the RocksDB library installed; conda-forge has a build of it.

The pangolin tests run small networks with fixed random weights on a synthetic genome, and the published weights of
Pangolin on an excerpt of GRCh38 chr22 in `packages/aberrant_expression/tests/pangolin/data/`, which git stores
directly. They download the published weights (35 MB) once into the pytest cache, or into the folder in
`ABEXP_PANGOLIN_MODELS_DIR` if set. Without network access, the tests with the published weights fail:
```bash
cd packages/aberrant_expression && pytest tests/pangolin
```

CI (`.github/workflows/packages.yml`) runs the tests of each area on pull requests and pushes to main that change
`packages/`. It installs only the extras of the area. For the enformer tests, it fetches the test data with the same
command and caches it. For the pangolin tests, it caches the published weights of Pangolin, with a key from
`abexp/pangolin/weights.py`, which holds their SHA-256 sums. A second job installs each extra alone and imports all
modules of its subpackages, so that a missing requirement fails. It skips the extra `rocksdb`, because pip cannot build
python-rocksdb there.

A push of new test data needs git-lfs in the pushing clone, so that its pre-push hook uploads the LFS objects.

## Testing a package change in the workflow

Install the package in editable mode into the conda environment that Snakemake created:
```bash
snakemake --list-conda-envs  # shows the location of each environment
<location>/bin/pip install -e ./packages/aberrant_expression
```
Snakemake does not notice the change. Force the jobs you want to test with `--forcerun <rule>`.
Snakemake builds a new environment when the environment file changes, e.g. with the next release. The new
environment has the released package again.

## Releases

`.github/workflows/release-please.yml` runs release-please on each push to main. It keeps one release PR open for
aberrant-expression. A second release PR covers the workflow itself, see "Releases" in the
[README](../README.md#releases). A commit counts for the package if it changes files in `packages/aberrant_expression`
and if its type appears in the changelog: `feat`, `fix`, `perf`, `deps`, `revert` or `docs`. Before 1.0.0, a
breaking change bumps the minor version.

The tag `aberrant-expression-v0.0.1` marks the first state of aberrant-expression, which the conda environments pin
until the first release. The tag has no GitHub release, so release-please finds no release of aberrant-expression.
The first release PR therefore counts all commits since the `bootstrap-sha` in `release-please-config.json` and
proposes 0.1.0.

The release PR bumps:
- the version in `pyproject.toml` and the `CHANGELOG.md` of the package
- the tag pins in the conda environments and in the package README. These are the lines marked
  `x-release-please-version` and the lines between `x-release-please-start-version` and `x-release-please-end`.

Merging the release PR creates the tag, e.g. `aberrant-expression-v0.1.0`, and the GitHub release. The pins in the
conda environments then point to the new tag, and Snakemake rebuilds the affected environments once. No jobs
rerun, because `software-env` is not a rerun trigger, see [When jobs rerun](../README.md#when-jobs-rerun). If
the new version changes the outputs of the workflow, bump the affected `OUTPUT_VERSION` keys in a follow-up commit.

The release PR moves the pins but does not change the workflow scripts. So a breaking change in the package must
keep the workflow working with both the old and the new release: keep the old name as a deprecated alias, and
update the workflow after the pin has moved.

After a release, the workflow builds the sdist and wheel of the package. PyPI rejects direct URL dependencies.
If the package has one, the build emits a warning. If the upload to PyPI is on, the build fails instead.
The extras `splicing`, `enformer` and `pangolin` depend on kipoiseq2 0.1 from PyPI. kipoiseq2 requires Python
3.12 or later, so aberrant-expression does too.

Notes:
- GitHub starts no workflows for pull requests that the `GITHUB_TOKEN` opens, so the tests do not run on release PRs.
- release-please needs the repository setting "Allow GitHub Actions to create and approve pull requests".

## Bioconda

Planned: after the first release, the conda environments install aberrant-expression from bioconda instead of the
git URL, e.g. `bioconda::aberrant-expression-polars=X.Y.Z` in `abexp-veff-py.yaml`, and
`aberrant-expression-enformer` and `aberrant-expression-rocksdb` in the TensorFlow environments. pip then installs
only SpliceAI there. Bioconda builds a version only after its tag exists, so the release PR cannot bump these pins.
Bump them by hand once bioconda has published the version.

The recipe would build one noarch package from the archive of the release tag on GitHub:
- `aberrant-expression-base` installs all code, with the base dependencies.
- The metapackages `aberrant-expression-models`, `-polars`, `-spark`, `-gff3`, `-splicing`, `-rocksdb`,
  `-enformer` and `-pangolin` add the dependencies of each extra, and pin `-base` to the same build. Like the
  extras, `-rocksdb` depends on `-splicing`, and `-enformer` and `-pangolin` depend on `-gff3`.
- `aberrant-expression` depends on all of them.

Before that, kipoiseq2 needs a bioconda recipe.

## Enabling PyPI

The upload to PyPI is prepared but off. To turn it on:
1. On PyPI, add a trusted publisher to the project `aberrant-expression`: owner `gagneurlab`, repository `AbExp`,
   workflow `release-please.yml`, environment `pypi`. For a new project, add a pending publisher.
2. Create the GitHub environment `pypi` in the repository settings.
3. Set the repository variable `PYPI_PUBLISH` to `true`.
4. Switch the install snippet in the package README to a version, e.g.
   `pip install "aberrant-expression[polars]==X.Y.Z"`.
