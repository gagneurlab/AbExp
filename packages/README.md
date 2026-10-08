# Python packages

tl;dr: The workflow installs these packages from this repository, pinned to their release tags. release-please
releases them; kipoi-enformer and abexp-splicing always release together, at one version. A change of a package
reaches the workflow only with its next release.

| directory        | distribution     | import                                                       | tag                         |
| ---------------- | ---------------- | ------------------------------------------------------------ | --------------------------- |
| `abexp_utils`    | `abexp-utils`    | `abexp_utils`                                                | `abexp-utils-v<version>`    |
| `kipoi_enformer` | `kipoi-enformer` | `kipoi_enformer`                                             | `kipoi-enformer-v<version>` |
| `abexp_splicing` | `abexp-splicing` | `abexp.mmsplice`, `abexp.absplice`, `abexp.spliceai_rocksdb` | `abexp-splicing-v<version>` |

## How the workflow installs the packages

The conda environments install the packages with pip from a git URL of this repository, pinned to a release tag,
e.g. `...@kipoi-enformer-v0.0.1#subdirectory=packages/kipoi_enformer`:
- `abexp-utils` in `workflow/envs/abexp-veff-py.yaml`
- `kipoi-enformer` and `abexp-splicing` in `workflow/modules/veff/envs/abexp-tensorflow-cpu.yaml` and
  `abexp-tensorflow-cuda.yaml`

A relative path to `packages/` does not work. Snakemake copies the environment files to `.snakemake/conda/`
before conda runs pip, so the path would resolve from there. Also, a workflow that imports AbExp as a module
has no copy of `packages/`.

Between releases, main installs the packages of the last release.

## Running the tests

The tests run across all cores by default, through pytest-xdist. Pass `-n0` to run them in one process, which a
debugger needs and which restores per-test output order.

```bash
pip install -e "./packages/abexp_utils[all,dev]"
cd packages/abexp_utils && pytest
```
The spark tests need Java.

The kipoi-enformer tests read example files from `packages/kipoi_enformer/tests/data/` and the chr22 sequence
from `example/chr22_hg19.fa`. git-lfs stores the files in `tests/data/`, and a clone fetches them only on request.
Install git-lfs and fetch them once:
```bash
git lfs pull --include="packages/kipoi_enformer/tests/data/**" --exclude=""
pip install -e "./packages/kipoi_enformer[dev]"
cd packages/kipoi_enformer && pytest
```
Without the files, pytest skips the tests that need them.

The fetch is opt-in, because `.lfsconfig` excludes `tests/data/` from the default fetch. The conda environments
pip-install the packages from a git URL of this repository. Without the exclude, each environment build on a machine
with git-lfs would download the test data and use up the LFS quota of the repository. Once the quota runs out, the
download fails, and so does the pip install. `--exclude=""` in the command above lifts the exclude.

The abexp-splicing tests read the files in `packages/abexp_splicing/tests/data/`, which git stores directly, and
the chr22 sequence from `example/chr22_hg38.fa`:
```bash
pip install -e "./packages/abexp_splicing[dev]"
cd packages/abexp_splicing && pytest
```
The SpliceAI-RocksDB test also needs the extra `rocksdb` and the hg38 chr22 database. Set
`ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22` to the path of `spliceAI_hg38_chr22.db`, otherwise pytest skips it. pip
builds python-rocksdb only with the RocksDB library installed; conda-forge has a build of it.

CI (`.github/workflows/packages.yml`) runs the tests on pull requests and pushes to main that change `packages/`.
It fetches the test data with the same command and caches it.

A push of new test data needs git-lfs in the pushing clone, so that its pre-push hook uploads the LFS objects.

## Testing a package change in the workflow

Install the package in editable mode into the conda environment that Snakemake created:
```bash
snakemake --list-conda-envs  # shows the location of each environment
<location>/bin/pip install -e ./packages/kipoi_enformer
```
Snakemake does not notice the change. Force the jobs you want to test with `--forcerun <rule>`.
Snakemake builds a new environment when the environment file changes, e.g. with the next release. The new
environment has the released package again.

## Releases

`.github/workflows/release-please.yml` runs release-please on each push to main. It keeps one release PR open for
abexp-utils, and one for kipoi-enformer and abexp-splicing together, see [Linked versions](#linked-versions). A
third release PR covers the workflow itself, see "Releases" in the [README](../README.md#releases). A
commit counts for a package if it changes files in the package directory and if its type appears in the changelog:
`feat`, `fix`, `perf`, `deps`, `revert` or `docs`. Before 1.0.0, a breaking change bumps the minor version. The
tags `abexp-utils-v0.0.1` and `kipoi-enformer-v0.0.1` mark the state imported from gagneurlab/AbExp-utils, and
`abexp-splicing-v0.0.1` marks the first state of abexp-splicing.

The release PR of a package bumps:
- the version in its `pyproject.toml` and its `CHANGELOG.md`
- the tag pins in the conda environments and in the package README. These are the lines marked
  `x-release-please-version` and the lines between `x-release-please-start-version` and `x-release-please-end`.

Merging the release PR creates the tag, e.g. `kipoi-enformer-v0.1.0`, and the GitHub release. The pins in the
conda environments then point to the new tag, and Snakemake rebuilds the affected environments once. No jobs
rerun, because `software-env` is not a rerun trigger, see [When jobs rerun](../README.md#when-jobs-rerun). If
the new version changes the outputs of the workflow, bump the affected `OUTPUT_VERSION` keys in a follow-up commit.

The release PR moves the pins but does not change the workflow scripts. So a breaking change in a package must
keep the workflow working with both the old and the new release: keep the old name as a deprecated alias, and
update the workflow after the pin has moved.

After a release, the workflow builds the sdist and wheel of the package. PyPI rejects direct URL dependencies.
If the package has one, the build emits a warning. If the upload to PyPI is on, the build fails instead.
kipoi-enformer and abexp-splicing depend on kipoiseq2 0.1 from PyPI, so pip can install them only after that
release.

### Linked versions

kipoi-enformer and abexp-splicing always release together, at the same version. Both pins in the TensorFlow
environment files carry the marker `x-release-please-version`, and the marker names no package. So the release PR
of one package would also move the pin of the other one to its version. The plugin `linked-versions` in
`release-please-config.json` (group `abexp-tensorflow`) avoids that:
- Both packages get one release PR, `chore: release abexp-tensorflow libraries`.
- Both get the higher of their two new versions, and both pins move to it. A package without changes gets a release
  whose changelog only says that it synchronizes the versions.
- `group-pull-request-title-pattern` lets release-please recognize the merged PR of the group. Without it,
  release-please creates no tags and no GitHub releases for the group.

The TensorFlow environments pin kipoi-enformer 0.0.1, which needs a kipoiseq2 commit from before 0.1, and
abexp-splicing, which needs kipoiseq2 0.1. conda can build these environments only after the first release of the
group has moved both pins.

Notes:
- GitHub starts no workflows for pull requests that the `GITHUB_TOKEN` opens, so the tests do not run on release PRs.
- release-please needs the repository setting "Allow GitHub Actions to create and approve pull requests".

## Enabling PyPI

The upload to PyPI is prepared but off. To turn it on:
1. Check that kipoiseq2 0.1 is on PyPI.
2. On PyPI, add a trusted publisher to each project (`abexp-utils`, `kipoi-enformer`, `abexp-splicing`): owner
   `gagneurlab`, repository `AbExp`, workflow `release-please.yml`, environment `pypi`. For a new project, add a
   pending publisher.
3. Create the GitHub environment `pypi` in the repository settings.
4. Set the repository variable `PYPI_PUBLISH` to `true`.
5. Switch the pins in the conda environments to versions, e.g. `kipoi-enformer==X.Y.Z  # x-release-please-version`,
   and the install snippets in the package READMEs.
