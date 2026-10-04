# AbExp-utils

This package contains utility functions for the [AbExp](https://github.com/gagneurlab/abexp) project. 
It contains functions to manipulates polars and spark dataframes, as well as scikit-learn models used by AbExp. 

## Installation

Install a release from the AbExp repository with pip.
Choose one of the extras `polars`, `spark` or `all`:

<!-- x-release-please-start-version -->
```bash
pip install "abexp-utils[polars] @ git+https://github.com/gagneurlab/AbExp.git@abexp-utils-v0.0.1#subdirectory=packages/abexp_utils"
```
<!-- x-release-please-end -->

For development, install the package in editable mode from a checkout of AbExp:

```bash
pip install -e "packages/abexp_utils[dev,all]"
```

## Tests

```bash
pytest packages/abexp_utils/tests
```

The tests run across all cores by default, through pytest-xdist. Pass `-n0` to run them in one process, which a
debugger needs and which restores per-test output order.

The Spark tests need a Java runtime. pyspark 4 requires Java 17 or newer.
