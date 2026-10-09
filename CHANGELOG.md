# Changelog

## [2.0.0](https://github.com/gagneurlab/AbExp/compare/v1.1.0...v2.0.0) (2026-10-09)


### ⚠ BREAKING CHANGES

* **vep:** support VEP 116 and drop VEP 99 and 105 ([#30](https://github.com/gagneurlab/AbExp/issues/30))
* merge the AbExp packages into aberrant-expression ([#18](https://github.com/gagneurlab/AbExp/issues/18))
* store the AbExp models as LightGBM text models ([#14](https://github.com/gagneurlab/AbExp/issues/14))
* remove the key `features` from custom entries of `system.models`, because config validation now rejects it.
* run MMSplice, AbSplice and SpliceAI-RocksDB on kipoiseq2 ([#8](https://github.com/gagneurlab/AbExp/issues/8))
* download resources in the modules, rerun on output versions ([#10](https://github.com/gagneurlab/AbExp/issues/10))
* move to Snakemake 9 and split the pipeline into modules ([#9](https://github.com/gagneurlab/AbExp/issues/9))
* Rename `gtf_file` to `gff3_file` and point it at the GENCODE GFF3 of the same release.

### Features

* **absplice2:** score Pangolin with abexp.pangolin ([#19](https://github.com/gagneurlab/AbExp/issues/19)) ([31313e6](https://github.com/gagneurlab/AbExp/commit/31313e68886f227f6d9ab70047d755c5e7498410))
* add release-please and options to turn off CADD and LOFTEE ([#13](https://github.com/gagneurlab/AbExp/issues/13)) ([e44a326](https://github.com/gagneurlab/AbExp/commit/e44a3264d3c2ba1082fa79e1e8f83c682bde2806))
* download resources in the modules, rerun on output versions ([#10](https://github.com/gagneurlab/AbExp/issues/10)) ([b0f8a03](https://github.com/gagneurlab/AbExp/commit/b0f8a031dbfa11c972e07384b61da9a7783e0a1c))
* **nmd_scanner:** run NMD-Scanner 0.5.0 on the GFF3 and build NMD features per tissue ([#16](https://github.com/gagneurlab/AbExp/issues/16)) ([f261606](https://github.com/gagneurlab/AbExp/commit/f261606557d80c525111e19b5e7a5940ac512254))
* read the genome annotation from a GENCODE GFF3 ([d093251](https://github.com/gagneurlab/AbExp/commit/d09325102d459d5b17c0afe97aaf9d478628511b))
* run MMSplice, AbSplice and SpliceAI-RocksDB on kipoiseq2 ([#8](https://github.com/gagneurlab/AbExp/issues/8)) ([e7f7a1c](https://github.com/gagneurlab/AbExp/commit/e7f7a1c7b7984779d1b8f8e95f56409a318c3e3f))
* store the AbExp models as LightGBM text models ([#14](https://github.com/gagneurlab/AbExp/issues/14)) ([3f64cb2](https://github.com/gagneurlab/AbExp/commit/3f64cb298214bf854820490affcf3bd598ee3e22))
* **veff:** add AbSplice2 as an optional splicing module ([#12](https://github.com/gagneurlab/AbExp/issues/12)) ([7967294](https://github.com/gagneurlab/AbExp/commit/7967294047d8e90228b9a5ee1068928325c1ba23))
* **veff:** add mehari as transcript consequence annotator ([d093251](https://github.com/gagneurlab/AbExp/commit/d09325102d459d5b17c0afe97aaf9d478628511b))
* **veff:** add NMD-Scanner as an NMD escape module ([d093251](https://github.com/gagneurlab/AbExp/commit/d09325102d459d5b17c0afe97aaf9d478628511b))
* **veff:** add reloftee as a loss-of-function module ([d093251](https://github.com/gagneurlab/AbExp/commit/d09325102d459d5b17c0afe97aaf9d478628511b))
* **vep:** support VEP 116 and drop VEP 99 and 105 ([#30](https://github.com/gagneurlab/AbExp/issues/30)) ([9f405f7](https://github.com/gagneurlab/AbExp/commit/9f405f78ea84b0cb5fb73c5bb20358c29a97eab0))


### Bug Fixes

* **deps:** require polars-bio 0.36.1 ([#36](https://github.com/gagneurlab/AbExp/issues/36)) ([8fe3fba](https://github.com/gagneurlab/AbExp/commit/8fe3fba3b42356f1a3514407e4ada0a7c7ead02b))
* feed abexp_v1.1_Enformer its features in training order ([3f64cb2](https://github.com/gagneurlab/AbExp/commit/3f64cb298214bf854820490affcf3bd598ee3e22))
* **tissue_specific_vep:** count each PAR transcript once ([b5bb939](https://github.com/gagneurlab/AbExp/commit/b5bb9391d48044faaed10ff93dabd760d9a1b65b))
* **tissue_specific_vep:** keep variants that end at a batch boundary ([f261606](https://github.com/gagneurlab/AbExp/commit/f261606557d80c525111e19b5e7a5940ac512254))


### Code Refactoring

* merge the AbExp packages into aberrant-expression ([#18](https://github.com/gagneurlab/AbExp/issues/18)) ([0235b79](https://github.com/gagneurlab/AbExp/commit/0235b79969df41982e9a6f74f30ed2da6fc64493))
* move to Snakemake 9 and split the pipeline into modules ([#9](https://github.com/gagneurlab/AbExp/issues/9)) ([d093251](https://github.com/gagneurlab/AbExp/commit/d09325102d459d5b17c0afe97aaf9d478628511b))
