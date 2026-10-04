# abexp-splicing

tl;dr: The parts of MMSplice that the AbExp pipeline runs, in one package on kipoiseq2. kipoi, kipoiseq 0.7 and
pyranges are no longer needed. The outputs match mmsplice, except for the row order and rounding in the last
digits, see [Differences](#differences-from-the-upstream-packages).

| import                    | ported from                                                                                     | what AbExp uses                                                     |
| ------------------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| `abexp.mmsplice`          | [mmsplice](https://github.com/gagneurlab/MMSplice_MTSplice) 2.4.0 (31513da)                     | `MMSplice` and the junction VCF dataloaders                         |

## Installation

Install a release from the AbExp repository with pip:

<!-- x-release-please-start-version -->
```bash
pip install "abexp-splicing @ git+https://github.com/gagneurlab/AbExp.git@abexp-splicing-v0.0.1#subdirectory=packages/abexp_splicing"
```
<!-- x-release-please-end -->

For development, install the package in editable mode from a checkout of AbExp:

```bash
pip install -e "packages/abexp_splicing[dev]"
```

## Differences from the upstream packages

- kipoiseq2 yields the variant-junction pairs in the order of the VCF file, and kipoiseq 0.7 in the order of
  pyranges. So the samples of the junction dataloaders come in another order.
- MMSplice therefore predicts batches with other samples. TensorFlow rounds them differently, so `delta_logit_psi`
  can differ in the last digits.
- A VCF file without variants yields no samples. mmsplice compared the contigs of the VCF header with the FASTA
  file, and abexp-splicing compares the chromosomes of the variants.
- Only the parts that AbExp uses are kept. Not included: MTSplice, the VEP plugin, the dataloaders of GTF exons
  and exon tables, and the VCF writers.

## Tests

The tests read the files in `tests/data/` and the hg38 chr22 sequence in `example/chr22_hg38.fa` of the AbExp
repository. Pytest skips the tests that need the sequence, if it is missing, e.g. in an sdist.

```bash
pytest packages/abexp_splicing/tests
```

The expected outputs in `tests/data/expected/` come from the upstream packages, see `tests/make_expected.py`.

## License

The code is under the MIT license, see [LICENSE](LICENSE). It lists the copyright notices of the upstream
projects.
