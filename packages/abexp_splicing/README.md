# abexp-splicing

tl;dr: The parts of MMSplice, AbSplice and SpliceAI-RocksDB that the AbExp pipeline runs, in one package on
kipoiseq2. kipoi, kipoiseq 0.7, cyvcf2 and pyranges are no longer needed. The outputs match the upstream packages,
except for the row order, rounding in the last digits and the choice among tied junctions, see
[Differences](#differences-from-the-upstream-packages).

| import                    | ported from                                                                                     | what AbExp uses                                                     |
| ------------------------- | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------- |
| `abexp.mmsplice`          | [mmsplice](https://github.com/gagneurlab/MMSplice_MTSplice) 2.4.0 (31513da)                     | `MMSplice` and the junction VCF dataloaders                         |
| `abexp.absplice`          | [absplice](https://github.com/gagneurlab/absplice) daad7b6, [splicemap](https://github.com/gagneurlab/splicemap) cf922eb | `SpliceOutlierDataloader`, `SpliceOutlier`, `SplicingOutlierResult`, `read_spliceai_vcf` |
| `abexp.spliceai_rocksdb`  | [spliceai_rocksdb](https://github.com/gagneurlab/spliceai_rocksdb) 3c40d6e                      | `SpliceAI`                                                          |

## Installation

Install a release from the AbExp repository with pip:

<!-- x-release-please-start-version -->
```bash
pip install "abexp-splicing[rocksdb] @ git+https://github.com/gagneurlab/AbExp.git@abexp-splicing-v0.0.1#subdirectory=packages/abexp_splicing"
```
<!-- x-release-please-end -->

The extra `rocksdb` installs python-rocksdb and SpliceAI for `abexp.spliceai_rocksdb`. The AbExp pipeline uses
the SpliceAI fork [hoeze/SpliceAI](https://github.com/hoeze/SpliceAI), which also predicts multi-nucleotide
variants. pip keeps an installed SpliceAI, so install the fork first or in the same pip command.

For development, install the package in editable mode from a checkout of AbExp:

```bash
pip install -e "packages/abexp_splicing[rocksdb,dev]"
```

## Usage

The example runs AbSplice-DNA on the test data of this package. Run it in the root of a checkout of AbExp.

```python
import pandas as pd
from abexp.absplice import SpliceOutlier, SpliceOutlierDataloader, SplicingOutlierResult

data_dir = 'packages/abexp_splicing/tests/data'

# MMSplice with SpliceMaps: delta PSI per variant, junction and tissue
dl = SpliceOutlierDataloader(
    'example/chr22_hg38.fa', f'{data_dir}/clinvar_chr22.vcf',
    splicemap5=[f'{data_dir}/Whole_Blood_splicemap_psi5.csv.gz'],
    splicemap3=[f'{data_dir}/Whole_Blood_splicemap_psi3.csv.gz'],
)
SpliceOutlier().predict_save(dl, 'mmsplice_splicemap.csv')

# AbSplice-DNA per variant, gene and tissue, from MMSplice and SpliceAI
result = SplicingOutlierResult(
    df_mmsplice=pd.read_csv('mmsplice_splicemap.csv'),
    df_spliceai=pd.read_csv(f'{data_dir}/expected/spliceai_vcf.csv'),
)
df = result.predict_absplice_dna()
```

`abexp.spliceai_rocksdb.SpliceAI(fasta, annotation='grch38', db_path={'22': <path>}).predict_save(vcf, csv)`
looks up the SpliceAI scores in SpliceAI-RocksDB and runs SpliceAI for the variants that are not in it.
`abexp.absplice.read_spliceai_vcf` reads a VCF file that the SpliceAI command line tool annotated.

## Differences from the upstream packages

- kipoiseq2 yields the variant-junction pairs in the order of the VCF file, and kipoiseq 0.7 in the order of
  pyranges. So the rows of the MMSplice table come in another order.
- MMSplice therefore predicts batches with other samples. TensorFlow rounds them differently, so `delta_logit_psi`
  can differ in the last digits.
- A variant can have the same delta PSI at several junctions of a gene. absplice daad7b6 picked any one of them,
  and the row order and the numpy version decided which one. abexp-splicing picks the first in the order of
  junction, event_type, splice_site and the other columns. So the reported junction can differ from absplice, and
  with its `median_n` also `splice_site_is_expressed` and `AbSplice_DNA`.
- `read_spliceai_vcf` names the column `acceptor_loss_position`, not `acceptor_loss_positiin`.
- `SpliceOutlier` works with pandas 3. absplice concatenated the psi5 and psi3 rows of each batch. With pandas 3,
  an empty part turned the PSI columns into objects, and the delta PSI failed. abexp-splicing skips the empty part.
- A VCF file without variants yields no samples. mmsplice compared the contigs of the VCF header with the FASTA
  file, and abexp-splicing compares the chromosomes of the variants.
- Only the parts that AbExp uses are kept. Not included: MTSplice, the VEP plugin, the dataloaders of GTF exons
  and exon tables, the VCF writers, the variant filters of absplice, AbSplice-RNA, the CADD-Splice and
  sample-based features, and the SpliceMap count tables. `predict_absplice_dna` reads only ONNX models.

## Tests

The tests read the files in `tests/data/` and the hg38 chr22 sequence in `example/chr22_hg38.fa` of the AbExp
repository. Pytest skips the tests that need the sequence, if it is missing, e.g. in an sdist.

```bash
pytest packages/abexp_splicing/tests
```

The expected outputs in `tests/data/expected/` come from the upstream packages, see `tests/make_expected.py`.
The SpliceAI-RocksDB test also needs the extra `rocksdb` and the hg38 chr22 database (2.6 GB). Set
`ABEXP_SPLICEAI_ROCKSDB_HG38_CHR22` to the path of `spliceAI_hg38_chr22.db`, otherwise pytest skips it.

## License

The code is under the MIT license, see [LICENSE](LICENSE). It lists the copyright notices of the upstream
projects.
