"""The inputs of the tests in test_absplice2.py, one constellation per test, and the stand-in for the AbSplice2 model.

make_expected.py runs the upstream AbSplice2 scripts on the same inputs. So this module needs only the standard
library, and runs on Python 3.8 like AbSplice2.

The genes and sites are made up. A SpliceMap site of a psi5 event is the junction start on the plus strand and the
junction end on the minus strand; for psi3 it is the other end. Only the minus-strand constellation depends on the
strand, because absplice2_dna does not compare it.
"""
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PangolinRow:
    """A row of abexp.pangolin, with unrounded scores."""
    variant: str
    gene_id: str
    gain_score: float
    gain_pos: int
    loss_score: float
    loss_pos: int


@dataclass(frozen=True)
class SpliceMapRow:
    """A row of the SpliceMap of `event_type` (psi5 or psi3) and `tissue`."""
    tissue: str
    event_type: str
    junctions: str
    gene_id: str
    splice_site: str
    ref_psi: float
    median_n: float


@dataclass(frozen=True)
class MMSpliceRow:
    """A row of MMSplice with SpliceMaps."""
    variant: str
    gene_id: str
    tissue: str
    ref_psi: float
    median_n: float
    delta_logit_psi: float
    delta_psi: float
    junction: str
    event_type: str
    splice_site: str


@dataclass(frozen=True)
class Constellation:
    """The inputs of one test. Each tissue has a psi5 and a psi3 SpliceMap, also without rows."""
    pangolin: tuple
    splicemap: tuple
    mmsplice: tuple = ()
    tissues: tuple = ('Whole_Blood',)


@dataclass(frozen=True)
class InputFiles:
    """The files that write_inputs writes."""
    splicemap5: list
    splicemap3: list
    mmsplice: str


# the new names of the SpliceMap tissues; Brain_Cortex keeps its name
TISSUE_MAPPING = {'Whole_Blood': 'Whole Blood'}

# The stand-in for the AbSplice2 model: AbSplice_DNA = 1 / (1 + exp(-(BIAS + sum of WEIGHTS[i] * feature i))), with
# the features in the order of abexp.absplice2.FEATURES. The weights read the features by position, so the tests pin
# their order. The gain weight is positive and the loss weight negative, so a larger gain or loss raises the score.
WEIGHTS = (0.5, 2.0, 3.0, -3.0, 0.01, 0.02)
BIAS = -4.0

MMSPLICE_COLUMNS = ['variant', 'gene_id', 'tissue', 'ref_psi', 'median_n', 'delta_logit_psi', 'delta_psi',
                    'junction', 'event_type', 'splice_site']


def write_splicemap(constellation, directory, tissue, event_type):
    """Write the SpliceMap of `tissue` and `event_type` into `directory`, and return its path."""
    path = Path(directory) / f'{tissue}_splicemap_{event_type}.csv.gz'
    with gzip.open(path, 'wt', newline='') as fd:
        fd.write(f'# name: {tissue}\n')
        writer = csv.writer(fd)
        # splicemap infers its method from the columns k and n; the scripts read neither
        writer.writerow(['junctions', 'gene_id', 'splice_site', 'ref_psi', 'k', 'n', 'median_n'])
        for r in constellation.splicemap:
            if r.tissue == tissue and r.event_type == event_type:
                writer.writerow([r.junctions, r.gene_id, r.splice_site, r.ref_psi, 1, 2, r.median_n])
    return str(path)


def write_inputs(constellation, directory):
    """Write the SpliceMaps and the MMSplice table of `constellation` into `directory`."""
    directory = Path(directory)
    splicemap5 = [write_splicemap(constellation, directory, tissue, 'psi5') for tissue in constellation.tissues]
    splicemap3 = [write_splicemap(constellation, directory, tissue, 'psi3') for tissue in constellation.tissues]
    mmsplice = directory / 'mmsplice_splicemap.csv'
    with open(mmsplice, 'w', newline='') as fd:
        writer = csv.writer(fd)
        writer.writerow(MMSPLICE_COLUMNS)
        for r in constellation.mmsplice:
            writer.writerow([getattr(r, c) for c in MMSPLICE_COLUMNS])
    return InputFiles(splicemap5=splicemap5, splicemap3=splicemap3, mmsplice=str(mmsplice))


# The gain 10 bp downstream of the variant is 2 bp from a minus-strand psi5 site, and 3 bp from a psi3 site.
GAIN_NEAR_MINUS_STRAND_SITE = Constellation(
    pangolin=(PangolinRow('chr1:1000:A>G', 'ENSG00000000001.4', 0.35, 10, 0.0, -50),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:900-1012:-', 'ENSG00000000001', 'chr1:1012:-', 0.3, 15.0),
        SpliceMapRow('Whole_Blood', 'psi3', 'chr1:1013-1100:-', 'ENSG00000000001', 'chr1:1013:-', 0.8, 40.0),
    ),
)

# The gain and the loss each match a site; the loss has the larger score.
GAIN_AND_LOSS_MATCHED = Constellation(
    pangolin=(PangolinRow('chr1:2000:C>T', 'ENSG00000000002.7', 0.2, 5, -0.45, -15),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:2005-2100:+', 'ENSG00000000002', 'chr1:2005:+', 0.1, 8.0),
        SpliceMapRow('Whole_Blood', 'psi3', 'chr1:1900-1985:+', 'ENSG00000000002', 'chr1:1985:+', 0.9, 30.0),
    ),
    mmsplice=(
        MMSpliceRow('chr1:2000:C>T', 'ENSG00000000002', 'Whole_Blood', 0.9, 30.0, -1.2, -0.25, 'chr1:1900-1985:+',
                    'psi3', 'chr1:1985:+'),
    ),
)

# The gain and the loss each match a site, with scores of the same size.
GAIN_AND_LOSS_EQUAL = Constellation(
    pangolin=(PangolinRow('chr1:2500:G>C', 'ENSG00000000003.1', 0.3, 4, -0.3, -6),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:2504-2600:+', 'ENSG00000000003', 'chr1:2504:+', 0.2, 9.0),
        SpliceMapRow('Whole_Blood', 'psi3', 'chr1:2400-2494:+', 'ENSG00000000003', 'chr1:2494:+', 0.85, 22.0),
    ),
)

# The gain of 0.004 rounds to 0, so its site is the variant, not the psi5 site 7 bp downstream. The loss of 0 also
# sits at the variant. Both match the psi3 site 1 bp upstream.
SCORE_ZERO = Constellation(
    pangolin=(PangolinRow('chr1:3000:G>A', 'ENSG00000000004.2', 0.004, 7, 0.0, -50),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:3007-3100:+', 'ENSG00000000004', 'chr1:3007:+', 0.5, 12.0),
        SpliceMapRow('Whole_Blood', 'psi3', 'chr1:2900-2999:+', 'ENSG00000000004', 'chr1:2999:+', 0.6, 20.0),
    ),
)

# The scores lie halfway between two rounded values: 12.5 and -37.5 hundredths.
ROUNDING_HALF_TO_EVEN = Constellation(
    pangolin=(PangolinRow('chr1:3500:T>G', 'ENSG00000000005.1', 0.125, 1, -0.375, -1),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:3600-3700:+', 'ENSG00000000005', 'chr1:3600:+', 0.4, 11.0),
    ),
)

# A gain of 0.9 counts as 0.7 for the model.
GAIN_ABOVE_CAP = Constellation(
    pangolin=(PangolinRow('chr1:4000:T>C', 'ENSG00000000006.3', 0.9, 3, -0.02, -30),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:4003-4200:+', 'ENSG00000000006', 'chr1:4003:+', 0.05, 25.0),
    ),
    mmsplice=(
        MMSpliceRow('chr1:4000:T>C', 'ENSG00000000006', 'Whole_Blood', 0.05, 25.0, 2.5, 0.4, 'chr1:4003-4200:+',
                    'psi5', 'chr1:4003:+'),
    ),
)

# A deletion of 3 bases; its loss 4 bp downstream of the first base matches a minus-strand psi3 site.
DELETION = Constellation(
    pangolin=(PangolinRow('chr1:5000:ACGT>A', 'ENSG00000000007.5', 0.15, -10, -0.6, 4),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi3', 'chr1:5004-5300:-', 'ENSG00000000007', 'chr1:5004:-', 0.95, 50.0),
    ),
    mmsplice=(
        MMSpliceRow('chr1:5000:ACGT>A', 'ENSG00000000007', 'Whole_Blood', 0.95, 50.0, -3.0, -0.5, 'chr1:5004-5300:-',
                    'psi3', 'chr1:5004:-'),
    ),
)

# Two MMSplice rows with the same features tie for the largest score; they differ in ref_psi and the junction.
TIED_ROWS = Constellation(
    pangolin=(PangolinRow('chr1:6000:A>T', 'ENSG00000000008.1', 0.1, -3, 0.0, -50),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:5997-6100:+', 'ENSG00000000008', 'chr1:5997:+', 0.7, 18.0),
    ),
    mmsplice=(
        MMSpliceRow('chr1:6000:A>T', 'ENSG00000000008', 'Whole_Blood', 0.7, 18.0, 0.8, 0.1, 'chr1:5997-6100:+',
                    'psi5', 'chr1:5997:+'),
        MMSpliceRow('chr1:6000:A>T', 'ENSG00000000008', 'Whole_Blood', 0.4, 18.0, 0.8, 0.1, 'chr1:5900-6100:+',
                    'psi3', 'chr1:6100:+'),
        MMSpliceRow('chr1:6000:A>T', 'ENSG00000000008', 'Whole_Blood', 0.2, 18.0, 0.1, 0.01, 'chr1:5950-6100:+',
                    'psi3', 'chr1:6100:+'),
    ),
)

# A gene in the pseudoautosomal region: its chrX and chrY copies have the same gene_id, and sites at the same
# positions with different values.
CHRX_AND_CHRY_PAR = Constellation(
    pangolin=(
        PangolinRow('chrX:300000:G>C', 'ENSG00000000009.3', 0.3, 4, 0.0, -50),
        PangolinRow('chrY:300000:G>C', 'ENSG00000000009.3', 0.5, 4, 0.0, -50),
    ),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chrX:300004-300500:+', 'ENSG00000000009', 'chrX:300004:+', 0.2, 10.0),
        SpliceMapRow('Whole_Blood', 'psi5', 'chrY:300004-300500:+', 'ENSG00000000009', 'chrY:300004:+', 0.7, 3.0),
    ),
)

# Two tissues: the gain matches a site in Whole_Blood only. Brain_Cortex has no SpliceMap row of the gene, but an
# MMSplice row, and keeps its name.
TWO_TISSUES = Constellation(
    pangolin=(PangolinRow('chr1:7000:C>G', 'ENSG00000000010.2', 0.25, 6, 0.0, -50),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:7006-7200:+', 'ENSG00000000010', 'chr1:7006:+', 0.15, 14.0),
        SpliceMapRow('Brain_Cortex', 'psi5', 'chr1:9000-9100:+', 'ENSG00000000099', 'chr1:9000:+', 0.5, 5.0),
    ),
    mmsplice=(
        MMSpliceRow('chr1:7000:C>G', 'ENSG00000000010', 'Brain_Cortex', 0.6, 7.0, 0.3, 0.05, 'chr1:6900-7100:+',
                    'psi3', 'chr1:7100:+'),
    ),
    tissues=('Whole_Blood', 'Brain_Cortex'),
)

# MMSplice scores the variant in a second gene, which Pangolin did not score.
GENE_ONLY_IN_MMSPLICE = Constellation(
    pangolin=(PangolinRow('chr1:8000:T>A', 'ENSG00000000011.1', 0.05, 2, -0.05, -2),),
    splicemap=(
        SpliceMapRow('Whole_Blood', 'psi5', 'chr1:8500-8600:+', 'ENSG00000000011', 'chr1:8500:+', 0.3, 6.0),
    ),
    mmsplice=(
        MMSpliceRow('chr1:8000:T>A', 'ENSG00000000012', 'Whole_Blood', 0.45, 16.0, 1.5, 0.3, 'chr1:7900-8050:-',
                    'psi5', 'chr1:8050:-'),
    ),
)

CONSTELLATIONS = {
    'gain_near_minus_strand_site': GAIN_NEAR_MINUS_STRAND_SITE,
    'gain_and_loss_matched': GAIN_AND_LOSS_MATCHED,
    'gain_and_loss_equal': GAIN_AND_LOSS_EQUAL,
    'score_zero': SCORE_ZERO,
    'rounding_half_to_even': ROUNDING_HALF_TO_EVEN,
    'gain_above_cap': GAIN_ABOVE_CAP,
    'deletion': DELETION,
    'tied_rows': TIED_ROWS,
    'chrx_and_chry_par': CHRX_AND_CHRY_PAR,
    'two_tissues': TWO_TISSUES,
    'gene_only_in_mmsplice': GENE_ONLY_IN_MMSPLICE,
}
