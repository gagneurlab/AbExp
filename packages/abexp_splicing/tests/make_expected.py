"""Write the expected test outputs in tests/data with the upstream packages.

The tests compare abexp-splicing with these files. Run this script only to change the test data. Each command
needs an environment with the upstream packages, e.g. the AbExp environments before abexp-splicing:

- `python make_expected.py mmsplice`: mmsplice 2.4.0

All commands read the chr22 sequence of the AbExp example, example/chr22_hg38.fa.
"""
import json
import sys
from pathlib import Path

import pandas as pd

TESTS = Path(__file__).resolve().parent
DATA = TESTS / 'data'
EXPECTED = DATA / 'expected'
FASTA = TESTS.parents[2] / 'example' / 'chr22_hg38.fa'
VCF = DATA / 'clinvar_chr22.vcf'
SPLICEMAP5 = DATA / 'Whole_Blood_splicemap_psi5.csv.gz'
SPLICEMAP3 = DATA / 'Whole_Blood_splicemap_psi3.csv.gz'
# sequences with overhang for the MMSplice modules; the first one is from the mmsplice tests
MODULAR_SEQS = [
    ('ATGCGACGTACCCAGTAAAT', (4, 4)),
    ('CTTTCTTTCCCCACAGGTTCAGCTGCAGGTGAGCATCAGGAGGTAAGTTCTGGACTTTGGATAGCTGACTAGCTCATTGTCTAGAGG', (16, 20)),
]


def read_junctions(path):
    """The junctions of a SpliceMap file, as absplice gives them to the junction dataloaders of mmsplice."""
    df = pd.read_csv(path, skiprows=1)
    return df[['junctions', 'Chromosome', 'Start', 'End', 'Strand']] \
        .drop_duplicates(subset='junctions').set_index('junctions')


def mmsplice():
    import numpy as np
    from mmsplice import MMSplice
    from mmsplice.junction_dataloader import JunctionPSI5VCFDataloader, JunctionPSI3VCFDataloader

    model = MMSplice()
    scores = [model.predict_on_seq(seq, overhang) for seq, overhang in MODULAR_SEQS]
    pd.DataFrame(np.array(scores), columns=['acceptor_intron', 'acceptor', 'exon', 'donor', 'donor_intron']) \
        .to_csv(EXPECTED / 'mmsplice_modular_scores.csv', index=False)

    # samples of the junction dataloaders
    with open(EXPECTED / 'junction_samples.jsonl', 'w') as f:
        for event_type, cls, path in [('psi5', JunctionPSI5VCFDataloader, SPLICEMAP5),
                                      ('psi3', JunctionPSI3VCFDataloader, SPLICEMAP3)]:
            for row in cls(read_junctions(path), str(FASTA), str(VCF), encode=False):
                meta = row['metadata']
                f.write(json.dumps({
                    'event_type': event_type,
                    'variant': meta['variant']['annotation'],
                    'region': meta['variant']['region'],
                    'junction': meta['exon']['junction'],
                    'exon': meta['exon']['annotation'],
                    'left_overhang': int(meta['exon']['left_overhang']),
                    'right_overhang': int(meta['exon']['right_overhang']),
                    'seq': row['inputs']['seq'],
                    'mut_seq': row['inputs']['mut_seq'],
                }) + '\n')


if __name__ == '__main__':
    EXPECTED.mkdir(exist_ok=True)
    {
        'mmsplice': mmsplice,
    }[sys.argv[1]](*sys.argv[2:])
