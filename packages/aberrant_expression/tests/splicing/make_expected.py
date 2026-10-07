"""Write the expected test outputs in tests/splicing/data with the upstream packages.

The tests compare aberrant-expression with these files. Run this script only to change the test data. Each
command needs an environment with the upstream packages, e.g. the AbExp environments before abexp-splicing, the
former name of aberrant-expression:

- `python make_expected.py mmsplice`: mmsplice 2.4.0
- `python make_expected.py spliceai_vcf`: the SpliceAI fork hoeze/SpliceAI a1583bd. It runs the SpliceAI command
  line tool on clinvar_chr22.vcf and writes clinvar_chr22.spliceai.vcf.
- `python make_expected.py absplice`: mmsplice 2.4.0, absplice daad7b6 and splicemap cf922eb, with pandas 2.
  absplice daad7b6 fails with pandas 3. Run it after spliceai_vcf.
- `python make_expected.py spliceai_rocksdb <spliceAI_hg38_chr22.db>`: spliceai_rocksdb 3c40d6e and the SpliceAI
  fork.

All commands read the chr22 sequence of the AbExp example, example/chr22_hg38.fa.
"""
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

TESTS = Path(__file__).resolve().parent
DATA = TESTS / 'data'
EXPECTED = DATA / 'expected'
FASTA = TESTS.parents[3] / 'example' / 'chr22_hg38.fa'
VCF = DATA / 'clinvar_chr22.vcf'
SPLICEAI_VCF = DATA / 'clinvar_chr22.spliceai.vcf'
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


def spliceai_vcf():
    subprocess.run(['spliceai', '-I', VCF, '-O', SPLICEAI_VCF, '-R', FASTA, '-A', 'grch38'], check=True)


def absplice():
    import onnxruntime
    from absplice import SpliceOutlierDataloader, SpliceOutlier, SplicingOutlierResult
    from absplice.result import ABSPLICE_DNA, _load_features_from_model_file
    from absplice.utils import read_spliceai_vcf

    # MMSplice with SpliceMaps
    dl = SpliceOutlierDataloader(str(FASTA), str(VCF), splicemap5=[str(SPLICEMAP5)], splicemap3=[str(SPLICEMAP3)])
    SpliceOutlier().predict_save(dl, EXPECTED / 'mmsplice_splicemap.csv')

    # SpliceAI VCF to CSV, as AbExp's spliceai_vcf_to_csv.py did
    df = read_spliceai_vcf(str(SPLICEAI_VCF))
    df = df.rename(columns={'acceptor_loss_positiin': 'acceptor_loss_position'})
    df.to_csv(EXPECTED / 'spliceai_vcf.csv', index=False)

    # AbSplice-DNA; abexp.absplice reads the model inputs with onnxruntime instead of onnx
    session_inputs = [i.name for i in onnxruntime.InferenceSession(ABSPLICE_DNA).get_inputs()]
    assert session_inputs == _load_features_from_model_file(ABSPLICE_DNA), session_inputs
    df_mmsplice = pd.read_csv(EXPECTED / 'mmsplice_splicemap.csv')
    # one entry per variant and gene, see test_predict_absplice_dna
    df_spliceai = pd.read_csv(EXPECTED / 'spliceai_vcf.csv').drop_duplicates(subset=['variant', 'gene_name'])
    result = SplicingOutlierResult(df_mmsplice=df_mmsplice, df_spliceai=df_spliceai)
    result.predict_absplice_dna().reset_index().to_csv(EXPECTED / 'absplice_dna.csv', index=False)


def spliceai_rocksdb(db_path):
    from spliceai_rocksdb.spliceAI import SpliceAI

    model = SpliceAI(str(FASTA), annotation='grch38', db_path={'22': db_path})
    model.predict_save(str(VCF), EXPECTED / 'spliceai_rocksdb.csv', batch_size=1000)


if __name__ == '__main__':
    EXPECTED.mkdir(exist_ok=True)
    {
        'mmsplice': mmsplice,
        'spliceai_vcf': spliceai_vcf,
        'absplice': absplice,
        'spliceai_rocksdb': spliceai_rocksdb,
    }[sys.argv[1]](*sys.argv[2:])
