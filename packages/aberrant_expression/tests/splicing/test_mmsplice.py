import json

import numpy as np
import polars as pl
import pytest

from abexp.mmsplice import JunctionPSI3VCFDataloader, JunctionPSI5VCFDataloader, MMSplice, batch_iter, encodeDNA, \
    numpy_collate
from abexp.mmsplice.dataloader import SeqSpliter
from conftest import EXPECTED_DIR, SPLICEMAP3, SPLICEMAP5, VCF, read_junctions

# the sequences of make_expected.py
MODULAR_SEQS = [
    ('ATGCGACGTACCCAGTAAAT', (4, 4)),
    ('CTTTCTTTCCCCACAGGTTCAGCTGCAGGTGAGCATCAGGAGGTAAGTTCTGGACTTTGGATAGCTGACTAGCTCATTGTCTAGAGG', (16, 20)),
]


def test_encodeDNA():
    np.testing.assert_array_equal(
        encodeDNA(['ACGTN']),
        np.array([[[1., 0., 0., 0.],
                   [0., 1., 0., 0.],
                   [0., 0., 1., 0.],
                   [0., 0., 0., 1.],
                   [0., 0., 0., 0.]]])
    )

    # padded with N at the end to the longest sequence
    encoded = encodeDNA(['AA', 'ATT', 'ATTCGG'])
    assert encoded.shape == (3, 6, 4)
    np.testing.assert_array_equal(encoded[0, :2], [[1., 0., 0., 0.], [1., 0., 0., 0.]])
    np.testing.assert_array_equal(encoded[0, 2:], np.zeros((4, 4)))
    np.testing.assert_array_equal(encoded[1, 3:], np.zeros((3, 4)))
    np.testing.assert_array_equal(encoded[2, 5], [0., 0., 1., 0.])


def test_numpy_collate():
    samples = [
        {'inputs': {'seq': 'AC'}, 'metadata': {'pos': 1, 'score': 0.5, 'region': None, 'start': np.int64(3)}},
        {'inputs': {'seq': 'GTT'}, 'metadata': {'pos': 2, 'score': 1.5, 'region': None, 'start': np.int64(4)}},
    ]
    batch = numpy_collate(samples)
    np.testing.assert_array_equal(batch['inputs']['seq'], np.array(['AC', 'GTT']))
    np.testing.assert_array_equal(batch['metadata']['pos'], np.array([1, 2]))
    np.testing.assert_array_equal(batch['metadata']['score'], np.array([0.5, 1.5]))
    np.testing.assert_array_equal(batch['metadata']['start'], np.array([3, 4]))
    assert list(batch['metadata']['region']) == [None, None]

    batches = list(batch_iter(iter(samples * 3), batch_size=4))
    assert [len(b['metadata']['pos']) for b in batches] == [4, 2]


def test_mmsplice_modular_scores():
    model = MMSplice()
    spliter = SeqSpliter()
    scores = []
    for seq, overhang in MODULAR_SEQS:
        batch = {k: encodeDNA([v]) for k, v in spliter.split(seq, overhang).items()}
        scores.append(model.predict_modular_scores_on_batch(batch)[0])
    expected = pl.read_csv(EXPECTED_DIR / 'mmsplice_modular_scores.csv')
    np.testing.assert_allclose(np.array(scores), expected.to_numpy(), rtol=0, atol=1e-5)



def test_read_junction():
    # from test_JunctionVCFDataloader_read_junction of mmsplice
    junctions = pl.DataFrame({
        'Chromosome': ['17', '17'],
        'Start': [41276032, 41276002],
        'End': [41279742, 41279042],
        'Strand': ['-', '+'],
    })

    df = JunctionPSI3VCFDataloader._read_junction(junctions, 'psi3', overhang=(100, 100), exon_len=50)
    assert df.height == 2
    assert df.select('Start', 'End').rows() == [(41279742 - 100, 41279742 + 50), (41276002 - 50, 41276002 + 100)]
    assert df['junction'].to_list() == ['17:41276032-41279742:-', '17:41276002-41279042:+']

    df = JunctionPSI5VCFDataloader._read_junction(junctions, 'psi5', overhang=(100, 100), exon_len=50)
    assert df.select('Start', 'End').rows() == [(41276032 - 50, 41276032 + 100), (41279042 - 100, 41279042 + 50)]


@pytest.mark.parametrize('event_type, cls, splicemap', [
    ('psi5', JunctionPSI5VCFDataloader, SPLICEMAP5),
    ('psi3', JunctionPSI3VCFDataloader, SPLICEMAP3),
])
def test_junction_vcf_dataloader(fasta_file, event_type, cls, splicemap):
    # the sequences and metadata of each variant-junction pair, as in mmsplice 2.4.0
    samples = []
    for row in cls(read_junctions(splicemap), fasta_file, str(VCF), encode=False):
        meta = row['metadata']
        samples.append({
            'event_type': event_type,
            'variant': meta['variant']['annotation'],
            'region': meta['variant']['region'],
            'junction': meta['exon']['junction'],
            'exon': meta['exon']['annotation'],
            'left_overhang': int(meta['exon']['left_overhang']),
            'right_overhang': int(meta['exon']['right_overhang']),
            'seq': row['inputs']['seq'],
            'mut_seq': row['inputs']['mut_seq'],
        })
    with open(EXPECTED_DIR / 'junction_samples.jsonl') as f:
        expected = [s for s in map(json.loads, f) if s['event_type'] == event_type]

    def key(s):
        return s['variant'], s['junction']

    # kipoiseq2 yields the pairs in the order of the VCF file, kipoiseq in the order of pyranges
    assert len(samples) == len(expected)
    assert sorted(samples, key=key) == sorted(expected, key=key)
