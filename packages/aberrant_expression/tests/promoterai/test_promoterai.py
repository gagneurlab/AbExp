import numpy as np
import polars as pl
import pytest

import synthetic as s
from abexp.promoterai.dataloader import VCFTSSDataloader
from abexp.promoterai.promoterai import PromoterAI, PromoterAIVeff, human_output
from abexp.promoterai.utils import RandomModel

# the scores of the test below and of PromoterAI differ in the last float32 digits
ABS_TOL = 1e-6


def make_dataloader(genome, vcf_path, **kwargs):
    args = {
        'fasta_file': genome / 'genome.fa',
        'genome_annotation': genome / 'genes.gff3',
        'vcf_file': vcf_path,
        'seq_length': s.SEQ_LENGTH,
        'variant_upstream_tss': s.UPSTREAM_TSS,
        'variant_downstream_tss': s.DOWNSTREAM_TSS,
    }
    return VCFTSSDataloader(**(args | kwargs))


def expected_scores(dl):
    """
    The scores of RandomModel for the samples of a dataloader: the tanh of the rounded mean difference of the
    outputs. The weights of RandomModel depend only on its seed, so this model is the one of PromoterAI.
    """
    weights = RandomModel().track_weights.numpy()
    scores = []
    for sample in dl:
        if not sample['scorable']:
            scores.append(float('nan'))
            continue
        difference = sample['sequences']['alt'] @ weights - sample['sequences']['ref'] @ weights
        scores.append(float(np.tanh(np.round(difference.mean(), 4))))
    return scores


@pytest.fixture
def predict(genome, vcf_file, tmp_path):
    """Scores the given VCF records with RandomModel; returns the dataloader, the scores and their file."""
    counter = 0

    def predict(*records, batch_size=2, dataloader_kwargs=None, **model_kwargs):
        nonlocal counter
        counter += 1
        dl = make_dataloader(genome, vcf_file(*records), **(dataloader_kwargs or {}))
        path = tmp_path / f'scores_{counter}.parquet'
        PromoterAI(is_random=True, **model_kwargs).predict(dl, batch_size=batch_size, filepath=path)
        return dl, pl.read_parquet(path), path
    return predict


def veff(annotation, score_path, output_path, aggregation_mode='canonical', **kwargs):
    PromoterAIVeff(genome_annotation=annotation). \
        run(score_path, output_path, aggregation_mode=aggregation_mode, **kwargs)
    return pl.read_parquet(output_path)


def test_model_path_is_required():
    with pytest.raises(ValueError, match='model_path is required'):
        PromoterAI()


@pytest.mark.parametrize('outputs, expected', [
    ({'human': 1, 'nonhuman': 2}, 1),
    # a legacy model returns a tuple or a list with the human head first
    ((1, 2), 1),
    ([1, 2], 1),
    # the outputs of a restored SavedModel are named, e.g. output0 for the human head
    ({'output0': 1, 'output1': 2}, 1),
    (1, 1),
])
def test_human_output(outputs, expected):
    assert human_output(outputs) == expected


@pytest.mark.parametrize('batch_size', [1, 2, 3, 10])
def test_predict(predict, batch_size):
    dl, df, _ = predict(s.SNV_PLUS, s.SNV_BOTH_TRANSCRIPTS, s.SNV_MINUS, batch_size=batch_size)
    assert df.columns == ['score', *dl.pyarrow_metadata_schema.names]
    assert df.schema['score'] == pl.Float32
    assert df.height == 4
    assert df['score'].to_list() == pytest.approx(expected_scores(dl), abs=ABS_TOL)
    # PromoterAI scores lie between -1 and 1
    assert df['score'].abs().max() <= 1


@pytest.mark.parametrize('sequences_per_call', [1, 2])
def test_predict_splits_a_large_batch_into_several_calls(predict, monkeypatch, sequences_per_call):
    # the bound lets `sequences_per_call` sequences into one call of the model, so a batch of 3 needs several calls
    monkeypatch.setattr('abexp.promoterai.promoterai.MAX_POSITIONS_PER_CALL', sequences_per_call * s.SEQ_LENGTH + 1)
    calls = []
    call_model = PromoterAI._call_model
    monkeypatch.setattr(PromoterAI, '_call_model', lambda self, x: calls.append(len(x)) or call_model(self, x))
    dl, df, _ = predict(s.SNV_PLUS, s.SNV_BOTH_TRANSCRIPTS, s.SNV_MINUS, batch_size=3)
    assert max(calls) == sequences_per_call
    assert df['score'].to_list() == pytest.approx(expected_scores(dl), abs=ABS_TOL)


@pytest.mark.parametrize('output_type', ['dict', 'tuple', 'tensor'])
def test_predict_reads_the_human_head(predict, output_type):
    _, df, _ = predict(s.SNV_PLUS, output_type=output_type)
    _, expected, _ = predict(s.SNV_PLUS, output_type='dict')
    assert df['score'].to_list() == pytest.approx(expected['score'].to_list(), abs=ABS_TOL)


def test_predict_scores_the_strand_of_the_transcript(predict):
    # the same variant in the promoters of PLUS3 and of MINUS2
    _, df, _ = predict(s.SNV_BOTH_STRANDS)
    assert sorted(df['gene_id'].to_list()) == ['MINUS2', 'PLUS3']
    # the sequences of a transcript on the minus strand are reverse complemented, so its score differs
    assert df['score'][0] != df['score'][1]


def test_predict_scores_the_transcripts_of_a_gene_alike(predict):
    # the two transcripts of PLUS1 lie on the same strand, so the variant has the same score for both
    _, df, _ = predict(s.SNV_BOTH_TRANSCRIPTS)
    assert df.height == 2
    assert df['score'][0] == df['score'][1]


@pytest.mark.parametrize('variant', [s.UNKNOWN_ALLELE, s.REF_MISMATCH, s.CHROM_START])
def test_predict_writes_nan_for_a_variant_that_cannot_be_scored(predict, variant):
    _, df, _ = predict(variant)
    assert df.height == 1
    assert np.isnan(df['score'][0])


def test_predict_on_an_empty_dataloader(predict):
    dl, df, _ = predict(s.FAR_FROM_TSS)
    assert len(dl) == 0
    assert df.height == 0
    assert df.columns == ['score', *dl.pyarrow_metadata_schema.names]


@pytest.mark.parametrize('aggregation_mode', ['canonical', 'median'])
def test_veff(genome, predict, tmp_path, aggregation_mode):
    _, scores, score_path = predict(s.SNV_BOTH_TRANSCRIPTS)
    df = veff(genome / 'genes.gff3', score_path, tmp_path / 'veff.parquet', aggregation_mode)
    # one row per gene and variant; both transcripts of PLUS1 have the same score
    assert df.schema == pl.Schema(PromoterAIVeff.OUTPUT_SCHEMA)
    assert df.height == 1
    assert df['gene_id'].to_list() == ['PLUS1']
    assert df['veff_score'][0] == pytest.approx(scores['score'][0], abs=ABS_TOL)


def test_veff_sorts_the_rows_by_gene_and_variant(genome, predict, tmp_path):
    _, _, score_path = predict(s.SNV_MINUS, s.SNV_PLUS, s.SNV_LNCRNA)
    df = veff(genome / 'genes.gff3', score_path, tmp_path / 'veff.parquet', 'median')
    assert df.select('chrom', 'strand', 'gene_id', 'variant_start').rows() == [
        ('chr1', '+', 'LNC1', s.SNV_LNCRNA[1] - 1),
        ('chr1', '+', 'PLUS1', s.SNV_PLUS[1] - 1),
        ('chr1', '-', 'MINUS1', s.SNV_MINUS[1] - 1),
    ]


def test_veff_keeps_only_the_canonical_transcripts(genome, predict, tmp_path):
    # the lncRNA gene LNC1 has no canonical protein coding transcript
    _, _, score_path = predict(s.SNV_LNCRNA, s.SNV_PLUS)
    canonical = veff(genome / 'genes.gff3', score_path, tmp_path / 'canonical.parquet', 'canonical')
    median = veff(genome / 'genes.gff3', score_path, tmp_path / 'median.parquet', 'median')
    assert canonical['gene_id'].to_list() == ['PLUS1']
    assert median['gene_id'].to_list() == ['LNC1', 'PLUS1']


def test_veff_fails_with_several_canonical_transcripts(genome, predict, tmp_path):
    _, _, score_path = predict(s.SNV_BOTH_TRANSCRIPTS)
    # a second canonical transcript of PLUS1
    transcripts = tuple(
        (*t[:7], 'basic,Ensembl_canonical') if t[2] == 'PLUS1-T2' else t for t in s.TRANSCRIPTS
    )
    annotation = tmp_path / 'two_canonical.gff3'
    s.write_gff3(annotation, transcripts)
    with pytest.raises(ValueError, match='Multiple canonical transcripts'):
        veff(annotation, score_path, tmp_path / 'veff.parquet', 'canonical')


@pytest.mark.parametrize('kwargs, expected_rows', [
    # the variant lies 4 bases upstream of the TSS of PLUS1-T1 and 19 upstream of the one of PLUS1-T2
    ({}, 1),
    ({'upstream_tss': 5}, 1),  # only PLUS1-T1
    ({'upstream_tss': 2}, 0),
    ({'downstream_tss': 0}, 1),
])
def test_veff_filters_the_distance_to_the_tss(genome, predict, tmp_path, kwargs, expected_rows):
    _, scores, score_path = predict(s.SNV_PLUS, dataloader_kwargs={'variant_upstream_tss': 20})
    assert sorted(scores['transcript_id'].to_list()) == ['PLUS1-T1', 'PLUS1-T2']
    df = veff(genome / 'genes.gff3', score_path, tmp_path / 'veff.parquet', 'median', **kwargs)
    # the two transcripts of PLUS1 have the same score, so the kept ones give one row
    assert df.height == expected_rows


def test_veff_drops_the_variants_that_cannot_be_scored(genome, predict, tmp_path):
    _, _, score_path = predict(s.UNKNOWN_ALLELE, s.SNV_PLUS)
    df = veff(genome / 'genes.gff3', score_path, tmp_path / 'veff.parquet', 'median')
    assert df['variant_start'].to_list() == [s.SNV_PLUS[1] - 1]


def test_veff_on_empty_scores(genome, predict, tmp_path):
    _, _, score_path = predict(s.FAR_FROM_TSS)
    df = veff(genome / 'genes.gff3', score_path, tmp_path / 'veff.parquet', 'median')
    assert df.height == 0
    assert df.schema == pl.Schema(PromoterAIVeff.OUTPUT_SCHEMA)


def test_veff_strips_the_id_versions(genome, predict, tmp_path):
    # GENCODE gives the genes and the transcripts a version, which the output drops
    transcripts = tuple((t[0], t[1] + '.2', t[2] + '.3', *t[3:]) for t in s.TRANSCRIPTS)
    annotation = tmp_path / 'versioned.gff3'
    s.write_gff3(annotation, transcripts)
    _, scores, score_path = predict(s.SNV_PLUS, dataloader_kwargs={'genome_annotation': annotation})
    assert scores['gene_id'].to_list() == ['PLUS1.2']
    df = veff(annotation, score_path, tmp_path / 'veff.parquet', 'canonical')
    assert df['gene_id'].to_list() == ['PLUS1']
