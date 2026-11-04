import numpy as np
import pytest

import synthetic as s
from abexp.promoterai.dataloader import VCFTSSDataloader


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


def records(dl):
    """The metadata of the samples of a dataloader, in order."""
    return [sample['metadata'] for sample in dl]


def test_variant_in_one_promoter(genome, vcf_file):
    dl = make_dataloader(genome, vcf_file(s.SNV_PLUS))
    assert len(dl) == 1
    metadata, = records(dl)
    assert metadata == {
        'seq_start': 295 - s.SEQ_LENGTH // 2,  # 0-based start of the variant minus half the sequence
        'seq_end': 295 - s.SEQ_LENGTH // 2 + s.SEQ_LENGTH,
        'tss': 299,  # 0-based, the first base of PLUS1-T1
        'chrom': 'chr1',
        'strand': '+',
        'gene_id': 'PLUS1',
        'transcript_id': 'PLUS1-T1',
        'transcript_start': 299,
        'transcript_end': 600,
        'variant_start': 295,  # 0-based
        'variant_end': 296,  # 1-based
        'ref': s.SNV_PLUS[2],
        'alt': s.SNV_PLUS[3],
    }


def test_variant_on_the_minus_strand(genome, vcf_file):
    dl = make_dataloader(genome, vcf_file(s.SNV_MINUS))
    metadata, = records(dl)
    # the TSS of a transcript on the minus strand is its last base, 1-based 1000
    assert (metadata['strand'], metadata['tss'], metadata['gene_id']) == ('-', 999, 'MINUS1')


def test_variant_in_two_promoters_of_one_gene(genome, vcf_file):
    dl = make_dataloader(genome, vcf_file(s.SNV_BOTH_TRANSCRIPTS))
    assert len(dl) == 2
    assert sorted(m['transcript_id'] for m in records(dl)) == ['PLUS1-T1', 'PLUS1-T2']


def test_variant_outside_of_the_promoters(genome, vcf_file):
    dl = make_dataloader(genome, vcf_file(s.FAR_FROM_TSS))
    assert len(dl) == 0
    assert records(dl) == []


def test_canonical_only(genome, vcf_file):
    dl = make_dataloader(genome, vcf_file(s.SNV_BOTH_TRANSCRIPTS), canonical_only=True)
    assert [m['transcript_id'] for m in records(dl)] == ['PLUS1-T1']


def test_protein_coding_only(genome, vcf_file):
    path = vcf_file(s.SNV_LNCRNA)
    assert [m['gene_id'] for m in records(make_dataloader(genome, path))] == ['LNC1']
    assert records(make_dataloader(genome, path, protein_coding_only=True)) == []


def test_gene_ids(genome, vcf_file):
    path = vcf_file(s.SNV_PLUS, s.SNV_MINUS)
    assert [m['gene_id'] for m in records(make_dataloader(genome, path, gene_ids=['MINUS1']))] == ['MINUS1']


def test_size(genome, vcf_file):
    path = vcf_file(s.SNV_BOTH_TRANSCRIPTS)
    dl = make_dataloader(genome, path, size=1)
    assert len(dl) == 1
    assert len(records(dl)) == 1


@pytest.mark.parametrize('variant', [s.UNKNOWN_ALLELE, s.REF_MISMATCH, s.CHROM_START])
def test_variant_that_cannot_be_scored(genome, vcf_file, variant):
    dl = make_dataloader(genome, vcf_file(variant))
    sample, = list(dl)
    # the record is there, but it holds no sequences
    assert not sample['scorable']
    assert sample['metadata']['variant_start'] == variant[1] - 1
    for name in ('ref', 'alt'):
        np.testing.assert_array_equal(sample['sequences'][name], np.zeros((s.SEQ_LENGTH, 4), dtype=np.float32))


def test_batch_iter(genome, vcf_file):
    dl = make_dataloader(genome, vcf_file(s.SNV_PLUS, s.SNV_BOTH_TRANSCRIPTS, s.SNV_MINUS))
    assert len(dl) == 4
    batches = list(dl.batch_iter(batch_size=3))
    assert [len(batch['scorable']) for batch in batches] == [3, 1]
    for batch in batches:
        assert batch['sequences']['ref'].shape == (len(batch['scorable']), s.SEQ_LENGTH, 4)
        assert batch['sequences']['alt'].shape == (len(batch['scorable']), s.SEQ_LENGTH, 4)
        assert len(batch['metadata']['gene_id']) == len(batch['scorable'])


def test_pyarrow_metadata_schema(genome, vcf_file):
    dl = make_dataloader(genome, vcf_file(s.SNV_PLUS))
    schema = dl.pyarrow_metadata_schema
    assert schema.names == ['seq_start', 'seq_end', 'tss', 'chrom', 'strand', 'gene_id', 'transcript_id',
                            'transcript_start', 'transcript_end', 'variant_start', 'variant_end', 'ref', 'alt']
    assert schema.metadata == {
        b'seq_length': str(s.SEQ_LENGTH).encode(),
        b'variant_upstream_tss': str(s.UPSTREAM_TSS).encode(),
        b'variant_downstream_tss': str(s.DOWNSTREAM_TSS).encode(),
    }
