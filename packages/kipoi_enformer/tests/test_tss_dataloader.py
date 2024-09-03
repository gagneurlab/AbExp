import pytest

from kipoi_enformer.dataloader import VCFTSSDataloader, RefTSSDataloader
from kipoi_enformer.dataloader.dataloader import get_tss_from_genome_annotation
from kipoiseq.transforms.functional import one_hot2string

UPSTREAM_TSS = 10
DOWNSTREAM_TSS = 10


@pytest.fixture
def variants():
    return {
        # Positive strand
        # SNP
        # ===|TSS|===|Var|=======================
        'chr22:16364873:G>A:ENST00000438441.1': {
            'chrom': 'chr22',
            'strand': '+',
            'tss': 16364866,  # 0-based
            'ref_start': 16364856,  # 0-based
            'ref_stop': 16364877,  # 1-based
            'var_start': 16364872,  # 0-based
            'var_stop': 16364873,  # 1-based
            'ref': 'G',
            'alt': 'A',
            'ref_seq': 'ACTGGCTGGCCATGCCGTCCC',
            'alt_seq': 'ACTGGCTGGCCATGCCATCCC',
        },
        # Positive strand
        # SNP
        # ===|Var|===|TSS|=======================
        'chr22:17565895:G>C:ENST00000694950.1_1': {
            'chrom': 'chr22',
            'strand': '+',
            'tss': 17565901,  # 0-based
            'ref_start': 17565891,  # 0-based
            'ref_stop': 17565912,  # 1-based
            'var_start': 17565894,  # 0-based
            'var_stop': 17565895,  # 1-based
            'ref': 'G',
            'alt': 'C',
            'ref_seq': 'CTCGAACTCCACCGCGGAAAA',
            'alt_seq': 'CTCCAACTCCACCGCGGAAAA',
        },
        # Negative strand
        # SNP
        # ===|Var|===|TSS|=======================
        'chr22:16570002:C>T:ENST00000583607.1': {
            'chrom': 'chr22',
            'strand': '-',
            'tss': 16570005,  # 0-based
            'ref_start': 16569995,  # 0-based
            'ref_stop': 16570016,  # 1-based
            'var_start': 16570001,  # 0-based
            'var_stop': 16570002,  # 1-based
            'ref': 'C',
            'alt': 'T',
            # complement: CTGCAACGAGGGTCTGCATGT
            'ref_seq': 'ACATGCAGACCCTCGTTGCAG',
            # complement: CTGCAATGAGGGTCTGCATGT
            'alt_seq': 'ACATGCAGACCCTCATTGCAG',
        },
        # Negative strand
        # Deletion
        # ===|Var|===|TSS|=======================
        'chr22:29130718:CAAA>C:ENST00000403642.5_3': {
            'chrom': 'chr22',
            'strand': '-',
            'tss': 29130708,  # 0-based
            'ref_start': 29130698,  # 0-based
            'ref_end': 29130719,  # 1-based
            'var_start': 29130717,  # 0-based
            'var_end': 29130721,  # 1-based
            'ref': 'CAAA',
            'alt': 'C',
            # complement: TCCCGAGACATCACGACCTCA
            'ref_seq': 'TGAGGTCGTGATGTCTCGGGA',
            # complement: TCCCGAGACATCACGACCTCA
            'alt_seq': 'TGAGGTCGTGATGTCTCGGGA',
        },
        # Negative strand
        # Insertion
        # ===|Var|===|TSS|=======================
        'chr22:19109971:T>TCCCGCCC:ENST00000545799.5_4': {
            'chrom': 'chr22',
            'strand': '-',
            'tss': 19109966,  # 0-based
            'ref_start': 19109956,  # 0-based
            'ref_end': 19109977,  # 1-based
            'var_start': 19109970,  # 0-based
            'var_end': 19109971,  # 1-based
            'ref': 'T',
            'alt': 'TCCCGCCC',
            # complement: CGCCCCGCCCCGCCTCCCGCC
            'ref_seq': 'GGCGGGAGGCGGGGCGGGGCG',
            # complement: CGCCCCGCCCCGCCTCCCGCC
            'alt_seq': 'GGCGGGAGGCGGGGCGGGGCG',
        },
        # special case: the TSS is within the variant's interval; we take the downstream TSS
        # Negative strand
        # Deletion
        # ===|Var|===|TSS|=======================
        'chr22:18359465:GTTATGGAGGTTAGGGAGGTTATGGAGGTTAGGGAGC>G:ENST00000462645.1_3': {
            'chrom': 'chr22',
            'strand': '-',
            'tss': 18359468,
            'ref_start': 18359458,
            'ref_end': 18359479,
            'var_start': 18359464,
            'var-stop': 18359501,
            'ref': 'GTTATGGAGGTTAGGGAGGTTATGGAGGTTAGGGAGC',
            'alt': 'G',
            # complement 'TGCAGGGTTATGGAGGTTAGG',
            'ref_seq': 'CCTAACCTCCATAACCCTGCA',
            # complement 'ACA TGCAGGG TTATGGAGGTT',
            'alt_seq': 'AACCTCCATAACCCTGCATGT',
        },
    }


@pytest.fixture()
def references():
    return {
        # Positive strand
        # SNP
        # ===|TSS|===|Var|=======================
        'chr22:ENST00000438441.1': {
            'chrom': 'chr22',
            'strand': '+',
            'tss': 16364866,  # 0-based
            'start': 16364856,  # 0-based
            'end': 16364877,  # 1-based
            'seq': 'ACTGGCTGGCCATGCCGTCCC',
        },
        # Positive strand
        # SNP
        # ===|Var|===|TSS|=======================
        'chr22:ENST00000694950.1_1': {
            'chrom': 'chr22',
            'strand': '+',
            'tss': 17565901,  # 0-based
            'start': 17565891,  # 0-based
            'end': 17565912,  # 1-based
            'seq': 'CTCGAACTCCACCGCGGAAAA',
        },
        # Negative strand
        # SNP
        # ===|Var|===|TSS|=======================
        'chr22:ENST00000583607.1': {
            'chrom': 'chr22',
            'strand': '-',
            'tss': 16570005,  # 0-based
            'start': 16569995,  # 0-based
            'end': 16570016,  # 1-based
            # complement: CTGCAACGAGGGTCTGCATGT
            'seq': 'ACATGCAGACCCTCGTTGCAG',
        },
        # Negative strand
        # Deletion
        # ===|Var|===|TSS|=======================
        'chr22:ENST00000403642.5_3': {
            'chrom': 'chr22',
            'strand': '-',
            'tss': 29130708,  # 0-based
            'start': 29130698,  # 0-based
            'end': 29130719,  # 1-based
            # complement: TCCCGAGACATCACGACCTCA
            'seq': 'TGAGGTCGTGATGTCTCGGGA',
        },
        # Negative strand
        # Insertion
        # ===|Var|===|TSS|=======================
        'chr22:ENST00000545799.5_4': {
            'chrom': 'chr22',
            'strand': '-',
            'tss': 19109966,  # 0-based
            'start': 19109956,  # 0-based
            'end': 19109977,  # 1-based
            # complement: CGCCCCGCCCCGCCTCCCGCC
            'seq': 'GGCGGGAGGCGGGGCGGGGCG',
        },
        # special case: the TSS is within the variant's interval; we take the downstream TSS
        # Negative strand
        # Deletion
        # ===|Var|===|TSS|=======================
        'chr22:ENST00000462645.1_3': {
            'chrom': 'chr22',
            'strand': '-',
            'tss': 18359468,
            'start': 18359458,
            'end': 18359479,
            # complement 'TGCAGGGTTATGGAGGTTAGG',
            'seq': 'CCTAACCTCCATAACCCTGCA',
        },
    }


def test_get_tss_from_genome_annotation(chr22_example_files):
    roi = get_tss_from_genome_annotation(chr22_example_files['gtf'], chromosome='chr22', protein_coding_only=False,
                                         canonical_only=False)

    # check the number of transcripts
    # grep -v '^#' annot.chr22.gtf | cut -f 3 | grep transcript | wc -l
    assert len(roi) == 5279

    # are the extracted ROIs correct
    # criteria are the strand and transcript start and end
    # roi start is zero based
    # roi end is 1 based
    def test_tss(row):
        # make tss zero-based
        if row.Strand == '-':
            tss = row.transcript_end - 1
        else:
            tss = row.transcript_start

        assert row.Start == tss, \
            f'Transcript {row.transcript_id}; Strand {row.Strand}: {row.Start} != {tss}'
        assert row.End == tss + 1, \
            f'Transcript {row.transcript_id}; Strand {row.Strand}: {row.End} != ({tss} + 1)'

    roi.apply(test_tss, axis=1)

    # check the extracted ROI for a negative strand transcript
    roi_i = roi.set_index('transcript_id').loc['ENST00000448070.1']
    assert roi_i.Start == (16076172 - 1)
    assert roi_i.End == 16076172

    # check the extracted ROI for a positive strand transcript
    roi_i = roi.set_index('transcript_id').loc['ENST00000424770.1']
    assert roi_i.Start == (16062157 - 1)
    assert roi_i.End == 16062157


def test_genome_annotation_protein_canonical(chr22_example_files):
    # Ground truth bash:
    # grep -E '^\S+\s+\S+\s+transcript' annot.chr22.gtf | grep -E 'tag\s+"Ensembl_canonical"' |
    # grep -E 'transcript_type\s+"protein_coding"' | wc -l
    chromosome = 'chr22'
    # not all genes have a canonical transcript if the genome annotation is not current (e.g. GRCh37 is not current)
    roi = get_tss_from_genome_annotation(chr22_example_files['gtf'], chromosome=chromosome, protein_coding_only=True,
                                         canonical_only=True)
    assert len(roi) == 441
    roi = get_tss_from_genome_annotation(chr22_example_files['gtf'], chromosome=chromosome, protein_coding_only=True,
                                         canonical_only=True, gene_ids=['ENSG00000172967'])
    assert roi.gene_id.iloc[0] == 'ENSG00000172967.8_5'
    roi = get_tss_from_genome_annotation(chr22_example_files['gtf'], chromosome=chromosome, protein_coding_only=False,
                                         canonical_only=False)
    assert len(roi) == 5279
    roi = get_tss_from_genome_annotation(chr22_example_files['gtf'], chromosome=chromosome, protein_coding_only=True,
                                         canonical_only=False)
    assert len(roi) == 3623
    roi = get_tss_from_genome_annotation(chr22_example_files['gtf'], chromosome=chromosome, protein_coding_only=False,
                                         canonical_only=True)
    assert len(roi) == 1212


def test_vcf_dataloader(chr22_example_files, variants):
    dl = VCFTSSDataloader(
        fasta_file=chr22_example_files['fasta'],
        gtf=chr22_example_files['gtf'],
        vcf_file=chr22_example_files['vcf'],
        variant_downstream_tss=10,
        variant_upstream_tss=10,
        seq_length=21,
        shifts=[0],
    )
    total = 0
    checked_variants = dict()
    for i in dl:
        total += 1
        metadata = i['metadata']
        # example: chr22:16364873:G>A_
        var_id = (f'{metadata["chrom"]}:{metadata["variant_start"] + 1}:'
                  f'{metadata["ref"]}>{metadata["alt"]}:'
                  f'{metadata["transcript_id"]}')
        variant = variants.get(var_id, None)
        if variant is not None:
            # assert one_hot2string(i['sequences']['shift:0'][None, :, :])[0] == variant['ref_seq']
            assert one_hot2string(i['sequences'])[0] == variant['alt_seq']
            checked_variants[var_id] = 2

    # check that all variants in my list were found and checked
    assert set(checked_variants.keys()) == set(variants.keys())
    print(total)


def test_ref_dataloader(chr22_example_files, references):
    dl = RefTSSDataloader(
        fasta_file=chr22_example_files['fasta'],
        gtf=chr22_example_files['gtf'],
        seq_length=21,
        shifts=[0],
        chromosome='chr22',
    )
    total = 0
    checked_refs = dict()
    for i in dl:
        total += 1
        metadata = i['metadata']
        # example: chr22:16364873:G>A_
        ref_id = f'chr22:{metadata["transcript_id"]}'
        ref = references.get(ref_id)
        if ref is not None:
            assert one_hot2string(i['sequences'])[0] == ref['seq']
            checked_refs[ref_id] = 2

    # check that all variants in my list were found and checked
    assert set(checked_refs.keys()) == set(references.keys())
    print(total)
