import polars as pl
import pytest

from abexp.enformer.constants import AlleleType
from abexp.enformer.dataloader import TSSDataloader, VCFTSSDataloader, RefTSSDataloader
from abexp.enformer.dataloader.dataloader import get_tss_from_genome_annotation
from abexp.enformer.enformer import EnformerVeff
from abexp.enformer.utils import read_gff3
from kipoiseq2.transforms.functional import one_hot2string

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
    roi = get_tss_from_genome_annotation(chr22_example_files['genome_annotation'], chromosome='chr22',
                                         protein_coding_only=False, canonical_only=False)

    # check the number of transcripts
    # zcat chr22.gencode.v40lift37.annotation.gff3.gz | cut -f 3 | grep -cx transcript
    assert len(roi) == 5279

    # are the extracted ROIs correct
    # criteria are the strand and transcript start and end
    # roi start is zero based
    # roi end is 1 based
    for row in roi.iter_rows(named=True):
        # make tss zero-based
        if row['Strand'] == '-':
            tss = row['transcript_end'] - 1
        else:
            tss = row['transcript_start']

        assert row['Start'] == tss, \
            f"Transcript {row['transcript_id']}; Strand {row['Strand']}: {row['Start']} != {tss}"
        assert row['End'] == tss + 1, \
            f"Transcript {row['transcript_id']}; Strand {row['Strand']}: {row['End']} != ({tss} + 1)"

    # check the extracted ROI for a negative strand transcript
    roi_i = roi.row(by_predicate=pl.col('transcript_id') == 'ENST00000448070.1', named=True)
    assert roi_i['Start'] == (16076172 - 1)
    assert roi_i['End'] == 16076172

    # check the extracted ROI for a positive strand transcript
    roi_i = roi.row(by_predicate=pl.col('transcript_id') == 'ENST00000424770.1', named=True)
    assert roi_i['Start'] == (16062157 - 1)
    assert roi_i['End'] == 16062157


def test_genome_annotation_from_pandas(chr22_example_files):
    # a pandas DataFrame with pyranges-style columns gives the same TSS as the GFF3 file
    annotation = read_gff3(chr22_example_files['genome_annotation'])
    annotation_pandas = annotation.to_pandas()
    # pyranges stores the chromosome as a category
    annotation_pandas['Chromosome'] = annotation_pandas['Chromosome'].astype('category')

    columns = ['Chromosome', 'Start', 'End', 'Strand', 'gene_id', 'transcript_id', 'tag', 'tss',
               'transcript_start', 'transcript_end']
    args = dict(chromosome='chr22', protein_coding_only=True, canonical_only=True)
    roi = get_tss_from_genome_annotation(annotation, **args).select(columns)
    roi_from_pandas = get_tss_from_genome_annotation(annotation_pandas, **args).select(columns)
    assert len(roi) == 441
    assert roi_from_pandas.equals(roi)


def test_genome_annotation_protein_canonical(chr22_example_files):
    # Ground truth bash:
    # zcat chr22.gencode.v40lift37.annotation.gff3.gz | awk -F'\t' '$3 == "transcript"' |
    # grep 'tag=[^;]*Ensembl_canonical' | grep -c 'gene_type=protein_coding;'
    chromosome = 'chr22'
    genome_annotation = chr22_example_files['genome_annotation']
    # not all genes have a canonical transcript if the genome annotation is not current (e.g. GRCh37 is not current)
    roi = get_tss_from_genome_annotation(genome_annotation, chromosome=chromosome, protein_coding_only=True,
                                         canonical_only=True)
    assert len(roi) == 441
    roi = get_tss_from_genome_annotation(genome_annotation, chromosome=chromosome, protein_coding_only=True,
                                         canonical_only=True, gene_ids=['ENSG00000172967'])
    assert roi['gene_id'][0] == 'ENSG00000172967.8_5'
    roi = get_tss_from_genome_annotation(genome_annotation, chromosome=chromosome, protein_coding_only=False,
                                         canonical_only=False)
    assert len(roi) == 5279
    roi = get_tss_from_genome_annotation(genome_annotation, chromosome=chromosome, protein_coding_only=True,
                                         canonical_only=False)
    assert len(roi) == 3623
    roi = get_tss_from_genome_annotation(genome_annotation, chromosome=chromosome, protein_coding_only=False,
                                         canonical_only=True)
    assert len(roi) == 1212



def test_read_gff3(tmp_path):
    gff3 = tmp_path / 'annotation.gff3'
    gff3.write_text('\n'.join([
        '##gff-version 3',
        'chr22\tHAVANA\ttranscript\t11\t20\t.\t-\t.\tID=ENST01.1;Parent=ENSG01.1;gene_id=ENSG01.1;'
        'transcript_id=ENST01.1;gene_type=protein_coding;tag=basic,Ensembl_canonical',
        # "%2541" is an escaped "%41", which must not become "A"
        'chrY\tHAVANA\ttranscript\t101\t200\t.\t+\t.\tID=ENST02.1_PAR_Y;Parent=ENSG02.1_PAR_Y;gene_id=ENSG02.1;'
        'transcript_id=ENST02.1;gene_type=lncRNA;tag=a%3Bb,c%2541',
        'chrY\tHAVANA\ttranscript\t301\t400\t.\t+\t.\tID=ENSTR03.1;Parent=ENSGR03.1;gene_id=ENSG03.1;'
        'transcript_id=ENST03.1;gene_type=lncRNA',
    ]) + '\n')
    annotation = read_gff3(gff3)
    assert annotation.select('Chromosome', 'Feature', 'Start', 'End', 'Strand').rows() == [
        ('chr22', 'transcript', 10, 20, '-'),
        ('chrY', 'transcript', 100, 200, '+'),
        ('chrY', 'transcript', 300, 400, '+'),
    ]
    # the chrY PAR copies get the suffix "_PAR_Y", as in the GENCODE GTF
    assert annotation['transcript_id'].to_list() == ['ENST01.1', 'ENST02.1_PAR_Y', 'ENST03.1_PAR_Y']
    assert annotation['gene_id'].to_list() == ['ENSG01.1', 'ENSG02.1_PAR_Y', 'ENSG03.1_PAR_Y']
    assert annotation['tag'].to_list() == ['basic,Ensembl_canonical', 'a;b,c%41', None]
    roi = get_tss_from_genome_annotation(gff3, protein_coding_only=True, canonical_only=True)
    assert roi['transcript_id'].to_list() == ['ENST01.1']
    with pytest.raises(ValueError, match='GTF'):
        read_gff3(tmp_path / 'annotation.gtf.gz')


def test_gtf_is_a_deprecated_alias(chr22_example_files):
    # the workflow scripts pass the genome annotation as gtf, also in a dict of dataloader arguments
    genome_annotation = pl.DataFrame({
        'Chromosome': ['chr22', 'chr22'],
        'Feature': ['transcript', 'transcript'],
        'Start': [20_000_000, 30_000_000],
        'End': [20_001_000, 30_001_000],
        'Strand': ['+', '-'],
        'gene_id': ['ENSG01.1', 'ENSG02.1'],
        'transcript_id': ['ENST01.1', 'ENST02.1'],
        'gene_type': ['protein_coding', 'protein_coding'],
        'tag': ['basic,Ensembl_canonical', 'basic'],
    })
    with pytest.warns(DeprecationWarning, match="'gtf'"):
        roi = get_tss_from_genome_annotation(gtf=genome_annotation)
    assert roi.equals(get_tss_from_genome_annotation(genome_annotation=genome_annotation))

    dl_args = {'fasta_file': chr22_example_files['fasta'], 'chromosome': 'chr22', 'gtf': genome_annotation}
    with pytest.warns(DeprecationWarning, match="'gtf'"):
        dl = TSSDataloader.from_allele_type(AlleleType.REF, **dl_args)
    assert len(dl) == 2

    with pytest.warns(DeprecationWarning, match="'gtf'"):
        veff = EnformerVeff(gtf=genome_annotation)
    assert veff.canonical_transcripts.to_list() == ['ENST01']

    with pytest.raises(TypeError, match="'gtf'"):
        get_tss_from_genome_annotation(genome_annotation=genome_annotation, gtf=genome_annotation)


def test_vcf_lazy_is_deprecated(chr22_example_files):
    # the workflow script passes vcf_lazy=True; kipoiseq2 0.1 has no eager VCF reading
    with pytest.warns(DeprecationWarning, match='vcf_lazy'):
        dl = VCFTSSDataloader(
            fasta_file=chr22_example_files['fasta'],
            genome_annotation=chr22_example_files['genome_annotation'],
            vcf_file=chr22_example_files['vcf'],
            vcf_lazy=True,
            seq_length=21,
            shifts=[0],
        )
    assert len(dl) > 0

def test_variant_regions_are_extended_strand_aware(chr22_example_files):
    genome_annotation = pl.DataFrame({
        'Chromosome': ['chr22', 'chr22'],
        'Feature': ['transcript', 'transcript'],
        'Start': [1000, 2000],
        'End': [1500, 2500],
        'Strand': ['+', '-'],
        'gene_id': ['g1', 'g2'],
        'transcript_id': ['t1', 't2'],
    })
    dl = VCFTSSDataloader(
        fasta_file=chr22_example_files['fasta'],
        genome_annotation=genome_annotation,
        vcf_file=chr22_example_files['vcf'],
        variant_upstream_tss=5,
        variant_downstream_tss=20,
        seq_length=21,
        shifts=[0],
    )
    regions = dl._get_single_variant_matcher().intervals
    # + strand: the TSS is 1000, upstream is towards lower positions
    # - strand: the TSS is 2499, upstream is towards higher positions
    assert regions.select('chrom', 'start', 'end', 'strand', 'tss').rows() == [
        ('chr22', 995, 1021, '+', 1000),
        ('chr22', 2479, 2505, '-', 2499),
    ]


def test_vcf_dataloader(chr22_example_files, variants):
    dl = VCFTSSDataloader(
        fasta_file=chr22_example_files['fasta'],
        genome_annotation=chr22_example_files['genome_annotation'],
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
        # seq_end is the 1-based stop, so the window has 21 bases
        assert metadata['seq_end'] - metadata['seq_start'] == 21
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
        genome_annotation=chr22_example_files['genome_annotation'],
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
            assert (metadata['seq_start'], metadata['seq_end']) == (ref['start'], ref['end'])
            checked_refs[ref_id] = 2

    # check that all variants in my list were found and checked
    assert set(checked_refs.keys()) == set(references.keys())
    print(total)
