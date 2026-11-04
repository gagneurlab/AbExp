import polars as pl
import pyarrow as pa
from kipoiseq2.extractors import SingleVariantMatcher

from abexp.promoterai.constants import DEFAULT_INPUT_LENGTH
from abexp.promoterai.logger import logger
from .dataloader import Dataloader, extract_variant_sequences, get_tss_from_genome_annotation

__all__ = ['VCFTSSDataloader']


class VCFTSSDataloader(Dataloader):
    """
    The variants of a VCF file that lie in the promoter of a transcript, with the sequences that PromoterAI scores.

    One record per variant and transcript: a variant in the promoters of several transcripts is scored once per
    transcript. The sequences depend only on the variant and on the strand of the transcript, see
    `abexp.promoterai.dataloader.dataloader.extract_variant_sequences`, so the records of two transcripts on the
    same strand hold the same sequences.
    """

    def __init__(self, fasta_file, genome_annotation, vcf_file,
                 variant_upstream_tss: int = 500, variant_downstream_tss: int = 500,
                 seq_length: int = DEFAULT_INPUT_LENGTH, size: int = None,
                 canonical_only: bool = False, protein_coding_only: bool = False,
                 gene_ids: list | None = None, *args, **kwargs):
        """
        :param fasta_file: Fasta file with the reference genome
        :param genome_annotation: GFF3 file or DataFrame with the genome annotation,
            see `abexp.promoterai.utils.genome_annotation_to_polars`.
        :param vcf_file: VCF file with variants
        :param variant_upstream_tss: The number of bases upstream the TSS to look for variants
        :param variant_downstream_tss: The number of bases downstream the TSS to look for variants
        :param seq_length: The length of the sequences to return, the input length of the model.
        :param size: The number of samples to return. If None, all samples are returned.
        :param canonical_only: If True, only Ensembl canonical transcripts are extracted from the genome annotation
        :param protein_coding_only: If True, only protein coding transcripts are extracted from the genome annotation
        :param gene_ids: If provided, only the genes with these IDs are extracted from the genome annotation
        """
        super().__init__(fasta_file=fasta_file, size=size, *args, **kwargs)
        assert seq_length > 0, f'seq_length must be positive but got {seq_length}'

        self._seq_length = seq_length
        logger.debug(f"Loading genome annotation")
        self._genome_annotation = get_tss_from_genome_annotation(genome_annotation, chromosome=None,
                                                                 canonical_only=canonical_only,
                                                                 protein_coding_only=protein_coding_only,
                                                                 gene_ids=gene_ids)
        self.vcf_file = vcf_file
        self.variant_upstream_tss = variant_upstream_tss
        self.variant_downstream_tss = variant_downstream_tss
        self.metadata = {'seq_length': str(self._seq_length),
                         'variant_upstream_tss': str(self.variant_upstream_tss),
                         'variant_downstream_tss': str(self.variant_downstream_tss)}
        logger.debug(f"Dataloader is ready")

    def _sample_gen(self):
        for interval, variant in self._get_single_variant_matcher():
            try:
                attrs = interval.attrs
                sequences, seq_interval, scorable = extract_variant_sequences(
                    interval.chrom, interval.strand, variant, self._seq_length,
                    ref_seq_extractor=self._reference_sequence)
                metadata = {
                    "seq_start": seq_interval.start,  # 0-based start of the input sequence
                    "seq_end": seq_interval.end,  # 1-based stop of the input sequence
                    "tss": attrs['tss'],  # 0-based position of the TSS
                    "chrom": interval.chrom,
                    "strand": interval.strand,
                    "gene_id": attrs['gene_id'],
                    "transcript_id": attrs['transcript_id'],
                    "transcript_start": attrs['transcript_start'],  # 0-based
                    "transcript_end": attrs['transcript_end'],  # 1-based
                    "variant_start": variant.start,  # 0-based
                    "variant_end": variant.end,  # 1-based
                    "ref": variant.ref,
                    "alt": variant.alt,
                }
                yield metadata, sequences, scorable
            except Exception as e:
                logger.error(f"Error processing variant-interval")
                logger.error(f"Interval: {interval}")
                logger.error(f"Variant: {variant}")
                raise e

    def __len__(self):
        if self._genome_annotation is None or len(self._genome_annotation) == 0:
            return 0
        total = self._get_single_variant_matcher().pairs().height
        if self._size:
            return min(self._size, total)
        return total

    def _get_single_variant_matcher(self):
        if self._genome_annotation is None or len(self._genome_annotation) == 0:
            return iter([])
        interval_attrs = ['gene_id', 'transcript_id', 'tss', 'transcript_start', 'transcript_end']
        for attr in interval_attrs:
            assert attr in self._genome_annotation.columns, f"attr must be in {self._genome_annotation.columns}"

        # extend the TSS by variant_upstream_tss bases upstream and variant_downstream_tss bases downstream
        # (0-based start, 1-based end)
        minus = pl.col('Strand') == '-'
        upstream = self.variant_upstream_tss
        downstream = self.variant_downstream_tss
        regions = self._genome_annotation.select(
            pl.col('Chromosome').alias('chrom'),
            (pl.col('Start') - pl.when(minus).then(downstream).otherwise(upstream)).clip(lower_bound=0).alias('start'),
            (pl.col('End') + pl.when(minus).then(upstream).otherwise(downstream)).alias('end'),
            pl.col('Strand').alias('strand'),
            *interval_attrs,
        )
        return SingleVariantMatcher(
            vcf_file=self.vcf_file,
            intervals=regions,
            interval_attrs=interval_attrs
        )

    @property
    def pyarrow_metadata_schema(self):
        """
        Get the pyarrow schema for the metadata.
        :return: PyArrow schema for the metadata
        """
        return pa.schema(
            [
                ('seq_start', pa.int64()),
                ('seq_end', pa.int64()),
                ('tss', pa.int64()),
                ('chrom', pa.string()),
                ('strand', pa.string()),
                ('gene_id', pa.string()),
                ('transcript_id', pa.string()),
                ('transcript_start', pa.int64()),
                ('transcript_end', pa.int64()),
                ('variant_start', pa.int64()),
                ('variant_end', pa.int64()),
                ('ref', pa.string()),
                ('alt', pa.string()),
            ], metadata=self.metadata)
