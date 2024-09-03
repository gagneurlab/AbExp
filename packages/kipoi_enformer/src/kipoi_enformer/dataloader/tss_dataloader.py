from abc import ABC, abstractmethod

import pandas as pd
from kipoiseq.extractors import VariantSeqExtractor, SingleVariantMatcher, FastaStringExtractor
import pyranges as pr
from kipoiseq.extractors import MultiSampleVCF
import pyarrow as pa
import numpy as np
from .dataloader import Dataloader, get_tss_from_genome_annotation, extract_sequences_around_anchor
from kipoi_enformer.constants import AlleleType
from kipoi_enformer.logger import logger

__all__ = ['TSSDataloader', 'RefTSSDataloader', 'VCFTSSDataloader']

# length of sequence which enformer gets as input
# ═════┆═════┆════════════════════════┆═════┆═════
ENFORMER_SEQUENCE_LENGTH = 393_216


class TSSDataloader(Dataloader):
    def __init__(self, allele_type: AlleleType, fasta_file, gtf: pd.DataFrame | str, chromosome: str | None = None,
                 seq_length: int = ENFORMER_SEQUENCE_LENGTH, shifts: list[int] = (-43, 0, 43), size: int = None,
                 canonical_only: bool = False,
                 protein_coding_only: bool = False, gene_ids: list | None = None,
                 *args, **kwargs):
        """

        :param fasta_file: Fasta file with the reference genome
        :param gtf: GTF file with genome annotation or DataFrame with genome annotation
        :param chromosome: The chromosome to filter for. If None, all chromosomes are used.
        :param seq_length: The length of the sequence to return.
        :param shifts: The shifts in relation to the TSS.
        :param size: The number of samples to return. If None, all samples are returned.
        :param canonical_only: If True, only Ensembl canonical transcripts are extracted from the genome annotation
        :param protein_coding_only: If True, only protein coding transcripts are extracted from the genome annotation
        :param gene_ids: If provided, only the gene with this ID is extracted from the genome annotation
        """

        super().__init__(fasta_file=fasta_file, size=size, *args, **kwargs)
        for shift in shifts:
            assert abs(shift) < seq_length, f"shift must be smaller than seq_length but got {shift} >= {seq_length}"

        self._canonical_only = canonical_only
        self._protein_coding_only = protein_coding_only
        self._seq_length = seq_length
        self.chromosome = chromosome
        logger.debug(f"Loading genome annotation")
        self._genome_annotation = get_tss_from_genome_annotation(gtf, chromosome=self.chromosome,
                                                                 canonical_only=canonical_only,
                                                                 protein_coding_only=protein_coding_only,
                                                                 gene_ids=gene_ids)
        self._shifts = shifts
        self.metadata = {'shifts': ';'.join([str(x) for x in self._shifts]), 'allele_type': allele_type.value,
                         'seq_length': str(self._seq_length)}

    @classmethod
    def from_allele_type(cls, allele_type: AlleleType, *args, **kwargs):
        if allele_type == AlleleType.REF:
            return RefTSSDataloader(*args, **kwargs)
        elif allele_type == AlleleType.ALT:
            return VCFTSSDataloader(*args, **kwargs)
        else:
            raise ValueError(f"Unknown allele type: {allele_type}")


class RefTSSDataloader(TSSDataloader):
    def __init__(self, fasta_file, gtf: pd.DataFrame | str, chromosome: str,
                 seq_length: int = ENFORMER_SEQUENCE_LENGTH, shifts: list[int] = (-43, 0, 43), size: int = None,
                 canonical_only: bool = False,
                 protein_coding_only: bool = False, gene_ids: list | None = None, *args, **kwargs):
        """
        :param fasta_file: Fasta file with the reference genome
        :param gtf: GTF file with genome annotation or DataFrame with genome annotation
        :param chromosome: The chromosome to filter for
        :param seq_length: The length of the sequence to return.
        :param shifts: The shifts in relation to the TSS.
        :param size: The number of samples to return. If None, all samples are returned.
        :param canonical_only: If True, only Ensembl canonical transcripts are extracted from the genome annotation
        :param protein_coding_only: If True, only protein coding transcripts are extracted from the genome annotation
        :param gene_id: If provided, only the gene with this ID is extracted from the genome annotation
        """
        assert chromosome is not None, 'A chromosome should be provided'
        super().__init__(AlleleType.REF, chromosome=chromosome, fasta_file=fasta_file, gtf=gtf,
                         seq_length=seq_length, shifts=shifts, size=size, canonical_only=canonical_only,
                         protein_coding_only=protein_coding_only, gene_ids=gene_ids, *args, **kwargs)
        logger.debug(f"Dataloader is ready for chromosome {chromosome}")

    def _sample_gen(self):
        for _, row in self._genome_annotation.iterrows():
            try:
                chromosome = row['Chromosome']
                strand = row.get('Strand', '.')
                tss = row['tss']

                sequences, interval = extract_sequences_around_anchor(self._shifts, chromosome, strand, tss,
                                                                      self._seq_length,
                                                                      ref_seq_extractor=self._reference_sequence)

                metadata = {
                    "seq_start": interval.start,  # 0-based start of the input sequence
                    "seq_end": interval.end + 1,  # 1-based stop of the input sequence
                    "tss": tss,  # 0-based position of the TSS
                    "strand": strand,
                    "gene_id": row['gene_id'],
                    "transcript_id": row['transcript_id'],
                    "transcript_start": row['transcript_start'],  # 0-based
                    "transcript_end": row['transcript_end'],  # 1-based
                }

                yield metadata, np.stack(sequences)
            except Exception as e:
                logger.error(f"Error processing row: {row}")
                raise e

    def __len__(self):
        if self._genome_annotation is None:
            return 0
        return len(self._genome_annotation) if self._size is None else min(self._size, len(self._genome_annotation))

    @property
    def pyarrow_metadata_schema(self):
        """
        Get the pyarrow schema for the metadata.
        :return: PyArrow schema for the metadata
        """
        columns = [
            ('seq_start', pa.int64()),
            ('seq_end', pa.int64()),
            ('tss', pa.int64()),
            ('strand', pa.string()),
            ('gene_id', pa.string()),
            ('transcript_id', pa.string()),
            ('transcript_start', pa.int64()),
            ('transcript_end', pa.int64()), ]

        return pa.schema(columns, metadata=self.metadata)


class VCFTSSDataloader(TSSDataloader):
    def __init__(self, fasta_file, gtf: pd.DataFrame | str, vcf_file, vcf_lazy=True,
                 variant_upstream_tss: int = 10, variant_downstream_tss: int = 10,
                 seq_length: int = ENFORMER_SEQUENCE_LENGTH, shifts: list[int] = (-43, 0, 43),
                 size: int = None, canonical_only: bool = False, protein_coding_only: bool = False,
                 gene_ids: list | None = None, *args, **kwargs):
        """

        :param fasta_file: Fasta file with the reference genome
        :param gtf: GTF file with genome annotation or DataFrame with genome annotation
        :param vcf_file: VCF file with variants
        :param vcf_lazy: If True, the VCF file is read lazily
        :param variant_upstream_tss: The number of bases upstream the TSS to look for variants
        :param variant_downstream_tss: The number of bases downstream the TSS to look for variants
        :param seq_length: The length of the sequence to return.
        :param shifts: The shifts in relation to the TSS.
        :param size: The number of samples to return. If None, all samples are returned.
        :param canonical_only: If True, only Ensembl canonical transcripts are extracted from the genome annotation
        :param protein_coding_only: If True, only protein coding transcripts are extracted from the genome annotation
        :param gene_id: If provided, only the gene with this ID is extracted from the genome annotation
        """

        super().__init__(AlleleType.ALT, fasta_file=fasta_file, gtf=gtf, chromosome=None,
                         seq_length=seq_length, shifts=shifts, size=size, canonical_only=canonical_only,
                         protein_coding_only=protein_coding_only, gene_ids=gene_ids, *args, **kwargs)
        for shift in shifts:
            assert abs(shift) < variant_downstream_tss + variant_upstream_tss + 1, \
                (f"shift must be smaller than downstream_tss + upstream_tss + 1 but got "
                 f"{shift} >= {variant_downstream_tss + variant_upstream_tss + 1}")

        self._variant_seq_extractor = VariantSeqExtractor(reference_sequence=self._reference_sequence)
        self.vcf_file = vcf_file
        self.vcf_lazy = vcf_lazy
        self.variant_upstream_tss = variant_upstream_tss
        self.variant_downstream_tss = variant_downstream_tss
        logger.debug(f"Dataloader is ready")

    def _sample_gen(self):
        for interval, variant in self._get_single_variant_matcher(self.vcf_lazy):
            try:
                attrs = interval.attrs
                tss = attrs['tss']
                chromosome = interval.chrom
                strand = interval.strand

                sequences, interval = extract_sequences_around_anchor(self._shifts, chromosome, strand, tss,
                                                                      self._seq_length,
                                                                      ref_seq_extractor=self._reference_sequence,
                                                                      variant_extractor=self._variant_seq_extractor,
                                                                      variant=variant)
                metadata = {
                    "seq_start": interval.start,  # 0-based start of the input sequence
                    "seq_end": interval.end + 1,  # 1-based stop of the input sequence
                    "tss": tss,  # 0-based position of the TSS
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
                yield metadata, np.stack(sequences)
            except Exception as e:
                logger.error(f"Error processing variant-interval")
                logger.error(f"Interval: {interval}")
                logger.error(f"Variant: {variant}")
                raise e

    def __len__(self):
        if self._genome_annotation is None or len(self._genome_annotation) == 0:
            return 0
        tmp_matcher = self._get_single_variant_matcher(vcf_lazy=False)
        total = sum(1 for _, _ in tmp_matcher)
        if self._size:
            return min(self._size, total)
        return total

    def _get_single_variant_matcher(self, vcf_lazy=True):
        if self._genome_annotation is None or len(self._genome_annotation) == 0:
            return iter([])
        # reads the genome annotation
        # start and end are transformed to 0-based and 1-based respectively
        roi = pr.PyRanges(self._genome_annotation)
        roi = roi.extend(ext={"5": self.variant_upstream_tss, "3": self.variant_downstream_tss})
        # todo do assert length of roi

        interval_attrs = ['gene_id', 'transcript_id', 'tss', 'transcript_start', 'transcript_end']
        for attr in interval_attrs:
            assert attr in roi.columns, f"attr must be in {roi.columns}"
        variants = MultiSampleVCF(self.vcf_file, lazy=vcf_lazy)

        return SingleVariantMatcher(
            variant_fetcher=variants,
            pranges=roi,
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
