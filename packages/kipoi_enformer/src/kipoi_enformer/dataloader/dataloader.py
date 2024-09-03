from abc import ABC, abstractmethod

from kipoi.data import SampleGenerator
import pyarrow as pa
import math
import pandas as pd
from kipoiseq.extractors import VariantSeqExtractor, FastaStringExtractor
from kipoiseq import Interval, Variant
from kipoiseq.transforms.functional import one_hot_dna
from kipoi_enformer.utils import gtf_to_pandas


class Dataloader(SampleGenerator, ABC):
    def __init__(self, fasta_file, size: int = None, *args, **kwargs):
        """

        :param fasta_file: Fasta file with the reference genome
        :param gtf: GTF file with genome annotation or DataFrame with genome annotation
        :param chromosome: The chromosome to filter for. If None, all chromosomes are used.
        :param seq_length: The length of the sequence to return. This should be the length of the Enformer input sequence.
        :param shift: For each sequence, we have 3 shifts, -shift, 0, shift, in relation to a reference point.
        :param size: The number of samples to return. If None, all samples are returned.
        :param canonical_only: If True, only Ensembl canonical transcripts are extracted from the genome annotation
        :param protein_coding_only: If True, only protein coding transcripts are extracted from the genome annotation
        :param gene_ids: If provided, only the gene with this ID is extracted from the genome annotation
        """

        super().__init__(*args, **kwargs)
        self._reference_sequence = FastaStringExtractor(fasta_file, use_strand=True)
        if not self._reference_sequence.use_strand:
            raise ValueError(
                "Reference sequence fetcher does not use strand but this is needed to obtain correct sequences!")
        self._size = size

    @abstractmethod
    def __len__(self):
        raise NotImplementedError("The length of the dataset is not known.")

    @abstractmethod
    def _sample_gen(self):
        """
        Generate samples for the dataset. The generator should return a tuple of metadata and sequences.
        :return:
        """
        raise NotImplementedError("The sample generator is not implemented.")

    @property
    @abstractmethod
    def pyarrow_metadata_schema(self) -> pa.schema:
        """
        Get the pyarrow schema for the metadata.
        :return: PyArrow schema for the metadata
        """
        raise NotImplementedError("The metadata schema is not implemented.")

    def __iter__(self):
        """
        Iterate over the dataset.

        :return: Iterator over the dataset. Each item is a dictionary with the following
            keys:
            - sequences:
            - metadata:
        """
        counter = 0
        for metadata, sequences in self._sample_gen():
            # check if we reached the end of the dataset
            if self._size is not None and counter == self._size:
                break
            counter += 1

            yield {
                "sequences": sequences,
                "metadata": metadata
            }


def get_tss_from_genome_annotation(gtf: pd.DataFrame | str, chromosome: str | None = None,
                                   protein_coding_only: bool = False, canonical_only: bool = False,
                                   gene_ids: list | None = None):
    """
    Get TSS from genome annotation
    :return: genome_annotation with additional columns tss (0-based), transcript_start (0-based), transcript_end (1-based)
    """
    roi = get_roi_from_genome_annotation(gtf, chromosome, protein_coding_only, canonical_only, gene_ids)

    def adjust_row(row):
        if row.Strand == '-':
            # convert 1-based to 0-based
            tss = row.End - 1
        else:
            tss = row.Start

        row.Start = tss
        row.End = tss + 1
        return row
    if len(roi) > 0:
        roi = roi.apply(adjust_row, axis=1)
        roi['tss'] = roi["Start"]
    return roi


def get_roi_from_genome_annotation(gtf: pd.DataFrame | str, chromosome: str | None = None,
                                   protein_coding_only: bool = False, canonical_only: bool = False,
                                   gene_ids: list | None = None):
    """
    Get ROI from genome annotation
    :return: filtered genome_annotation
    """
    if not isinstance(gtf, pd.DataFrame):
        genome_annotation = gtf_to_pandas(gtf)
    else:
        genome_annotation = gtf.copy()
    if gene_ids is not None:
        genome_annotation = genome_annotation[genome_annotation['gene_id'].str.contains('|'.join(gene_ids))]
    if chromosome is not None:
        genome_annotation = genome_annotation.query("`Chromosome` == @chromosome")
    roi = genome_annotation.query("`Feature` == 'transcript'")
    if protein_coding_only:
        roi = roi.query("`gene_type` == 'protein_coding'")
    if canonical_only:
        # check if Ensembl_canonical is in the set of tags
        roi = roi[roi['tag'].apply(lambda x: False if pd.isna(x) else ('Ensembl_canonical' in x.split(',')))]
    if len(roi) > 0:
        roi = roi.assign(
            transcript_start=roi["Start"],
            transcript_end=roi["End"],
        )

    return roi


def construct_interval(chrom, strand, anchor, seq_length):
    # input interval without shift
    # if the sequence length is even, the tss is closer to the end of the sequence
    # if the sequence length is odd, the tss is in the middle of the sequence
    five_end_len = math.floor(seq_length / 2)
    three_end_len = math.ceil(seq_length / 2)

    # WARNING: kipoiseq.Interval has a 0-based end!
    interval = Interval(chrom=chrom,
                        start=anchor - five_end_len,
                        end=anchor + three_end_len,
                        strand=strand)

    assert (interval.width()) == seq_length, \
        f"interval width must be {seq_length} but got {interval.width()}"
    assert (anchor - interval.start) == seq_length // 2, \
        f"tss must be in the middle of the interval but got {anchor - interval.start}"
    return interval


def extract_sequences_around_anchor(shifts, chromosome, strand, anchor, seq_length,
                                    ref_seq_extractor: FastaStringExtractor,
                                    variant_extractor: VariantSeqExtractor | None = None,
                                    variant: Variant | None = None):
    assert variant_extractor is None or (variant is not None and variant_extractor is not None), \
        "variant_extractor must be provided if variant is not None"
    chrom_len = len(ref_seq_extractor.fasta.records[chromosome])

    # WARNING: kipoiseq.Interval has a 0-based end!
    interval = construct_interval(chromosome, strand, anchor, seq_length)
    sequences = []
    # shift intervals and extract sequences
    for shift in shifts:
        shifted_interval = interval.shift(shift, use_strand=True)
        five_end_pad = 0
        three_end_pad = 0
        if shifted_interval.start < 0:
            five_end_pad = abs(shifted_interval.start)
        if shifted_interval.end >= chrom_len:
            three_end_pad = shifted_interval.end - chrom_len + 1
        if five_end_pad > 0 or three_end_pad > 0:
            shifted_interval = shifted_interval.truncate(chrom_len)

        if variant is not None:
            seq = variant_extractor.extract(shifted_interval,
                                            [variant],
                                            anchor=anchor,
                                            fixed_length=True,
                                            is_padding=True,
                                            chrom_len=chrom_len,
                                            )
        else:
            seq = ref_seq_extractor.extract(shifted_interval)

        if five_end_pad > 0:
            seq = 'N' * five_end_pad + seq

        if three_end_pad > 0:
            seq = seq + 'N' * three_end_pad

        assert len(seq) == seq_length, \
            f"interval width must be {seq_length} but got {len(seq)}"

        sequences.append(one_hot_dna(seq))
    return sequences, interval
