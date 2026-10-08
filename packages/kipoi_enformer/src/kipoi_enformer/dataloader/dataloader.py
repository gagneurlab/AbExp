from abc import ABC, abstractmethod
from collections.abc import Mapping

import numpy as np
import pyarrow as pa
import math
import polars as pl
from kipoiseq2.extractors import VariantSeqExtractor, FastaStringExtractor
from kipoiseq2 import Interval, Variant
from kipoiseq2.transforms.functional import one_hot_dna
from kipoi_enformer.utils import genome_annotation_to_polars


def numpy_collate(samples: list):
    """
    Collate a list of samples into one batch.
    Dictionaries are collated per key, numpy arrays are stacked along a new first axis,
    and all other values (numbers, strings) are converted to a numpy array.
    :param samples: List of samples with the same structure
    :return: The batch, with the same structure as a single sample
    """
    first = samples[0]
    if isinstance(first, Mapping):
        return {key: numpy_collate([sample[key] for sample in samples]) for key in first}
    if isinstance(first, np.ndarray):
        return np.stack(samples)
    return np.asarray(samples)


class Dataloader(ABC):
    def __init__(self, fasta_file, size: int = None, *args, **kwargs):
        """

        :param fasta_file: Fasta file with the reference genome
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

    def batch_iter(self, batch_size: int):
        """
        Iterate over the dataset in batches.

        :param batch_size: The number of samples per batch. The last batch can be smaller.
        :return: Iterator over the batches. Each batch has the same keys as a sample,
            with the values of all samples in the batch collated by `numpy_collate`.
        """
        batch = []
        for sample in self:
            batch.append(sample)
            if len(batch) == batch_size:
                yield numpy_collate(batch)
                batch = []
        if len(batch) > 0:
            yield numpy_collate(batch)


def get_tss_from_genome_annotation(genome_annotation, chromosome: str | None = None,
                                   protein_coding_only: bool = False, canonical_only: bool = False,
                                   gene_ids: list | None = None) -> pl.DataFrame:
    """
    Get TSS from genome annotation
    :param genome_annotation: GFF3 file or DataFrame with the genome annotation, see `genome_annotation_to_polars`
    :return: genome_annotation with Start and End set to the TSS
        and the additional columns tss (0-based), transcript_start (0-based), transcript_end (1-based)
    """
    roi = get_roi_from_genome_annotation(genome_annotation, chromosome, protein_coding_only, canonical_only, gene_ids)
    # the TSS of a transcript on the minus strand is its last base
    tss = pl.when(pl.col('Strand') == '-').then(pl.col('End') - 1).otherwise(pl.col('Start'))
    return roi.with_columns(Start=tss, End=tss + 1, tss=tss)


def get_roi_from_genome_annotation(genome_annotation, chromosome: str | None = None,
                                   protein_coding_only: bool = False, canonical_only: bool = False,
                                   gene_ids: list | None = None) -> pl.DataFrame:
    """
    Get ROI from genome annotation
    :param genome_annotation: GFF3 file or DataFrame with the genome annotation, see `genome_annotation_to_polars`
    :return: the transcripts of the genome annotation, filtered,
        with the additional columns transcript_start (0-based), transcript_end (1-based)
    """
    roi = genome_annotation_to_polars(genome_annotation)
    if gene_ids is not None:
        roi = roi.filter(pl.col('gene_id').str.contains('|'.join(gene_ids)))
    if chromosome is not None:
        roi = roi.filter(pl.col('Chromosome') == chromosome)
    roi = roi.filter(pl.col('Feature') == 'transcript')
    if protein_coding_only:
        roi = roi.filter(pl.col('gene_type') == 'protein_coding')
    if canonical_only:
        # check if Ensembl_canonical is in the set of tags
        roi = roi.filter(pl.col('tag').str.split(',').list.contains('Ensembl_canonical').fill_null(False))
    return roi.with_columns(transcript_start=pl.col('Start'), transcript_end=pl.col('End'))


def construct_interval(chrom, strand, anchor, seq_length):
    # input interval without shift
    # if the sequence length is even, the tss is closer to the end of the sequence
    # if the sequence length is odd, the tss is in the middle of the sequence
    five_end_len = math.floor(seq_length / 2)
    three_end_len = math.ceil(seq_length / 2)

    # kipoiseq2.Interval is 0-based and half-open
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

    interval = construct_interval(chromosome, strand, anchor, seq_length)
    sequences = []
    # shift intervals and extract sequences
    for shift in shifts:
        shifted_interval = interval.shift(shift, use_strand=True)
        five_end_pad = 0
        three_end_pad = 0
        if shifted_interval.start < 0:
            five_end_pad = abs(shifted_interval.start)
        # the interval is half-open, so it may end at chrom_len
        if shifted_interval.end > chrom_len:
            three_end_pad = shifted_interval.end - chrom_len
        if five_end_pad > 0 or three_end_pad > 0:
            shifted_interval = shifted_interval.truncate(chrom_len)

        if variant is not None:
            seq = variant_extractor.extract(shifted_interval,
                                            [variant],
                                            anchor=anchor,
                                            fixed_len=True,
                                            is_padding=True,
                                            chrom_len=chrom_len,
                                            )
        else:
            seq = ref_seq_extractor.extract(shifted_interval)

        # the extractors reverse complement minus-strand sequences,
        # so the padding of the interval start belongs to the end of the sequence
        if strand == '-':
            seq = 'N' * three_end_pad + seq + 'N' * five_end_pad
        else:
            seq = 'N' * five_end_pad + seq + 'N' * three_end_pad

        assert len(seq) == seq_length, \
            f"interval width must be {seq_length} but got {len(seq)}"

        sequences.append(one_hot_dna(seq))
    return sequences, interval
