from abc import ABC, abstractmethod
from collections.abc import Mapping

import numpy as np
import polars as pl
import pyarrow as pa
from kipoiseq2 import Interval, Variant
from kipoiseq2.extractors import FastaStringExtractor
from kipoiseq2.transforms.functional import one_hot_dna

from abexp.promoterai.logger import logger
from abexp.promoterai.utils import genome_annotation_to_polars

# the bases that PromoterAI one-hot encodes; it encodes every other base, e.g. N, as zeros
BASES = 'ACGT'
# translation tables of `reverse_complement` and `one_hot`: a base other than ACGT, e.g. N or an IUPAC code, becomes N
_COMPLEMENT = {i: dict(zip(BASES, 'TGCA')).get(chr(i), 'N') for i in range(256)}
_SANITIZE = {i: chr(i) if chr(i) in BASES else 'N' for i in range(256)}


def numpy_collate(samples: list):
    """
    Collate a list of samples into one batch.
    Dictionaries are collated per key, numpy arrays are stacked along a new first axis,
    and all other values (numbers, strings, booleans) are converted to a numpy array.
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
        :param size: The number of samples to return. If None, all samples are returned.
        """
        super().__init__(*args, **kwargs)
        # PromoterAI builds the sequence of the plus strand and reverse complements it itself, see
        # extract_variant_sequences, so the extractor must not use the strand of the interval.
        self._reference_sequence = FastaStringExtractor(fasta_file, use_strand=False)
        self._size = size

    @abstractmethod
    def __len__(self):
        raise NotImplementedError("The length of the dataset is not known.")

    @abstractmethod
    def _sample_gen(self):
        """
        Generate samples for the dataset. The generator should return a tuple of metadata, sequences and a flag
        whether the variant could be scored.
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
            - sequences: dict of the one-hot encoded "ref" and "alt" sequences
            - scorable: whether the sequences hold the variant, see `extract_variant_sequences`
            - metadata: dict of the metadata of the record
        """
        counter = 0
        for metadata, sequences, scorable in self._sample_gen():
            # check if we reached the end of the dataset
            if self._size is not None and counter == self._size:
                break
            counter += 1

            yield {
                "sequences": sequences,
                "scorable": scorable,
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


def reverse_complement(sequence: str) -> str:
    """The reverse complement of a DNA sequence; a base other than ACGT becomes N."""
    return sequence.translate(_COMPLEMENT)[::-1]


def one_hot(sequence: str) -> np.ndarray:
    """
    One-hot encode a DNA sequence in the order ACGT; a base other than ACGT, e.g. N, becomes zeros, as in PromoterAI.

    :param sequence: the sequence
    :return: float32 array of shape (len(sequence), 4)
    """
    return np.asarray(one_hot_dna(sequence.translate(_SANITIZE)), dtype=np.float32)


def construct_interval(chrom, strand, variant_start, seq_length):
    """
    The input interval of a variant: `seq_length` bases with the first base of the variant in the middle.

    PromoterAI centers the sequence on the first base of the variant, also for an indel: with `pos` the 1-based
    position of the variant, it reads the bases `pos - seq_length // 2` to `pos + seq_length // 2 - 1`. So the
    variant starts at the index `seq_length // 2` of the sequence.

    :param chrom: chromosome of the variant
    :param strand: strand of the TSS; the interval covers the same bases on both strands
    :param variant_start: 0-based start of the variant
    :param seq_length: the input length of the model
    :return: kipoiseq2.Interval, 0-based and half-open
    """
    five_end_len = seq_length // 2
    interval = Interval(chrom=chrom,
                        start=variant_start - five_end_len,
                        end=variant_start - five_end_len + seq_length,
                        strand=strand)
    assert interval.width() == seq_length, f"interval width must be {seq_length} but got {interval.width()}"
    return interval


def extract_variant_sequences(chromosome: str, strand: str, variant: Variant, seq_length: int,
                              ref_seq_extractor: FastaStringExtractor):
    """
    Extract the reference and the alternative sequence around a variant, as PromoterAI does.

    PromoterAI reads `seq_length` bases of the reference around the first base of the variant, see
    `construct_interval`, and replaces the reference allele in the middle by the alternative allele. An insertion
    pushes the bases at the 3' end out of the sequence, and a deletion pads the 3' end with N, so that both
    sequences have `seq_length` bases. For a TSS on the minus strand, both sequences are reverse complemented,
    which is what the models of PromoterAI see.

    A variant is not scorable, and both sequences are all zeros, if
    - its input interval reaches beyond the chromosome: PromoterAI skips such a variant,
    - its reference allele does not match the reference genome, e.g. in a VCF of another assembly,
    - or its alternative allele holds a base other than ACGT, e.g. a symbolic allele such as <DEL>.

    :param chromosome: chromosome of the variant
    :param strand: strand of the TSS, "+" or "-"
    :param variant: the variant, with a 0-based start
    :param seq_length: the input length of the model
    :param ref_seq_extractor: extractor of the reference genome, with `use_strand=False`
    :return: tuple of the sequences, a dict of the one-hot encoded "ref" and "alt" sequence, the input interval,
        and whether the variant is scorable
    """
    interval = construct_interval(chromosome, strand, variant.start, seq_length)
    chrom_len = len(ref_seq_extractor.fasta.records[chromosome])
    ref_allele = variant.ref.upper()
    alt_allele = variant.alt.upper()

    def zeros():
        return {name: np.zeros((seq_length, len(BASES)), dtype=np.float32) for name in ('ref', 'alt')}

    # the interval is half-open, so it may end at chrom_len
    if interval.start < 0 or interval.end > chrom_len:
        logger.warning(f'Skipping {variant}: its input interval {interval} reaches beyond the chromosome')
        return zeros(), interval, False
    if not set(alt_allele).issubset(set(BASES)):
        logger.warning(f'Skipping {variant}: its alternative allele is not a sequence of {BASES}')
        return zeros(), interval, False

    # the extractor does not use the strand, so the sequence is that of the plus strand
    ref_seq = ref_seq_extractor.extract(interval).upper()
    assert len(ref_seq) == seq_length, f"sequence length must be {seq_length} but got {len(ref_seq)}"

    variant_idx = seq_length // 2
    if ref_seq[variant_idx:variant_idx + len(ref_allele)] != ref_allele:
        logger.warning(f'Skipping {variant}: its reference allele does not match the reference genome')
        return zeros(), interval, False

    # an insertion drops the bases beyond the 3' end of the interval, a deletion pads it with N
    alt_seq = ref_seq[:variant_idx] + alt_allele + ref_seq[variant_idx + len(ref_allele):]
    alt_seq = alt_seq.ljust(seq_length, 'N')[:seq_length]

    if strand == '-':
        ref_seq = reverse_complement(ref_seq)
        alt_seq = reverse_complement(alt_seq)

    sequences = {'ref': one_hot(ref_seq), 'alt': one_hot(alt_seq)}
    return sequences, interval, True
