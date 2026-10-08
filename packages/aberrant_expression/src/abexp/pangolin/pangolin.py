import logging
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import polars as pl
from kipoiseq2 import Interval, Variant
from kipoiseq2.extractors import FastaStringExtractor, SingleVariantMatcher, VariantSeqExtractor, scan_vcf_variants
from kipoiseq2.transforms.functional import one_hot_dna

logger = logging.getLogger(__name__)

NO_SITES_WARNING = 'NoAnnotatedSitesToMaskForThisGene'


@dataclass(frozen=True)
class _Window:
    """The sequences around a variant, as the model of one strand sees them."""
    variant_idx: int
    strand: str
    ref: str
    alt: str


def align_alt(alt, ref_len, distance):
    """Align the outputs of the alt sequence of an indel to the positions of the ref sequence, as Pangolin does.

    The variant starts at index `distance`. For a deletion, zeros fill the deleted positions after the first base.
    For an insertion, the first base and the inserted bases collapse into their maximum.

    Args:
      alt: outputs of the alt sequence in genomic order, array with the positions in the last axis
      ref_len: number of outputs of the ref sequence
      distance: number of positions on each side of the variant

    Returns:
      array with `ref_len` positions in the last axis
    """
    alt_len = alt.shape[-1]
    if alt_len < ref_len:
        zeros = np.zeros((*alt.shape[:-1], ref_len - alt_len), dtype=alt.dtype)
        return np.concatenate([alt[..., :distance + 1], zeros, alt[..., distance + 1:]], axis=-1)
    if alt_len > ref_len:
        inserted = alt[..., distance:distance + alt_len - ref_len + 1].max(axis=-1, keepdims=True)
        return np.concatenate([alt[..., :distance], inserted, alt[..., distance + alt_len - ref_len + 1:]], axis=-1)
    return alt


def splice_scores(ref, alt, distance):
    """Pangolin's loss and gain of splice site usage at each position.

    The change of a head is the mean over its replicates of alt minus ref. The loss is the smallest change of the
    heads, the gain the largest.

    Args:
      ref: outputs of the ref sequence in genomic order, array of shape (heads, replicates, positions)
      alt: outputs of the alt sequence in genomic order, array of shape (heads, replicates, alt positions)
      distance: number of positions on each side of the variant

    Returns:
      loss and gain, arrays with the positions of `ref`
    """
    change = (align_alt(alt, ref.shape[-1], distance) - ref).mean(axis=1)
    return change.min(axis=0), change.max(axis=0)


def mask_scores(loss, gain, sites):
    """Pangolin's masking: no gain at annotated splice sites, no loss elsewhere.

    Args:
      loss, gain: arrays of the positions
      sites: indices of the annotated splice sites in the arrays

    Returns:
      the masked loss and gain, as new arrays
    """
    annotated = np.zeros(len(loss), dtype=bool)
    annotated[sites] = True
    return np.where(annotated, loss, np.maximum(loss, 0)), np.where(annotated, np.minimum(gain, 0), gain)


class Pangolin:
    """Pangolin's splice scores of variants: per gene, the largest gain and loss of splice site usage near the variant.

    The scores are those of Pangolin (https://github.com/tkzeng/Pangolin) for the options `-d <distance>` and
    `-m <mask>`. A variant gets a score for each gene that its ref allele overlaps. The models predict the splice site
    usage of the ref and alt sequence on the strand of the gene. With `mask`, gains at the annotated splice sites of
    the gene and losses elsewhere count as 0.

    Every ALT allele of a VCF record is a variant. Pangolin skips a variant if its ref or alt has no A, C, G or T,
    if it is neither an SNV, an MNV, an insertion nor a deletion, if its ref is longer than `2 * distance`, or if its
    ref differs from the FASTA file. Beyond the ends of a chromosome, the sequence is N.

    Args:
      fasta_file: genome FASTA file. A VCF chromosome missing in the FASTA file gains or loses the prefix "chr".
      genes: the genes and their splice sites, as `read_gff3_genes` returns them, with the chromosome names of the
        FASTA file
      models: `PangolinModels`
      distance: number of bases on each side of the variant to score
      mask: mask the scores with the annotated splice sites
      batch_size: number of sequences per model call
    """

    # the columns of the output
    SCHEMA = {
        'variant': pl.String,
        'gene_id': pl.String,
        'gain_score': pl.Float32,
        'gain_pos': pl.Int64,
        'loss_score': pl.Float32,
        'loss_pos': pl.Int64,
        'warnings': pl.List(pl.String),
    }

    def __init__(self, fasta_file, genes, models, distance=50, mask=True, batch_size=8):
        self.fasta = FastaStringExtractor(fasta_file, force_upper=True)
        self.variant_seq_extractor = VariantSeqExtractor(reference_sequence=self.fasta)
        self.chrom_lengths = {chrom: len(record) for chrom, record in self.fasta.fasta.items()}
        self.genes = genes.select('chrom', 'start', 'end', 'strand', 'gene_id')
        self.sites = [np.array(sites, dtype=np.int64) for sites in genes.get_column('sites')]
        self.models = models
        self.distance = distance
        self.mask = mask
        self.batch_size = batch_size

    def _fasta_chrom(self, chrom):
        """The chromosome name of the FASTA file for the VCF chromosome `chrom`, as Pangolin maps it."""
        if chrom not in self.chrom_lengths and 'chr' + chrom in self.chrom_lengths:
            return 'chr' + chrom
        if chrom not in self.chrom_lengths and chrom[3:] in self.chrom_lengths:
            return chrom[3:]
        return chrom

    def _fetch(self, chrom, start, end):
        """The FASTA sequence of the 0-based interval [start, end), with N beyond the chromosome ends."""
        chrom_len = self.chrom_lengths[chrom]
        seq = self.fasta.extract(Interval(chrom, max(start, 0), min(end, chrom_len)))
        return 'N' * max(-start, 0) + seq + 'N' * max(end - chrom_len, 0)

    def _skip_reason(self, chrom, pos, ref, alt):
        """Why Pangolin skips the variant, or None."""
        if not set('ACGT') & set(ref) or not set('ACGT') & set(alt):
            return 'no A, C, G or T in ref or alt'
        if not set(ref + alt) <= set('ACGTN'):
            return 'letters other than A, C, G, T and N in ref or alt'
        if len(ref) != 1 and len(alt) != 1 and len(ref) != len(alt):
            return 'neither an SNV, an MNV, an insertion nor a deletion'
        if len(ref) > 2 * self.distance:
            return f'deletion longer than {2 * self.distance}'
        if chrom not in self.chrom_lengths:
            return f'chromosome {chrom} is not in the FASTA file'
        fasta_ref = self._fetch(chrom, pos - 1, pos - 1 + len(ref))
        if fasta_ref != ref:
            return f'ref differs from the FASTA file ({fasta_ref})'
        return None

    def _sequences(self, chrom, pos, ref, alt):
        """The ref and alt sequence around a variant: `distance` bases plus the flank of the models on each side."""
        flank = self.models.flank + self.distance
        start = pos - 1 - flank
        end = pos - 1 + len(ref) + flank
        ref_seq = self._fetch(chrom, start, end)
        # the alt sequence within the chromosome, and the N beyond its ends from ref_seq
        chrom_start = max(start, 0)
        chrom_end = min(end, self.chrom_lengths[chrom])
        alt_seq = self.variant_seq_extractor.extract(
            Interval(chrom, chrom_start, chrom_end), [Variant(chrom, pos, ref, alt)], anchor=pos - 1, fixed_len=False
        )
        alt_seq = ref_seq[:chrom_start - start] + alt_seq + ref_seq[len(ref_seq) - (end - chrom_end):]
        return ref_seq, alt_seq

    @staticmethod
    def _encode(seq, strand):
        """One-hot encoding of shape (4, length) with N as 0, reverse complemented for the - strand."""
        x = one_hot_dna(seq, neutral_value=0).T
        if strand == '-':
            # reversing the channels A, C, G, T complements the bases
            return x[::-1, ::-1]
        return x

    def _predict(self, seqs):
        """The model outputs of each (sequence, strand) of `seqs`, in genomic order."""
        outputs = [None] * len(seqs)
        by_length = defaultdict(list)
        for i, (seq, strand) in enumerate(seqs):
            by_length[len(seq)].append(i)
        for indices in by_length.values():
            for batch_start in range(0, len(indices), self.batch_size):
                batch = indices[batch_start:batch_start + self.batch_size]
                predictions = self.models.predict(np.stack([self._encode(*seqs[i]) for i in batch]))
                for j, i in enumerate(batch):
                    output = predictions[:, :, j]
                    outputs[i] = output[..., ::-1] if seqs[i][1] == '-' else output
        return outputs

    def _predict_pairs(self, pairs):
        """The rows of `SCHEMA` for the variant-gene pairs of one `SingleVariantMatcher` batch."""
        variants = pairs.group_by('variant_idx', maintain_order=True).agg(
            pl.col('chrom', 'variant_pos', 'variant_ref', 'variant_alt', 'variant_vcf_chrom').first(),
            pl.col('strand').unique(maintain_order=True),
        )
        names = {}
        windows = []
        for idx, chrom, pos, ref, alt, vcf_chrom, strands in variants.iter_rows():
            name = str(Variant(vcf_chrom, pos, ref, alt))
            ref = ref.upper()
            alt = alt.upper()
            reason = self._skip_reason(chrom, pos, ref, alt)
            if reason is not None:
                logger.warning('Skipping variant %s: %s', name, reason)
                continue
            names[idx] = name
            ref_seq, alt_seq = self._sequences(chrom, pos, ref, alt)
            windows.extend(_Window(idx, strand, ref_seq, alt_seq) for strand in strands)

        outputs = self._predict([(seq, w.strand) for w in windows for seq in (w.ref, w.alt)])
        scores = {
            (w.variant_idx, w.strand): splice_scores(outputs[2 * i], outputs[2 * i + 1], self.distance)
            for i, w in enumerate(windows)
        }

        rows = []
        for idx, strand, interval_idx, gene_id, pos in pairs.select(
                'variant_idx', 'strand', 'interval_idx', 'gene_id', 'variant_pos').iter_rows():
            if idx not in names:
                continue
            loss, gain = scores[(idx, strand)]
            warnings = []
            if self.mask:
                # the index of position p is p - (pos - distance)
                sites = self.sites[interval_idx] - (pos - self.distance)
                if len(sites) == 0:
                    warnings.append(NO_SITES_WARNING)
                loss, gain = mask_scores(loss, gain, sites[(sites >= 0) & (sites < len(loss))])
            gain_idx = int(np.argmax(gain))
            loss_idx = int(np.argmin(loss))
            rows.append((names[idx], gene_id, float(gain[gain_idx]), gain_idx - self.distance, float(loss[loss_idx]),
                         loss_idx - self.distance, warnings))
        return pl.DataFrame(rows, schema=self.SCHEMA, orient='row')

    def iter_batches(self, vcf_file, variant_batch_size=256):
        """Score the variants of a VCF file, `variant_batch_size` variants at a time.

        Yields:
          polars DataFrames with the columns of `SCHEMA`, one row per variant and gene, in the order of the VCF file
          and then of the genes. `variant` is chrom:pos:ref>alt with the chromosome name of the VCF file. The
          positions of the gain and loss are relative to the variant position.
        """
        variants = scan_vcf_variants(vcf_file).select('chrom', 'pos', 'ref', 'alt')
        chroms = variants.select(pl.col('chrom').unique()).collect().get_column('chrom')
        variants = variants.with_columns(
            pl.col('chrom').alias('vcf_chrom'),
            pl.col('chrom').replace({chrom: self._fasta_chrom(chrom) for chrom in chroms}),
        )
        matcher = SingleVariantMatcher(variants=variants, intervals=self.genes, interval_attrs=['gene_id'],
                                       variant_batch_size=variant_batch_size)
        for pairs in matcher.iter_batches():
            yield self._predict_pairs(pairs)

    def predict_df(self, vcf_file, variant_batch_size=256):
        """The scores of the variants of a VCF file as one polars DataFrame, see `iter_batches`."""
        return pl.concat([pl.DataFrame(schema=self.SCHEMA), *self.iter_batches(vcf_file, variant_batch_size)])
