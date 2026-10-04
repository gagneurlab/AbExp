# Ported from mmsplice 2.4.0 (https://github.com/gagneurlab/MMSplice_MTSplice, commit 31513da):
# ExonVariantSeqExtrator, SeqSpliter and ExonSplicingMixin from mmsplice/exon_dataloader.py,
# SplicingVCFMixin from mmsplice/vcf_dataloader.py and the junction VCF dataloaders from
# mmsplice/junction_dataloader.py. kipoiseq2 replaces kipoiseq, kipoi and pyranges.
# MIT License, Copyright (c) 2018, Jun Cheng; see LICENSE.
import logging

import numpy as np
import pandas as pd
import polars as pl
from kipoiseq2 import Interval
from kipoiseq2.extractors import VariantSeqExtractor, SingleVariantMatcher, scan_vcf_variants

from abexp.mmsplice.utils import encodeDNA, region_annotate

logger = logging.getLogger('mmsplice')


class ExonVariantSeqExtrator:
    """
    Extracts sequence with the variant integrated. The lengths overhang
    are fixed irrelevant to variants, even if the variants are indels and
    is in introns, lengths overhang will adapt. If the variant is in the
    exon, the length of the alternative exon (with variant) might change
    for indels.
    """

    def __init__(self, fasta_file):
        self.variant_seq_extractor = VariantSeqExtractor(fasta_file)
        self.fasta = self.variant_seq_extractor.ref_seq_extractor

    def extract(self, interval, variants, overhang=(100, 100)):
        """
        Args:
          interval (kipoiseq2.Interval): zero-based interval of exon
            without overhang.
        """
        down_interval = Interval(
            interval.chrom, interval.start - overhang[0],
            interval.start, strand=interval.strand)
        up_interval = Interval(
            interval.chrom, interval.end,
            interval.end + overhang[1], strand=interval.strand)

        down_seq = self.variant_seq_extractor.extract(
            down_interval, variants, anchor=interval.start)
        up_seq = self.variant_seq_extractor.extract(
            up_interval, variants, anchor=interval.start)

        exon_seq = self.variant_seq_extractor.extract(
            interval, variants, anchor=0, fixed_len=False)

        if interval.strand == '-':
            down_seq, up_seq = up_seq, down_seq

        return down_seq + exon_seq + up_seq


class SeqSpliter:
    """
    Splits given seq for each modules.

    length arguments of the __init__ function refer to the prefered
    sequence  length of the models

    Args:
      exon_cut_l: number of bp to cut out at the begining of an exon
      exon_cut_r: number of bp to cut out at the end of an exon
        (cut out the part that is considered as acceptor site or donor site)
      acceptor_intron_cut: number of bp to cut out at the end of
        acceptor intron that consider as acceptor site
      donor_intron_cut: number of bp to cut out at the end of donor intron
        that consider as donor site
      acceptor_intron_len: length in acceptor intron to consider
        for acceptor site model
      acceptor_exon_len: length in acceptor exon to consider
        for acceptor site model
      donor_intron_len: length in donor intron to consider for donor site model
      donor_exon_len: length in donor exon to consider for donor site model
    """

    def __init__(self, exon_cut_l=0, exon_cut_r=0,
                 acceptor_intron_cut=6, donor_intron_cut=6,
                 acceptor_intron_len=50, acceptor_exon_len=3,
                 donor_exon_len=5, donor_intron_len=13,
                 pattern_warning=False):
        self.exon_cut_l = exon_cut_l
        self.exon_cut_r = exon_cut_r
        self.acceptor_intron_cut = acceptor_intron_cut
        self.donor_intron_cut = donor_intron_cut
        self.acceptor_intron_len = acceptor_intron_len
        self.acceptor_exon_len = acceptor_exon_len
        self.donor_exon_len = donor_exon_len
        self.donor_intron_len = donor_intron_len
        self.pattern_warning = pattern_warning

    def split(self, seq, overhang, exon_row='', pattern_warning=True):
        """
        Split seqeunce for each module.

        Args:
          seq: seqeunce to split.
          overhang: (intron_length acceptor side, intron_length donor side) of
                    the input sequence
        """
        pattern_warning = self.pattern_warning and pattern_warning

        intronl_len, intronr_len = overhang
        assert intronl_len <= len(seq), "Input sequence acceptor intron" \
            " length cannot be longer than the input sequence"
        assert intronr_len <= len(seq), "Input sequence donor intron length" \
            " cannot be longer than the input sequence"

        # need to pad N if left seq not enough long
        lackl = self.acceptor_intron_len - intronl_len
        if lackl >= 0:
            seq = "N" * (lackl + 1) + seq
            intronl_len += lackl + 1
        lackr = self.donor_intron_len - intronr_len
        if lackr >= 0:
            seq = seq + "N" * (lackr + 1)
            intronr_len += lackr + 1

        acceptor_intron = seq[:intronl_len - self.acceptor_intron_cut]

        acceptor_start = intronl_len - self.acceptor_intron_len
        acceptor_end = intronl_len + self.acceptor_exon_len
        acceptor = seq[acceptor_start: acceptor_end]

        exon_start = intronl_len + self.exon_cut_l
        exon_end = -intronr_len - self.exon_cut_r
        exon = seq[exon_start: exon_end]

        donor_start = -intronr_len - self.donor_exon_len
        donor_end = -intronr_len + self.donor_intron_len
        donor = seq[donor_start: donor_end]

        donor_intron = seq[-intronr_len + self.donor_intron_cut:]

        if not exon:
            exon = 'N'

        if pattern_warning:
            if donor[self.donor_exon_len:self.donor_exon_len + 2] != "GT" \
               and overhang[1]:
                logger.warning('None GT donor: %s' % str(exon_row))

            if acceptor[self.acceptor_intron_len - 2:self.acceptor_intron_len] != "AG" \
               and overhang[0]:
                logger.warning('None AG acceptor: %s' % str(exon_row))

        splits = {
            "acceptor_intron": acceptor_intron,
            "acceptor": acceptor,
            "exon": exon,
            "donor": donor,
            "donor_intron": donor_intron
        }

        return splits


class ExonSplicingMixin:
    """
    Builds the MMSplice inputs of variant-exon pairs.

    Args:
      fasta_file: fasta file to fetch exon sequences.
      split_seq: whether or not already split the sequence
        when loading the data.
      encode: if split sequence, should it be one-hot-encoded.
      overhang: overhang of exon to fetch flanking sequence of exon.
      seq_spliter: SeqSpliter class instance specific how to split seqs.
    """
    optional_metadata = ('exon_id', 'gene_id', 'gene_name',
                         'transcript_id', 'junction', 'side', 'region')

    def __init__(self, fasta_file, split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None):
        self.fasta_file = fasta_file
        self.split_seq = split_seq
        self.encode = encode
        self.overhang = overhang
        self.spliter = seq_spliter or SeqSpliter()
        self.vseq_extractor = ExonVariantSeqExtrator(fasta_file)
        self.fasta = self.vseq_extractor.fasta

    def _next(self, exon, variant, overhang=None, mask_module=None):
        overhang = overhang or self.overhang

        inputs = {
            'seq': self.fasta.extract(Interval(
                exon.chrom, exon.start - overhang[0],
                exon.end + overhang[1], strand=exon.strand), use_strand=True).upper(),
            'mut_seq': self.vseq_extractor.extract(
                exon, [variant], overhang=overhang).upper()
        }

        if exon.strand == '-':
            overhang = (overhang[1], overhang[0])

        if self.split_seq:
            inputs['seq'] = self.spliter.split(inputs['seq'], overhang, exon)
            inputs['mut_seq'] = self.spliter.split(inputs['mut_seq'], overhang,
                                                   exon, pattern_warning=False)
            if mask_module:
                for i in mask_module:
                    if i in inputs['seq']:
                        inputs['seq'][i] = 'N' * len(inputs['seq'][i])
                        inputs['mut_seq'][i] = 'N' * len(inputs['mut_seq'][i])
                    else:
                        raise ValueError('%s is not in mmsplice modules' % i)

            if self.encode:
                inputs = {k: self._encode_seq(v) for k, v in inputs.items()}

        return {
            'inputs': inputs,
            'metadata': {
                'variant': self._variant_to_dict(variant, exon),
                'exon': self._exon_to_dict(exon, overhang)
            }
        }

    def _encode_batch_seq(self, batch):
        return {k: encodeDNA(v.tolist()) for k, v in batch.items()}

    def _encode_seq(self, seq):
        return {k: encodeDNA([v]) for k, v in seq.items()}

    def _variant_to_dict(self, variant, exon):
        return {
            'chrom': variant.chrom,
            'pos': variant.pos,
            'ref': variant.ref,
            'alt': variant.alt,
            'annotation': str(variant),
            'region': region_annotate(variant, exon)
        }

    def _exon_to_dict(self, exon, overhang):
        return {
            'chrom': exon.chrom,
            'start': exon.start,
            'end': exon.end,
            'strand': exon.strand,
            'left_overhang': overhang[0],
            'right_overhang': overhang[1],
            'annotation': str(exon),
            **exon.attrs
        }


def _remove_chr(df):
    df = df.copy()
    df['Chromosome'] = df['Chromosome'].str.replace('chr', '')
    return df


def _add_chr(df):
    df = df.copy()
    df['Chromosome'] = 'chr' + df['Chromosome'].astype(str)
    return df


class SplicingVCFMixin(ExonSplicingMixin):
    """
    Matches the variants of a VCF file with exons.

    Args:
      exons: pandas DataFrame with the exons (with overhang), with the
        columns Chromosome, Start (0-based), End, Strand and the
        `interval_attrs` columns.
    """

    def __init__(self, exons, annotation, fasta_file, vcf_file,
                 split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None,
                 interval_attrs=tuple()):
        super().__init__(fasta_file, split_seq, encode, overhang, seq_spliter)
        self.exons = exons
        self.annotation = annotation
        self.vcf_file = vcf_file
        self.vcf_chroms = set(
            scan_vcf_variants(vcf_file).select('chrom').unique().collect().get_column('chrom'))
        self._check_chrom_annotation()
        intervals = pl.from_pandas(
            self.exons[['Chromosome', 'Start', 'End', 'Strand', *interval_attrs]].reset_index(drop=True)
        ).rename({'Chromosome': 'chrom', 'Start': 'start', 'End': 'end', 'Strand': 'strand'})
        self.matcher = SingleVariantMatcher(
            vcf_file=vcf_file, intervals=intervals,
            interval_attrs=list(interval_attrs)
        )
        self._generator = iter(self.matcher)

    def _check_chrom_annotation(self):
        vcf_chroms = self.vcf_chroms
        if not vcf_chroms:
            # a VCF file without variants matches nothing
            return

        fasta_chroms = set(self.fasta.fasta.keys())
        if not fasta_chroms.intersection(vcf_chroms):
            raise ValueError(
                'Fasta chrom names do not match with vcf chrom names')

        gtf_chroms = set(self.exons['Chromosome'])
        if not gtf_chroms.intersection(vcf_chroms):
            chr_annotaion = any(chrom.startswith('chr')
                                for chrom in vcf_chroms)
            if not chr_annotaion:
                self.exons = _remove_chr(self.exons)
            else:
                self.exons = _add_chr(self.exons)

        gtf_chroms = set(self.exons['Chromosome'])
        if not gtf_chroms.intersection(vcf_chroms):
            raise ValueError(
                'GTF chrom names do not match with vcf chrom names')


class _JunctionVCFDataloader(SplicingVCFMixin):
    """
    Load intron annotation (bed) file along with a vcf file,
      return reference sequence and alternative sequence.

    Args:
      intron_annotation: path of tabular file or pandas DataFrame with the
        intron (junction) annotation, with the colunms
        `'Chromosome', 'Start', 'End', 'Strand'`. (0-based)
      fasta_file: file path; Genome sequence
      vcf_file: vcf file, each line should contain one
        and only one variant, left-normalized
      split_seq: whether or not already split the sequence
        when loading the data. Otherwise it can be done in the model class.
      encode: if split sequence, should it be one-hot-encoded.
      overhang: overhang of exon to fetch flanking sequence of exon.
      seq_spliter: SeqSpliter class instance specific how to split seqs.
         if None, use the default arguments of SeqSpliter
    """

    def __init__(self, intron_annotation, fasta_file, vcf_file,
                 event_type, split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None, exon_len=100):
        self.event_type = event_type
        exons = self._read_junction(intron_annotation, event_type,
                                    overhang, exon_len)
        super().__init__(exons, intron_annotation, fasta_file, vcf_file,
                         split_seq, encode, overhang, seq_spliter,
                         interval_attrs=('junction',))

    @staticmethod
    def _read_junction(intron_annotation, event_type, overhang=(100, 100), exon_len=100):
        if type(intron_annotation) == str:
            df = pd.read_csv(intron_annotation, dtype={'Chromosome': str})
        else:
            df = intron_annotation
            df['Chromosome'] = df['Chromosome'].astype(str)

        df = pd.concat([df, _JunctionVCFDataloader._junction_to_exon(
            df, overhang, exon_len)], axis=1)

        # split donor-acceptor from horizontal dataframe to two df
        # and merge this two df as one vertical dataframes
        df = df.rename(columns={"Start": "Junction_Start",
                                "End": "Junction_End"})

        if event_type == 'psi5':
            df_exons = df[['Chromosome', 'Acceptor_Start', 'Acceptor_End',
                           'Strand', 'Junction_Start', 'Junction_End']] \
                .rename(columns={'Acceptor_Start': 'Start', 'Acceptor_End': 'End'})
        elif event_type == 'psi3':
            df_exons = df[['Chromosome', 'Donor_Start', 'Donor_End',
                           'Strand', 'Junction_Start', 'Junction_End']] \
                .rename(columns={'Donor_Start': 'Start', 'Donor_End': 'End'})
        else:
            raise ValueError('event_type should be "psi5" or "psi3"')

        df_exons['junction'] = junction_df_junction_str(df_exons)

        del df_exons['Junction_Start']
        del df_exons['Junction_End']
        return df_exons

    @staticmethod
    def _junction_to_exon(df, overhang=(100, 100), exon_len=100):
        # calculates donor-acceptor exon start end based
        # on given fixed exon lenght parameter and overhang (1-based)
        # --...-- represents junction (- exon), (. intron), (/ cut)
        mat = np.where(df['Strand'] == '-',
                       # --a./..d--/
                       (df['End'] - overhang[0],  df['End'] + exon_len,
                        # /--a../.d--
                        df['Start'] - exon_len, df['Start'] + overhang[1]),
                       # /--d../.a--
                       (df['Start'] - exon_len, df['Start'] + overhang[1],
                        # --d./..a--/
                        df['End'] - overhang[0], df['End'] + exon_len)
                       )
        return pd.DataFrame(
            mat.T,
            columns=['Donor_Start', 'Donor_End',
                     'Acceptor_Start', 'Acceptor_End'],
            index=df.index
        )

    def __next__(self):
        exon, variant = next(self._generator)

        # TODO: fix overhang will be reversed based on strand
        #   in `ExonSplicingMixin`
        # ---...*---
        if (self.event_type == 'psi3' and exon.strand == '-') \
           or (self.event_type == 'psi5' and exon.strand == '+'):
            overhang = (self.overhang[0], 0)
        # ---*...---
        elif (self.event_type == 'psi5' and exon.strand == '-') \
                or (self.event_type == 'psi3' and exon.strand == '+'):
            overhang = (0, self.overhang[1])

        # the exon without overhang
        exon = Interval(exon.chrom, exon.start + overhang[0], exon.end - overhang[1],
                        strand=exon.strand, attrs=exon.attrs)

        if self.event_type == 'psi5':
            mask = ['donor', 'donor_intron']
        else:
            mask = ['acceptor', 'acceptor_intron']

        row = self._next(exon, variant, overhang, mask)

        return row

    def __iter__(self):
        return self


def junction_df_junction_str(df):
    return df['Chromosome'] + \
        ':' + df['Junction_Start'].astype('str') + \
        '-' + df['Junction_End'].astype('str') + \
        ':' + df['Strand']


class JunctionPSI5VCFDataloader(_JunctionVCFDataloader):
    def __init__(self, intron_annotation, fasta_file, vcf_file,
                 split_seq=True, encode=True, overhang=(100, 100),
                 seq_spliter=None, exon_len=100, **kwargs):
        super().__init__(intron_annotation, fasta_file, vcf_file,
                         'psi5', split_seq=split_seq, encode=encode,
                         overhang=overhang, seq_spliter=seq_spliter,
                         exon_len=exon_len, **kwargs)


class JunctionPSI3VCFDataloader(_JunctionVCFDataloader):
    def __init__(self, intron_annotation, fasta_file, vcf_file,
                 split_seq=True, encode=True, overhang=(100, 100),
                 seq_spliter=None, exon_len=100, **kwargs):
        super().__init__(intron_annotation, fasta_file, vcf_file, 'psi3',
                         split_seq=split_seq, encode=encode,
                         overhang=overhang, seq_spliter=seq_spliter,
                         exon_len=exon_len, **kwargs)
