# Ported from absplice (https://github.com/gagneurlab/absplice, commit daad7b6), absplice/result.py:
# SplicingOutlierResult with only the AbSplice-DNA prediction from MMSplice and SpliceAI, without samples,
# CADD-Splice, gene TPMs and AbSplice-RNA.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
from importlib.resources import files

import numpy as np
import polars as pl

from abexp.absplice.utils import get_abs_max_rows, normalize_gene_annotation, read_csv

_PRECOMPUTED = files('abexp.absplice') / 'precomputed'
GENE_MAP = str(_PRECOMPUTED / 'GENE_MAP.tsv.gz')
ABSPLICE_DNA = str(_PRECOMPUTED / 'AbSplice_DNA.onnx')

dtype_columns = {
    'variant': pl.String,
    'gene_id': pl.String,
    'tissue': pl.String,
    'sample': pl.String,
    'Chromosome': pl.String,
    'Start': pl.Int64,
    'End': pl.Int64,
    'Strand': pl.String,
    'junction': pl.String,
    'event_type': pl.String,
    'splice_site': pl.String,
    'gene_name': pl.String,
    'delta_logit_psi': pl.Float64,
    'delta_psi': pl.Float64,
    'ref_psi': pl.Float64,
    'k': pl.Int64,
    'n': pl.Int64,
    'median_n': pl.Float64,
    'novel_junction': pl.Boolean,
    'weak_site_donor': pl.Boolean,
    'weak_site_acceptor': pl.Boolean,
    'delta_score': pl.Float64,
    'gene_name_spliceai': pl.String,
    'gene_tpm': pl.Float64,
    'tissue_cat': pl.String,
    'k_cat': pl.Int64,
    'n_cat': pl.Int64,
    'median_n_cat': pl.Float64,
    'psi_cat': pl.Float64,
    'ref_psi_cat': pl.Float64,
    'delta_logit_psi_cat': pl.Float64,
    'delta_psi_cat': pl.Float64,
    'PHRED': pl.Float64,
    'AbSplice_DNA': pl.Float64,
    'AbSplice_RNA': pl.Float64,
    'pValueGene_g_minus_log10': pl.Float64,
}

# the columns of the rows of AbSplice-DNA
GROUPBY = ['variant', 'gene_id', 'tissue']


def _load_features_from_model_file(path):
    import onnxruntime
    session = onnxruntime.InferenceSession(path)
    return [i.name for i in session.get_inputs()]


def _predict_onnx(model_path, data):
    import onnxruntime
    session = onnxruntime.InferenceSession(model_path)

    features = session.get_inputs()
    inputs = {
        f.name: data[f.name].to_numpy() for f in features
    }

    results = session.run(None, inputs)[0]
    return results


class SplicingOutlierResult:
    """Combines the MMSplice and SpliceAI scores of variants and predicts AbSplice-DNA.

    Args:
      df_mmsplice: MMSplice scores per variant, junction and tissue, e.g. from `SpliceOutlier.predict_save`.
      df_spliceai: SpliceAI scores per variant and gene name.
      gene_map: table with the columns gene_id and gene_name, which maps the SpliceAI gene names to gene IDs.
        Default is GENE_MAP.

    Each table is a polars DataFrame or the path of a CSV, TSV or parquet file.
    """

    def __init__(
            self,
            df_mmsplice=None,
            df_spliceai=None,
            gene_map=None,
    ):
        self.df_mmsplice = self.validate_df_mmsplice(df_mmsplice)
        self.gene_map = self.validate_df_gene_map(gene_map)
        self.df_spliceai = self.validate_df_spliceai(df_spliceai)
        self._absplice_dna_input = None
        self._absplice_dna = None
        self._df_spliceai_tissue = None
        self._df_mmsplice_agg = None
        self._df_spliceai_agg = None

    def _validate_df(self, df, columns):
        if not isinstance(df, pl.DataFrame):
            df = read_csv(df)
        assert set(columns).issubset(df.columns)
        return df

    def _validate_dtype(self, df):
        return df.cast({col: dtype for col, dtype in dtype_columns.items() if col in df.columns})

    def validate_df_mmsplice(self, df_mmsplice):
        if df_mmsplice is not None:
            df_mmsplice = self._validate_df(
                df_mmsplice,
                columns=['variant', 'gene_id', 'tissue',
                         'delta_psi', 'ref_psi', 'median_n'])
            df_mmsplice = self._validate_dtype(df_mmsplice)
        return df_mmsplice

    def validate_df_spliceai(self, df_spliceai):
        if df_spliceai is not None:
            df_spliceai = self._validate_df(
                df_spliceai,
                columns=['variant', 'gene_name', 'delta_score'])
            if self.gene_map is not None:
                df_spliceai = normalize_gene_annotation(
                    df_spliceai, self.gene_map, key='gene_name', value='gene_id')
            df_spliceai = self._validate_dtype(df_spliceai)
        return df_spliceai

    def validate_df_gene_map(self, gene_map):
        if gene_map is not None:
            gene_map = self._validate_df(
                gene_map,
                columns=['gene_id', 'gene_name'])
            gene_map = self._validate_dtype(gene_map)
        else:
            gene_map = self._validate_df(
                GENE_MAP,
                columns=['gene_id', 'gene_name'])
        return gene_map

    def _add_tissue_info_to_spliceai(self):
        """The SpliceAI scores, which do not depend on the tissue, for each tissue of `df_mmsplice`.

        Without `df_mmsplice`, the tissue is 'Not provided'.
        """
        if self.df_mmsplice is not None:
            tissues = self.df_mmsplice.select(pl.col('tissue').unique(maintain_order=True))
        else:
            tissues = pl.DataFrame({'tissue': ['Not provided']})
        self._df_spliceai_tissue = self.df_spliceai.drop('tissue', strict=False).join(tissues, how='cross')
        return self._df_spliceai_tissue

    def _get_maximum_effect(self, df, groupby, score):
        if len(set(groupby).difference(df.columns)) != 0:
            raise KeyError(" %s are not in columns" %
                           set(groupby).difference(df.columns))
        return get_abs_max_rows(df, groupby, score)

    @staticmethod
    def _empty(columns):
        return pl.DataFrame(schema={col: dtype_columns[col] for col in columns})

    def _mmsplice_agg(self, groupby):
        # MMSplice (SpliceMap)
        cols_mmsplice = [
            'junction', 'event_type',
            'splice_site', 'ref_psi', 'median_n',
            'gene_name',
            'delta_logit_psi', 'delta_psi',
        ]
        if self.df_mmsplice is not None:
            return self._get_maximum_effect(
                self.df_mmsplice, groupby, score='delta_psi')
        else:
            return self._empty([*groupby, *cols_mmsplice])

    def _spliceai_agg(self, groupby):
        # SpliceAI
        cols_spliceai = ['delta_score', 'gene_name']
        if self.df_spliceai is not None:
            df_spliceai = self._add_tissue_info_to_spliceai()
            return self._get_maximum_effect(
                df_spliceai, groupby, score='delta_score')
        else:
            return self._empty([*groupby, *cols_spliceai])

    @staticmethod
    def _join(df, other, how):
        # A missing gene_id, of a SpliceAI gene name without an entry in the gene map, matches a missing gene_id,
        # as in the pandas index joins of absplice daad7b6.
        return df.join(other, on=GROUPBY, how=how, coalesce=True, suffix='_spliceai', nulls_equal=True,
                       maintain_order='left')

    @property
    def absplice_dna_input(self):
        """The features of AbSplice-DNA per variant, gene and tissue: the strongest MMSplice and SpliceAI scores.

        The rows are sorted by variant, gene_id and tissue.
        """
        if self._absplice_dna_input is None:
            cols_mmsplice = [
                'junction', 'event_type',
                'splice_site', 'ref_psi', 'median_n',
                'gene_name',
                'delta_logit_psi', 'delta_psi',
            ]
            self._df_mmsplice_agg = self._mmsplice_agg(GROUPBY)

            cols_spliceai = ['delta_score', 'gene_name']
            self._df_spliceai_agg = self._spliceai_agg(GROUPBY)

            # Join MMSplice & SpliceAI
            self._absplice_dna_input = self._join(
                self._df_mmsplice_agg.select(*GROUPBY, *cols_mmsplice),
                self._df_spliceai_agg.select(*GROUPBY, *cols_spliceai), how='full',
            ).sort(GROUPBY, nulls_last=True)

        return self._absplice_dna_input

    def add_extra_info(self):
        mmsplice_splicemap_cols = [
            'junction',
            'event_type',
            'splice_site',
            'ref_psi',
            'median_n'
        ]
        spliceai_cols = [
            'acceptor_gain',
            'acceptor_loss',
            'donor_gain',
            'donor_loss',
            'acceptor_gain_position',
            'acceptor_loss_position',
            'donor_gain_position',
            'donor_loss_position'
        ]

        # get aggregated scores of SpliceAI and MMSplice + SpliceMap
        if self._df_mmsplice_agg is None:
            self._df_mmsplice_agg = self._mmsplice_agg(GROUPBY)

        if self._df_spliceai_agg is None:
            self._df_spliceai_agg = self._spliceai_agg(GROUPBY)

        if 'acceptor_loss_positiin' in self._df_spliceai_agg.columns:
            self._df_spliceai_agg = self._df_spliceai_agg.rename(
                {'acceptor_loss_positiin': 'acceptor_loss_position'})

        self._absplice_dna = self._join(
            self._join(self._absplice_dna, self._df_mmsplice_agg.select(*GROUPBY, *mmsplice_splicemap_cols), 'left'),
            self._df_spliceai_agg.select(*GROUPBY, *spliceai_cols), 'left')

        return self._absplice_dna

    def _predict_absplice(self, df, absplice_score, model_file, features, abs_features, median_n_cutoff):
        expressed = (pl.col('median_n') > median_n_cutoff).fill_null(False)
        df = df.with_columns(splice_site_is_expressed=expressed.cast(pl.Int64))
        df = df.select(*GROUPBY, pl.col(features).fill_null(0))
        if abs_features:
            df = df.with_columns(pl.col(features).abs())

        onnx_pred = _predict_onnx(model_file, df)
        return df.with_columns(pl.Series(absplice_score, np.asarray(onnx_pred, dtype="float32")[:, 1]))

    def predict_absplice_dna(
            self,
            model_file=ABSPLICE_DNA,
            features=None,
            abs_features=False,
            median_n_cutoff=10,
            extra_info=True
    ):
        """Predict AbSplice-DNA per variant, gene and tissue.

        Args:
          model_file: ONNX file of the AbSplice-DNA model.
          features: input columns of the model. Default is the sorted inputs of the ONNX model.
          abs_features: use the absolute values of the features.
          median_n_cutoff: a splice site with a median coverage above this counts as expressed.
          extra_info: add the MMSplice and SpliceAI details, e.g. the junction and the SpliceAI positions.

        Returns:
          polars DataFrame with the columns variant, gene_id and tissue, the features and AbSplice_DNA, sorted by
          variant, gene_id and tissue.
        """
        # Load model and extract features
        if features is None:
            features = sorted(_load_features_from_model_file(model_file))

        self._absplice_dna = self._predict_absplice(
            df=self.absplice_dna_input,
            absplice_score='AbSplice_DNA',
            model_file=model_file,
            features=features,
            abs_features=abs_features,
            median_n_cutoff=median_n_cutoff)

        if extra_info:
            self._absplice_dna = self.add_extra_info()

        return self._absplice_dna
