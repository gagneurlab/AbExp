# Ported from absplice (https://github.com/gagneurlab/absplice, commit daad7b6), absplice/result.py:
# SplicingOutlierResult with only the AbSplice-DNA prediction from MMSplice and SpliceAI, without samples,
# CADD-Splice, gene TPMs and AbSplice-RNA.
# MIT License, Copyright (c) 2023 Muhammed Hasan Çelik and Nils Wagner; see LICENSE.
from importlib.resources import files

import numpy as np
import pandas as pd

from abexp.absplice.utils import get_abs_max_rows, normalize_gene_annotation, read_csv

_PRECOMPUTED = files('abexp.absplice') / 'precomputed'
GENE_MAP = str(_PRECOMPUTED / 'GENE_MAP.tsv.gz')
ABSPLICE_DNA = str(_PRECOMPUTED / 'AbSplice_DNA.onnx')

dtype_columns = {
    'variant': pd.StringDtype(),
    'gene_id': pd.StringDtype(),
    'tissue': pd.StringDtype(),
    'sample': pd.StringDtype(),
    'Chromosome': pd.StringDtype(),
    'Start': 'Int64',
    'End': 'Int64',
    'Strand': pd.StringDtype(),
    'junction': pd.StringDtype(),
    'event_type': pd.StringDtype(),
    'splice_site': pd.StringDtype(),
    'gene_name': pd.StringDtype(),
    'delta_logit_psi': 'float64',
    'delta_psi': 'float64',
    'ref_psi': 'float64',
    'k': 'Int64',
    'n': 'Int64',
    'median_n': 'float64',
    'novel_junction': pd.BooleanDtype(),
    'weak_site_donor': pd.BooleanDtype(),
    'weak_site_acceptor': pd.BooleanDtype(),
    'delta_score': 'float64',
    'gene_name_spliceai': pd.StringDtype(),
    'gene_tpm': 'float64',
    'tissue_cat': pd.StringDtype(),
    'k_cat': 'Int64',
    'n_cat': 'Int64',
    'median_n_cat': 'float64',
    'psi_cat': 'float64',
    'ref_psi_cat': 'float64',
    'delta_logit_psi_cat': 'float64',
    'delta_psi_cat': 'float64',
    'PHRED': 'float64',
    'AbSplice_DNA': 'float64',
    'AbSplice_RNA': 'float64',
    'pValueGene_g_minus_log10': 'float64',
}


def _load_features_from_model_file(path):
    import onnxruntime
    session = onnxruntime.InferenceSession(path)
    return [i.name for i in session.get_inputs()]


def _predict_onnx(model_path, data):
    import onnxruntime
    session = onnxruntime.InferenceSession(model_path)

    features = session.get_inputs()
    inputs = {
        f.name: np.asarray(data[f.name].values) for f in features
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
        if not isinstance(df, pd.DataFrame):
            df = read_csv(df)
        df = df.reset_index()
        if 'index' in df.columns:
            df = df.drop(columns='index')
        assert pd.Series(columns).isin(df.columns).all()
        return df

    def _validate_dtype(self, df):
        for col in df.columns:
            if col in dtype_columns.keys():
                df = df.astype({col: dtype_columns[col]})
        return df

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
        """
        checks if self.df_spliceai has 'tissue' column.
        If self.df_mmsplice has 'tissue' column and self.df_spliceai does not have 'tissue' column,
        tissue independent spliceai predictions are copied for each tissue in self.df_mmsplice
        """
        df_spliceai = self.df_spliceai
        if self.df_mmsplice is not None:
            l = list()
            for tissue in self.df_mmsplice['tissue'].unique():
                _df = df_spliceai.copy()
                _df['tissue'] = tissue
                l.append(_df)
            self._df_spliceai_tissue = pd.concat(l)
        else:
            self._df_spliceai_tissue = df_spliceai.copy()
            self._df_spliceai_tissue['tissue'] = 'Not provided'
        return self._df_spliceai_tissue

    def _get_maximum_effect(self, df, groupby, score):
        df = df.reset_index()
        if 'index' in df.columns:
            df = df.drop(columns='index')
        if len(set(groupby).difference(df.columns)) != 0:
            raise KeyError(" %s are not in columns" %
                           set(groupby).difference(df.columns))
        return get_abs_max_rows(df.set_index(groupby), groupby, score)

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
            return pd.DataFrame(columns=[*cols_mmsplice, *groupby]).set_index(groupby)

    def _spliceai_agg(self, groupby):
        # SpliceAI
        cols_spliceai = ['delta_score', 'gene_name']
        if self.df_spliceai is not None:
            df_spliceai = self._add_tissue_info_to_spliceai()
            return self._get_maximum_effect(
                df_spliceai, groupby, score='delta_score')
        else:
            return pd.DataFrame(columns=[*cols_spliceai, *groupby]).set_index(groupby)

    @property
    def absplice_dna_input(self):
        """The features of AbSplice-DNA per variant, gene and tissue: the strongest MMSplice and SpliceAI scores."""
        if self._absplice_dna_input is None:
            groupby = ['variant', 'gene_id', 'tissue']

            cols_mmsplice = [
                'junction', 'event_type',
                'splice_site', 'ref_psi', 'median_n',
                'gene_name',
                'delta_logit_psi', 'delta_psi',
            ]
            self._df_mmsplice_agg = self._mmsplice_agg(groupby)

            cols_spliceai = ['delta_score', 'gene_name']
            self._df_spliceai_agg = self._spliceai_agg(groupby)

            # Join MMSplice & SpliceAI
            self._absplice_dna_input = self._df_mmsplice_agg[cols_mmsplice].join(
                self._df_spliceai_agg[cols_spliceai], how='outer', rsuffix='_spliceai')

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
        groupby = ['variant', 'gene_id', 'tissue']

        if self._df_mmsplice_agg is None:
            self._df_mmsplice_agg = self._mmsplice_agg(groupby)

        if self._df_spliceai_agg is None:
            self._df_spliceai_agg = self._spliceai_agg(groupby)

        if 'acceptor_loss_positiin' in self._df_spliceai_agg.columns:
            self._df_spliceai_agg = self._df_spliceai_agg.rename(
                columns={'acceptor_loss_positiin': 'acceptor_loss_position'})

        self._absplice_dna = self._absplice_dna.join(
            self._df_mmsplice_agg[mmsplice_splicemap_cols]).join(
            self._df_spliceai_agg[spliceai_cols])

        return self._absplice_dna

    def _predict_absplice(self, df, absplice_score, model_file, features, abs_features, median_n_cutoff):
        df['splice_site_is_expressed'] = (
                df['median_n'] > median_n_cutoff).astype(int)
        df = df[features].fillna(0)
        if abs_features:
            df = np.abs(df)

        onnx_pred = _predict_onnx(model_file, df)
        df[absplice_score] = np.asarray(onnx_pred, dtype="float32")[:, 1]

        return df

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
          DataFrame with the index variant, gene_id and tissue.
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
