import gzip
from dataclasses import dataclass

import polars as pl

# the largest distance in bp between Pangolin's site and a SpliceMap site that matches it
SLACK = 2
# the cap of gain_score for the model, as in AbSplice2: its training data had too few larger gains
GAIN_SCORE_CAP = 0.7
# the inputs of the AbSplice2 model, in its order
FEATURES = [
    'delta_logit_psi',
    'delta_psi',
    'gain_score',
    'loss_score',
    'median_n',
    'median_n_pangolin',
]
# the columns of the output after the key columns chrom, start, end, ref, alt, gene and tissue
OUTPUT_COLUMNS = [
    'gain_score',
    'gain_pos',
    'loss_score',
    'loss_pos',
    'pangolin_tissue_score',
    'ref_psi_pangolin',
    'median_n_pangolin',
    'ref_psi',
    'median_n',
    'delta_logit_psi',
    'delta_psi',
    'junction',
    'event_type',
    'splice_site',
    'AbSplice_DNA',
]

_SPLICEMAP_COLUMNS = {
    'junctions': pl.String,
    'gene_id': pl.String,
    'splice_site': pl.String,
    'ref_psi': pl.Float64,
    'median_n': pl.Float64,
}
_MMSPLICE_COLUMNS = {
    'variant': pl.String,
    'gene_id': pl.String,
    'tissue': pl.String,
    'ref_psi': pl.Float64,
    'median_n': pl.Float64,
    'delta_logit_psi': pl.Float64,
    'delta_psi': pl.Float64,
    'junction': pl.String,
    'event_type': pl.String,
    'splice_site': pl.String,
}


@dataclass(frozen=True)
class _SpliceMap:
    """The rows of one SpliceMap file for some genes, and the tissue in its first line."""
    tissue: str
    df: pl.DataFrame


def read_mmsplice_splicemap(path):
    """The columns of a CSV file of MMSplice with SpliceMaps that `absplice2_dna` needs."""
    return pl.read_csv(path, columns=list(_MMSPLICE_COLUMNS), schema_overrides=_MMSPLICE_COLUMNS)


def _pangolin_scores(pangolin):
    """The rows of abexp.pangolin, as AbSplice2 reads them from the VCF file of Pangolin.

    The gene id loses its version, as in the SpliceMaps. The scores are rounded to 2 decimals like the VCF text of
    Pangolin, which AbSplice2 was trained with. Rounding also turns scores near 0 into 0, which the matching of
    SpliceMap sites tests for. Pangolin rounds in float32, polars in float64. The two differ by 0.01 only for scores
    within float32 precision of a rounding boundary.
    """
    variant = pl.col('variant').str.split(':')
    df = (
        pangolin
        .select(
            'variant',
            variant.list.get(0).alias('chrom'),
            variant.list.get(1).cast(pl.Int64).alias('pos'),
            pl.col('gene_id').str.split('.').list.first(),
            pl.col('gain_score').cast(pl.Float64).round(2),
            'gain_pos',
            pl.col('loss_score').cast(pl.Float64).round(2),
            'loss_pos',
        )
        .unique(maintain_order=True)
    )
    assert not df.select(pl.struct('variant', 'gene_id').is_duplicated().any()).item(), (
        'Pangolin has several scores for one variant and gene'
    )
    return df


def _read_splicemap(path, event_type, gene_ids):
    """The rows of the genes `gene_ids` in a SpliceMap: a gzipped CSV file with the line `# name: <tissue>` first."""
    with gzip.open(path, 'rb') as fd:
        header = fd.readline().decode()
        assert header.startswith('# name: '), f'{path} has no SpliceMap name in its first line'
        tissue = header.split(':')[1].strip()
        df = pl.read_csv(fd, columns=list(_SPLICEMAP_COLUMNS), schema_overrides=_SPLICEMAP_COLUMNS)
    df = (
        df
        .filter(pl.col('gene_id').is_in(gene_ids.implode()))
        .with_columns(
            pl.lit(tissue).alias('tissue'),
            pl.lit(event_type).alias('event_type'),
        )
    )
    return _SpliceMap(tissue=tissue, df=df)


def _match_splicemap_sites(pangolin, sites, score_type):
    """The SpliceMap sites within `SLACK` bp of Pangolin's gain or loss site (`score_type`), per tissue."""
    site = (
        pl.when(pl.col(f'{score_type}_score').abs() > 0)
        .then(pl.col('pos') + pl.col(f'{score_type}_pos'))
        .otherwise(pl.col('pos'))
    )
    return (
        pangolin
        .select(
            'variant',
            'chrom',
            'gene_id',
            pl.int_ranges(site - SLACK, site + SLACK + 1).alias('site_pos'),
        )
        .explode('site_pos')
        .join(sites, on=['chrom', 'gene_id', 'site_pos'], how='inner')
        .select(
            'variant',
            'gene_id',
            'tissue',
            *[
                pl.col(c).alias(f'{c}_{score_type}')
                for c in ['junctions', 'splice_site', 'ref_psi', 'median_n', 'event_type']
            ],
        )
    )


def _pangolin_splicemap(pangolin, splicemap5, splicemap3):
    """Pangolin with the SpliceMap sites of its gain and loss, per tissue of the SpliceMaps.

    Pangolin's site of a gain or loss is the variant position plus its relative position, or the variant position if
    the score is 0. A SpliceMap site of the same gene matches if it is at most `SLACK` bp away; the strand is not
    compared. Pangolin sites without a SpliceMap site are dropped. The gain and loss matches of a variant, gene and
    tissue are joined in all combinations, and each combination takes the SpliceMap site of the larger score:

    - only gain or only loss matched: that site
    - both matched: the gain site if |gain_score| >= |loss_score|, else the loss site

    The score, `ref_psi` and `median_n` of that site are `pangolin_tissue_score`, `ref_psi_pangolin` and
    `median_n_pangolin`. Every variant and gene of Pangolin gets a row per tissue of the SpliceMaps, with nulls where
    no site matched. The tissues are those of all SpliceMaps, also of those without a row for Pangolin's genes.
    """
    gene_ids = pangolin['gene_id'].unique()
    splicemaps = (
        [_read_splicemap(path, 'psi5', gene_ids) for path in splicemap5]
        + [_read_splicemap(path, 'psi3', gene_ids) for path in splicemap3]
    )
    tissues = pl.DataFrame({'tissue': [s.tissue for s in splicemaps]}).unique(maintain_order=True)
    splicemap = pl.concat([s.df for s in splicemaps]).with_columns(pl.col('ref_psi', 'median_n').fill_nan(None))
    del splicemaps

    sites = splicemap.select(
        pl.col('splice_site').str.split(':').list.get(0).alias('chrom'),
        pl.col('splice_site').str.split(':').list.get(1).cast(pl.Int64).alias('site_pos'),
        'gene_id',
        'tissue',
        'junctions',
        'splice_site',
        'ref_psi',
        'median_n',
        'event_type',
    )
    matches = (
        _match_splicemap_sites(pangolin, sites, 'gain')
        .join(_match_splicemap_sites(pangolin, sites, 'loss'), on=['variant', 'gene_id', 'tissue'], how='full',
              coalesce=True)
        .join(pangolin.select('variant', 'gene_id', 'gain_score', 'loss_score'), on=['variant', 'gene_id'],
              how='left')
    )
    use_gain = pl.col('ref_psi_loss').is_null() | (
        pl.col('ref_psi_gain').is_not_null() & (pl.col('gain_score').abs() >= pl.col('loss_score').abs())
    )
    matches = matches.select(
        'variant',
        'gene_id',
        'tissue',
        pl.when(use_gain).then('gain_score').otherwise('loss_score').alias('pangolin_tissue_score'),
        pl.when(use_gain).then('ref_psi_gain').otherwise('ref_psi_loss').alias('ref_psi_pangolin'),
        pl.when(use_gain).then('median_n_gain').otherwise('median_n_loss').alias('median_n_pangolin'),
    )
    return (
        pangolin
        .select('variant', 'gene_id', 'gain_score', 'gain_pos', 'loss_score', 'loss_pos')
        .join(tissues, how='cross')
        .join(matches, on=['variant', 'gene_id', 'tissue'], how='left')
        .unique(maintain_order=True)
    )


def absplice2_dna(pangolin, splicemap5, splicemap3, mmsplice_splicemap, predict, tissue_mapping=None):
    """AbSplice2-DNA per variant, gene and tissue, from Pangolin, MMSplice with SpliceMaps, and the SpliceMaps.

    All rows of MMSplice and of Pangolin are joined in all combinations per variant, gene and tissue. For the model
    only, `gain_score` is capped at `GAIN_SCORE_CAP`, and missing inputs are 0. Per variant, gene and tissue, the
    row with the largest `AbSplice_DNA` stays. The tissue gets its new name from `tissue_mapping` first. Ties go to
    the first row in the order of `OUTPUT_COLUMNS`, so that the output does not depend on the row order. AbSplice2
    itself keeps any one of the tied rows, so its MMSplice columns can differ from these, with the same
    `AbSplice_DNA`.

      pangolin: the output of abexp.pangolin, one row per variant and gene, with unrounded scores
      splicemap5: paths of the psi5 SpliceMaps
      splicemap3: paths of the psi3 SpliceMaps
      mmsplice_splicemap: MMSplice with SpliceMaps, as `read_mmsplice_splicemap` reads it
      predict: a function that gives the AbSplice2-DNA score of each row of a DataFrame with the columns `FEATURES`,
        e.g. `lambda features: model.predict_proba(features.to_pandas())[:, 1]` for the AbSplice2 model
      tissue_mapping: a dict from the tissue names of the SpliceMaps to the output names; other names stay

    Returns the rows that AbSplice2 scores, i.e. the variants and genes that MMSplice or Pangolin scored, and the
    columns of AbSplice2's own output. The variant is split into the key columns `chrom`, `start` (0-based), `end`,
    `ref` and `alt`, and `gene_id` is named `gene`. The columns are the key columns, `tissue` and `OUTPUT_COLUMNS`.
    """
    pangolin = _pangolin_scores(pangolin)
    df = _pangolin_splicemap(pangolin, splicemap5, splicemap3).join(
        mmsplice_splicemap.unique(maintain_order=True), on=['variant', 'gene_id', 'tissue'], how='full',
        coalesce=True,
    )
    del pangolin, mmsplice_splicemap

    features = (
        df
        .select(FEATURES)
        .with_columns(pl.col('gain_score').clip(upper_bound=GAIN_SCORE_CAP))
        .with_columns(pl.all().fill_nan(None).fill_null(0))
    )
    if features.height > 0:
        scores = predict(features)
    else:
        scores = []
    df = df.with_columns(pl.Series('AbSplice_DNA', scores, dtype=pl.Float64))

    tie_break_columns = [c for c in OUTPUT_COLUMNS if c != 'AbSplice_DNA']
    variant = pl.col('variant').str.split(':')
    return (
        df
        .with_columns(pl.col('tissue').replace(tissue_mapping or {}))
        .sort(
            ['variant', 'gene_id', 'tissue', 'AbSplice_DNA', *tie_break_columns],
            descending=[False, False, False, True, *[False] * len(tie_break_columns)],
            nulls_last=True,
        )
        .unique(['variant', 'gene_id', 'tissue'], keep='first', maintain_order=True)
        .select(
            variant.list.get(0).alias('chrom'),
            (variant.list.get(1).cast(pl.Int64) - 1).alias('start'),
            variant.list.get(2).str.split('>').list.get(0).alias('ref'),
            variant.list.get(2).str.split('>').list.get(1).alias('alt'),
            pl.col('gene_id').alias('gene'),
            'tissue',
            *OUTPUT_COLUMNS,
        )
        .with_columns((pl.col('start') + pl.col('ref').str.len_chars()).alias('end'))
        .select('chrom', 'start', 'end', 'ref', 'alt', 'gene', 'tissue', *OUTPUT_COLUMNS)
        .sort(['chrom', 'start', 'end', 'ref', 'alt', 'gene', 'tissue'])
    )
