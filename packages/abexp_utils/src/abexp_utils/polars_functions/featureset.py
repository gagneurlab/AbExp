from typing import Dict, List, Union

from collections import defaultdict
from functools import reduce

import polars as pl
import polars.datatypes as t

from .reshape import select_nested_fields

import logging

log = logging.getLogger(__name__)


__all__ = [
    "transform_featureset",
    "join_featuresets",
]


def transform_featureset(
    fset_df: Union[pl.DataFrame, pl.LazyFrame],
    fset_name: str,
    variables: Dict[str, List[str]],
    index_cols: List[str],
):
    """
    Extracts selected features from some Spark dataframe and flattens them.

    :param fset_df: Some Polars dataframe
    :param fset_name: Name of the resulting DataFrame
    :param variables: (Nested) list/dictionary of columns which will be used as features.
        All feature columns will be flattened and renamed to "feature.{fset_name}@{alias}"
    :param index_cols: List of index columns which should be kept
    :returns: Spark dataframe
    """
    # create column selection from variables
    feature_columns = list(select_nested_fields(variables))
    # rename feature columns: add prefix
    feature_columns = [
        col.alias(f"feature.{fset_name}@{alias}") for alias, col in feature_columns
    ]

    # keep only relevant columns from dataframe
    fset_df = fset_df.select(
        [
            # index
            *[c for c in index_cols if c in fset_df.columns],
            # features
            *feature_columns,
        ]
    )

    return fset_df


def join_featuresets(
    dataframes: Dict[str, Union[pl.DataFrame, pl.LazyFrame]],
    variables: Dict[str, Union[List, Dict, str]],
    index_cols: List[str],
    fill_values: Dict[str, object] = None,
    join="outer_coalesce",
    initial_df: Union[pl.DataFrame, pl.LazyFrame] = None,
    broadcast_columns: Dict[str, Union[pl.DataFrame, pl.LazyFrame]]=None, 
    ignore_missing_columns=True,
):
    """
    Opens multiple parquet dataset, extracts selected features and joins them to a single dataframe.
    Returns a Polars dataframe.

    :param paths: Dictionary of {fset_name: path} pairs to some parquet datasets
    :param variables: Dictionary of (Nested) list/dictionary of columns which will be used as features.
        All feature columns will be flattened and renamed to "feature.{fset_name}@{alias}"
    :param fill_values: Dictionary of column names and fill values to fill missing values
    :param index_cols: List of index columns on which to join
    :param join: Type of the join.
    :param initial_df: dataframe to join on all loaded feature sets.
        Especially useful for left joins.
    :param broadcast_columns: A dictionary of column -> dataframe pairs where each dataframe has
        the distinct values which should be used for broadcasting the column
    :param ignore_missing_columns: Ignore if an index column is missing in all dataframes
    :returns: Polars dataframe
    """
    if broadcast_columns is None:
        broadcast_columns = []

    # first, transform the feature sets and decide for a join order
    partial_column_sets = defaultdict(list)
    full_column_sets = []
    found_columns = set()
    for (fset_name, fset_df) in dataframes.items():
        fset_variables = variables[fset_name]

        fset_df = transform_featureset(
            fset_df=fset_df,
            fset_name=fset_name,
            variables=fset_variables,
            index_cols=index_cols,
        )

        # determine join columns
        join_columns = frozenset({c for c in index_cols if c in fset_df.columns})

        partial_column_sets[join_columns].append(fset_df)

        found_columns = found_columns.union(join_columns)

    if initial_df is not None:
        join_columns = {c for c in index_cols if c in initial_df.columns}
        found_columns = found_columns.union(join_columns)

    if set(index_cols) != found_columns:
        if ignore_missing_columns:
            index_cols = list(found_columns)
        else:
            raise ValueError(
                f"Could not find all columns in the provided dataframes! Missing columns: '{set(index_cols).difference(found_columns)}'"
            )

    # perform partial joins
    joint_partial_column_sets = defaultdict(list)
    for (join_columns, fsets) in partial_column_sets.items():
        joint_fset = reduce(
            lambda left, right: left.join(right, on=list(join_columns), how="outer_coalesce"), fsets
        )
        joint_partial_column_sets[join_columns] = joint_fset

    # compute necessary cross-join columns
    cross_join_columns = {}
    for join_column_set in partial_column_sets.keys():
        cross_join_columns = {
            *cross_join_columns,
            *set(index_cols).difference(join_column_set),
        }

    # distinct_values = dict()
    # for col in cross_join_columns:
    #     # get all datasets that have the column
    #     dfs = []
    #     for idx, df in joint_partial_column_sets.items():
    #         if col in idx:
    #             # get distinct values from df
    #             df_distinct = df.select(pl.col(col)).unique()
    #             dfs.append(df_distinct)
    #
    #     distinct_values[col] = pl.concat(dfs).unique()

    # now perform all cross-joins
    joint_full_column_sets = []
    for idx, df in joint_partial_column_sets.items():
        df_cross_join_columns = set(index_cols).difference(idx)
        for col in df_cross_join_columns:
            if col in broadcast_columns:
                df = broadcast_columns[col].join(df, how="cross")

        joint_full_column_sets.append(df)

    # finally, join all full datasets
    if initial_df is None and join == "outer_coalesce":
        full_df = reduce(
            lambda left, right: left.join(right, on=[
                c for c in index_cols if c in left.columns and c in right.columns
            ], how="outer_coalesce"),
            joint_full_column_sets,
        )
    else:
        full_df = initial_df
        for other_df in joint_full_column_sets:
            full_df = full_df.join(
                other_df,
                on=[c for c in index_cols if c in full_df.columns and c in other_df.columns],
                how=join,
            )

    if fill_values is not None:
        fill_exprs = []
        for col, fval in fill_values.items():
            expr = pl.col(col).fill_null(fval)
            if full_df.schema[col] in t.FLOAT_DTYPES:
                expr = expr.fill_nan(fval)
            expr = expr.alias(col)

            fill_exprs.append(expr)

        full_df = full_df.with_columns(fill_exprs)

    return full_df
