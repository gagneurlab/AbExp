from typing import Dict, Union, List, Tuple, Iterable
import re

from itertools import chain

import logging

log = logging.getLogger(__name__)

import polars as pl
import polars.datatypes as t


def zip_explode_cols(
        df: Union[pl.DataFrame, pl.LazyFrame],
        cols: list,
        result_name: str,
        rename_fields: Dict[str, str] = None
):
    """
    Explode multiple equally-sized arrays into one struct by zipping all arrays into one `ArrayType[StructType]`

    Args:
        df: The input DataFrame
        cols: The array columns that should be zipped
        result_name: The name of the column that will contain the newly created struct
        rename_fields: dictionary mapping column names to new struct field names.
            Used to rename columns in the newly created struct.

    Returns: `df.with_column(result_name, zip(explode(cols)))`

    """
    if rename_fields is None:
        rename_fields = {}
    
    df = df.explode(cols)
    df = df.with_column(pl.struct([
        pl.col(c).alias(rename_fields[c]) if c in rename_fields else pl.col(c)
        for c in cols
    ]).alias(result_name))
    df = df.drop(cols)

    return df


def normalise_name(raw: str):
    """
    Returns a url-encoded version of a raw string
    """
    from urllib.parse import quote, unquote

    return quote(raw.strip())
    # return re.sub('[^A-Za-z0-9_]+', '.', raw.strip())

    
def denormalise_name(raw: str):
    """
    Returns a url-encoded version of a raw string
    """
    from urllib.parse import quote, unquote

    return unquote(raw.strip())


def __rename_nested_field__(in_field: t.DataType, fieldname_normaliser):
    if isinstance(in_field, t.List):
        inner_type = __rename_nested_field__(in_field.inner, fieldname_normaliser)
        dtype = t.List(inner_type)
    elif isinstance(in_field, t.Struct):
        fields = {}
        for field in in_field.fields:
            fields[fieldname_normaliser(field.name)] = __rename_nested_field__(field.dtype, fieldname_normaliser)
        dtype = t.Struct(fields)
    else:
        dtype = in_field
    return dtype


def normalise_fields_names(df: Union[pl.DataFrame, pl.LazyFrame], fieldname_normaliser=normalise_name):
    """
    Normalize all field names s.t. there are no special characters in the DataFrame schema.
    Uses URL-encoding of special characters by default.
    """
    return df.select([
        pl.col(c).cast(__rename_nested_field__(c_dtype, fieldname_normaliser))
            .alias(fieldname_normaliser(c)) for c, c_dtype in df.schema.items()
    ])

# def flatten(df: pyspark.sql.DataFrame, fieldname_normaliser=normalise_name):
#     """
#     Flatten all fields in dataframe s.t. it can be natively loaded without nested-type support (e.g. Pandas)
#     """
#     cols = []
#     for child in __get_fields_info__(df.schema):
#         if len(child) > 2:
#             ex = "x.{}".format(child[-1])
#             for seg in child[-2:0:-1]:
#                 if seg != '``':
#                     ex = "transform(x.{outer}, x -> {inner})".format(outer=seg, inner=ex)
#             ex = "transform({outer}, x -> {inner})".format(outer=child[0], inner=ex)
#         else:
#             ex = ".".join(child)
#         cols.append(f.expr(ex).alias(fieldname_normaliser("_".join(child).replace('`', ''))))
#     return df.select(cols)


# def rename_values(col, map_dict: dict, default=None):
#     """
#     Renames all values in a column according to `map_dict<old, new>`
#     """
#     if not isinstance(col, pyspark.sql.Column): # Allows either column name string or column instance to be passed
#         col = pl.col(col)
#     mapping_expr = pl.create_map([f.lit(x) for x in chain(*map_dict.items())])
#     if default is None:
#         return mapping_expr.getItem(col)
#     else:
#         return pl.coalesce(mapping_expr.getItem(col), default)


def _recursive_select(fields, c=None, prefix: str = "", sep="."):
    """
    Recursively select fields and yield Tuples of (alias, Polars expression) pairs
    :param fields: nested dictionary/list of columns.
        Example:

        ```python
        {
            "vep": {
                "any": [
                  "transcript_ablation.max",
                  "stop_gained.max",
                ]
            }
        }
        ```
    :param c: struct-type column for which we want to select fields, or None if `fields` is already the leaf
    :param prefix: current prefix of the column names
    :param sep: separator of prefix and column alias
    """
    from collections.abc import Iterable

    if isinstance(fields, str):
        if c is None:
            # 'fields' is leaf column
            alias = fields
            yield alias, pl.col(fields)
        else:
            # we want to select a single column from the 'fields'-struct
            alias=f"{prefix}{sep}{fields}"
            yield alias, c.struct.field(fields).alias(alias)
    elif isinstance(fields, dict):
        # we want to select multiple columns from the 'fields'-struct,
        # with the dictionary key as additional prefix
        for k, v in fields.items():
            if c is None:
                new_c = pl.col(k)
                new_prefix = k
            else:
                new_c = c.struct.field(k)
                new_prefix = f"{prefix}{sep}{k}"

            yield from _recursive_select(v, c=new_c, prefix=new_prefix)
    elif isinstance(fields, Iterable):
        # we want to select multiple columns from the 'fields'-struct
        for v in fields:
            yield from _recursive_select(v, c=c, prefix=prefix)
    else:
        raise ValueError(f"Unknown type: {type(fields)}")


def select_nested_fields(fields, sep="."):
    """
    Recursively select fields and yield Tuples of (alias, Polars expression) pairs
        Example:

        ```python
        select_nested_fields({
            "vep": {
                "any": [
                  "transcript_ablation.max",
                  "stop_gained.max",
                ]
            }
        }, sep=".")
        ```
        Result:
        ```python
        [
            ('vep.any.transcript_ablation.max', Column<'vep[any][transcript_ablation.max] AS `vep.any.transcript_ablation.max`'>),
            ('vep.any.stop_gained.max', Column<'vep[any][stop_gained.max] AS `vep.any.stop_gained.max`'>)
        ]
        ```
    :param fields: nested dictionary/list of columns.
    :param sep: separator of prefix and column alias
    """
    return _recursive_select(fields, sep=sep)

