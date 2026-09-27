from typing import Dict, Union
import logging
import polars as pl

log = logging.getLogger(__name__)


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
    df = df.with_columns(pl.struct([
        pl.col(c).alias(rename_fields[c]) if c in rename_fields else pl.col(c)
        for c in cols
    ]).alias(result_name))
    df = df.drop(cols)

    return df


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
            alias = f"{prefix}{sep}{fields}"
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
