from typing import Dict, Union, List, Tuple, Iterable
import re

from itertools import chain

import logging

log = logging.getLogger(__name__)

import pyspark.sql.types as t
import pyspark.sql.functions as f
import pyspark


def displayHead(df: pyspark.sql.DataFrame, nrows: int = 5):
    """
    returns the first `nrow` lines of a Spark dataframe as Pandas dataframe

    Args:
        df: The spark dataframe
        nrows: number of rows

    Returns: Pandas dataframe
    """
    return df.limit(nrows).toPandas()


def zip_explode_cols(
        df: pyspark.sql.dataframe.DataFrame,
        cols: list,
        result_name: str,
        rename_fields: Dict[str, str] = None
):
    """
    Explode multiple equally-sized arrays into one struct by zipping all arrays into one `ArrayType[StructType]`

    Args:
        df: The input Spark DataFrame
        cols: The array columns that should be zipped
        result_name: The name of the column that will contain the newly created struct
        rename_fields: dictionary mapping column names to new struct field names.
            Used to rename columns in the newly created struct.

    Returns: `df.withColumn(result_name, zip(explode(cols)))`

    """
    df = df.withColumn(result_name, f.explode(f.arrays_zip(*cols)))

    if rename_fields:  # create schema of new struct by simply renaming the top-level struct fields
        old_schema: t.StructType = df.schema[result_name].dataType

        # rename field if field ist in `old_schema.fieldNames()`
        new_field_names = [
            rename_fields[field] if field in rename_fields else field for field in old_schema.fieldNames()
        ]

        new_schema = t.StructType([
            t.StructField(name, field.dataType) for name, field in zip(new_field_names, old_schema.fields)
        ])

        df = df.withColumn(result_name, f.col(result_name).cast(new_schema))

        # # old method using withColumn and a new struct; breaks with PySpark 3.0
        # df = df.withColumn(target_struct, f.struct(*[
        #     f.col(target_struct + "." + actualName).alias(targetName)
        #     for targetName, actualName in zip(target_colnames, df.schema[target_struct].dataType.fieldNames())
        # ]))

    return df


def __rename_nested_field__(in_field: t.DataType, fieldname_normaliser):
    if isinstance(in_field, t.ArrayType):
        dtype = t.ArrayType(__rename_nested_field__(in_field.elementType, fieldname_normaliser), in_field.containsNull)
    elif isinstance(in_field, t.StructType):
        dtype = t.StructType()
        for field in in_field.fields:
            dtype.add(fieldname_normaliser(field.name), __rename_nested_field__(field.dataType, fieldname_normaliser))
    else:
        dtype = in_field
    return dtype


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


def __get_fields_info__(dtype: t.DataType, name: str = ""):
    ret = []
    if isinstance(dtype, t.StructType):
        for field in dtype.fields:
            for child in __get_fields_info__(field.dataType, field.name):
                wrapped_child = ["{prefix}{suffix}".format(
                    prefix=("" if name == "" else "`{}`.".format(name)), suffix=child[0])] + child[1:]
                ret.append(wrapped_child)
    elif isinstance(dtype, t.ArrayType) and (
            isinstance(dtype.elementType, t.ArrayType) or isinstance(dtype.elementType, t.StructType)):
        for child in __get_fields_info__(dtype.elementType):
            wrapped_child = ["`{}`".format(name)] + child
            ret.append(wrapped_child)
    else:
        return [["`{}`".format(name)]]
    return ret


def normalise_fields_names(df: pyspark.sql.DataFrame, fieldname_normaliser=normalise_name):
    """
    Normalize all field names s.t. there are no special characters in the DataFrame schema.
    Uses URL-encoding of special characters by default.
    """
    return df.select([
        f.col("`{}`".format(field.name)).cast(__rename_nested_field__(field.dataType, fieldname_normaliser))
            .alias(fieldname_normaliser(field.name)) for field in df.schema.fields
    ])


def flatten(df: pyspark.sql.DataFrame, fieldname_normaliser=normalise_name):
    """
    Flatten all fields in Spark dataframe s.t. it can be natively loaded without nested-type support (e.g. Pandas)
    """
    cols = []
    for child in __get_fields_info__(df.schema):
        if len(child) > 2:
            ex = "x.{}".format(child[-1])
            for seg in child[-2:0:-1]:
                if seg != '``':
                    ex = "transform(x.{outer}, x -> {inner})".format(outer=seg, inner=ex)
            ex = "transform({outer}, x -> {inner})".format(outer=child[0], inner=ex)
        else:
            ex = ".".join(child)
        cols.append(f.expr(ex).alias(fieldname_normaliser("_".join(child).replace('`', ''))))
    return df.select(cols)


def rename_values(col, map_dict: dict, default=None):
    """
    Renames all values in a column according to `map_dict<old, new>`
    """
    if not isinstance(col, pyspark.sql.Column): # Allows either column name string or column instance to be passed
        col = f.col(col)
    mapping_expr = f.create_map([f.lit(x) for x in chain(*map_dict.items())])
    if default is None:
        return  mapping_expr.getItem(col)
    else:
        return f.coalesce(mapping_expr.getItem(col), default)


def _recursive_select(fields, c=None, prefix: str = "", sep="."):
    """
    Recursively select fields and return tuple of string alias and pyspark.sql.column.Column
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
            alias=fields
            yield alias, f.col(fields)
        else:
            # we want to select a single column from the 'fields'-struct
            alias=f"{prefix}{sep}{fields}"
            yield alias, c[fields].alias(alias)
    elif isinstance(fields, dict):
        # we want to select multiple columns from the 'fields'-struct,
        # with the dictionary key as additional prefix
        for k, v in fields.items():
            if c is None:
                new_c = f.col(k)
                new_prefix = k
            else:
                new_c = c[k]
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
    Recursively select fields and return tuple of string alias and pyspark.sql.column.Column
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


def melt(
    df: pyspark.sql.DataFrame,
    id_vars: Iterable[str],
    value_vars: Iterable[str],
    var_name: str = "variable",
    value_name: str = "value"
) -> pyspark.sql.DataFrame:
    """
    Convert :class:`DataFrame` from wide to long format.

    :param df: The Dataframe on which the operation will be carried out.
    :param id_vars: Array of columns which will be the index to which the values of the columns will be matched to. 
    :param value_vars: Array of columns that contain the actual values to extract.
    :param var_name: The name of the variable column in the resulting DataFrame.
    :param value_name: The name of the value variable in the resulting DataFrame.
    :returns: Melted DataFrame
    """

    # Create array<struct<variable: str, value: ...>>
    _vars_and_vals = f.array(*(
        f.struct(
            f.lit(c).alias(var_name), 
            f.col(c).alias(value_name)
        ) for c in value_vars
    ))

    # Add to the DataFrame and explode
    _tmp = df.withColumn("_vars_and_vals", f.explode(_vars_and_vals))

    cols = id_vars + [
        f.col("_vars_and_vals")[x].alias(x) for x in [var_name, value_name]
    ]
    return _tmp.select(*cols)


