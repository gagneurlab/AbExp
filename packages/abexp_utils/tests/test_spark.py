import abexp_utils.spark_functions as sp_funcs
from pyspark.sql import SparkSession, Row
from pyspark.sql.types import StructType, StructField, IntegerType, MapType, StringType
import pyspark.sql.functions as f
import pandas as pd

spark = SparkSession.builder.master("local").appName("example").getOrCreate()


def test_displayHead():
    # Create Sample DataFrame
    data = [("Alice", 1), ("Bob", 2), ("Cathy", 3)]
    df = spark.createDataFrame(data, ["name", "id"])

    # Define Expected Output
    expected_df = pd.DataFrame({"name": ["Alice", "Bob"], "id": [1, 2]})

    # Apply Function
    result_df = sp_funcs.displayHead(df, nrows=2)

    # Assert Equality
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_flatten():
    # Create Sample DataFrame
    df = spark.createDataFrame(
        data=[
            (1, {"nested": {"field1": 10, "field2": 20, "field3": 30}}),
            (2, {"nested": {"field1": 30, "field2": 40, "field3": 50}})
        ],
        schema=StructType([
            StructField("id", IntegerType(), True),
            StructField("data", StructType([
                StructField("nested", StructType([
                    StructField("field1", IntegerType(), True),
                    StructField("field2", IntegerType(), True),
                    StructField("field3", IntegerType(), True)
                ]), True)
            ]), True)
        ])
    )

    # Define Expected Output
    expected_df = spark.createDataFrame(
        data=[
            (1, 10, 20, 30),
            (2, 30, 40, 50)
        ],
        schema=StructType([
            StructField("id", IntegerType(), True),
            StructField("data_nested_field1", IntegerType(), True),
            StructField("data_nested_field2", IntegerType(), True),
            StructField("data_nested_field3", IntegerType(), True)
        ])
    )

    # Apply Function
    result_df = sp_funcs.flatten(df)

    # Assert Equality
    assert result_df.collect() == expected_df.collect(), f"Expected:\n{expected_df}\nGot:\n{result_df}"


def test_zip_explode_cols():
    # Create Sample DataFrame
    data = [(1, [4, 5], [7, 8]), (2, [6, 7], [9, 10])]
    df = spark.createDataFrame(data, ["A", "B", "C"])

    # Define Expected Output
    expected_data = [
        (1, Row(B=4, C=7)),
        (1, Row(B=5, C=8)),
        (2, Row(B=6, C=9)),
        (2, Row(B=7, C=10))
    ]
    expected_df = spark.createDataFrame(expected_data, ["A", "new_name"])

    # Apply Function
    result_df = sp_funcs.zip_explode_cols(df, cols=["B", "C"], result_name="new_name")

    # Assert Equality
    assert result_df.collect() == expected_df.collect(), f"Expected:\n{expected_df.collect()}\nGot:\n{result_df.collect()}"


def test_normalise_field_names():
    # Define df
    df = spark.createDataFrame(data=[
        (1, {"C D": 3, "E&F": 4}, {"H I": 7, "J&K": 8}),
        (2, {"C D": 5, "E&F": 6}, {"H I": 9, "J&K": 10})
    ], schema=StructType([
        StructField("A 1", IntegerType(), True),
        StructField("B 2", StructType([
            StructField("C D", IntegerType(), True),
            StructField("E&F", IntegerType(), True)
        ]), True),
        StructField("G 3", StructType([
            StructField("H I", IntegerType(), True),
            StructField("J&K", IntegerType(), True)
        ]), True),
    ]))

    expected_df = spark.createDataFrame(data=[
        (1, {"C%20D": 3, "E%26F": 4}, {"H%20I": 7, "J%26K": 8}),
        (2, {"C%20D": 5, "E%26F": 6}, {"H%20I": 9, "J%26K": 10})
    ], schema=StructType([
        StructField("A%201", IntegerType(), True),
        StructField("B%202", StructType([
            StructField("C%20D", IntegerType(), True),
            StructField("E%26F", IntegerType(), True)
        ]), True),
        StructField("G%203", StructType([
            StructField("H%20I", IntegerType(), True),
            StructField("J%26K", IntegerType(), True)
        ]), True),
    ]))

    # Apply Function
    result_df = sp_funcs.normalise_fields_names(df)

    # Assert Equality
    assert result_df.collect() == expected_df.collect(), f"Expected:\n{expected_df}\nGot:\n{result_df}"


def test_rename_values():
    # Create Sample DataFrame
    df = spark.createDataFrame(
        data=[("Alice", 1), ("Bob", 2), ("Cathy", 3)],
        schema=StructType([
            StructField("name", StringType(), True),
            StructField("id", IntegerType(), True)
        ]))

    # Define Mapping Dictionary
    map_dict = {"Alice": "Alicia", "Bob": "Robert", "Cathy": "Cathy"}

    # Define Expected Output
    expected_df = spark.createDataFrame(
        data=[("Alicia", 1), ("Robert", 2), ("Cathy", 3)],
        schema=StructType([
            StructField("name", StringType(), True),
            StructField("id", IntegerType(), True)
        ]))

    # Apply Function
    result_df = df.withColumn("name", sp_funcs.rename_values("name", map_dict))

    # Assert Equality
    assert result_df.collect() == expected_df.collect(), f"Expected:\n{expected_df}\nGot:\n{result_df}"


def test_select_nested_fields():
    # Create Sample DataFrame
    df = spark.createDataFrame(
        data=[
            ({"any": {"transcript_ablation.max": 1, "stop_gained.max": 2}},),
            ({"any": {"transcript_ablation.max": 3, "stop_gained.max": 4}},)
        ], schema=StructType([
            StructField("vep", StructType([
                StructField("any", StructType([
                    StructField("transcript_ablation.max", IntegerType(), True),
                    StructField("stop_gained.max", IntegerType(), True)
                ]), True)
            ]), True)
        ])
    )

    # Define Expected Output
    expected_output = [
        ('vep.any.transcript_ablation.max',
         f.col('vep')['any']['transcript_ablation.max'].alias('vep.any.transcript_ablation.max')),
        ('vep.any.stop_gained.max',
         f.col('vep')['any']['stop_gained.max'].alias('vep.any.stop_gained.max'))
    ]

    # Apply Function
    fields = {
        "vep": {
            "any": [
                "transcript_ablation.max",
                "stop_gained.max"
            ]
        }
    }
    result = list(sp_funcs.select_nested_fields(fields))

    for i, x in enumerate(expected_output):
        assert result[i][0] == x[0], f"Expected: {x[0]}, Got: {result[i][0]}"
        assert str(result[i][1]) == str(x[1]), f"Expected: {x[1].name}, Got: {result[i][1].name}"
        assert df.select(result[i][1]).collect() == df.select(x[1]).collect()


def test_melt():
    # Create Sample DataFrame
    df = spark.createDataFrame(
        data=[
            (1, 10, 20),
            (2, 30, 40)
        ],
        schema=StructType([
            StructField("id", IntegerType(), True),
            StructField("A", IntegerType(), True),
            StructField("B", IntegerType(), True)
        ])
    )

    # Define Expected Output
    expected_df = spark.createDataFrame(
        data=[
            (1, "A", 10),
            (1, "B", 20),
            (2, "A", 30),
            (2, "B", 40)
        ],
        schema=StructType([
            StructField("id", IntegerType(), True),
            StructField("variable", StringType(), True),
            StructField("value", IntegerType(), True)
        ])
    )

    # Apply Function
    result_df = sp_funcs.melt(df, id_vars=["id"], value_vars=["A", "B"], var_name="variable", value_name="value")

    # Assert Equality
    assert result_df.collect() == expected_df.collect(), f"Expected:\n{expected_df.collect()}\nGot:\n{result_df.collect()}"


def test_transform_featureset():
    # Create Sample DataFrame
    df = spark.createDataFrame(
        data=[
            (1, {"nested": {"field1": 10, "field2": 20, "field3": 30}}),
            (2, {"nested": {"field1": 30, "field2": 40, "field3": 50}})
        ],
        schema=StructType([
            StructField("id", IntegerType(), True),
            StructField("data", StructType([
                StructField("nested", StructType([
                    StructField("field1", IntegerType(), True),
                    StructField("field2", IntegerType(), True),
                    StructField("field3", IntegerType(), True)
                ]), True)
            ]), True)
        ])
    )

    # Define Variables
    variables = {
        "data": {
            "nested": ["field1", "field3"]
        }
    }

    # Define Index Columns
    index_cols = ["id"]

    # Create Spark DataFrame
    expected_df = spark.createDataFrame(data=[
        (1, 10, 30),
        (2, 30, 50)
    ], schema=StructType([
        StructField("id", IntegerType(), True),
        StructField("feature.example@data.nested.field1", IntegerType(), True),
        StructField("feature.example@data.nested.field3", IntegerType(), True)
    ]))

    # Apply Function
    result_df = sp_funcs.transform_featureset(df, fset_name="example", variables=variables, index_cols=index_cols)

    # Assert Equality
    assert result_df.collect() == expected_df.collect(), f"Expected:\n{expected_df}\nGot:\n{result_df}"


def test_join_featuresets():
    # Create Sample DataFrames
    df1 = spark.createDataFrame(
        data=[
            (1, 3, {"nested": {"field1": 10, "field2": 20}}),
            (2, 3, {"nested": {"field1": 30, "field2": 40}})
        ],
        schema=StructType([
            StructField("id", IntegerType(), True),
            StructField("random", IntegerType(), True),
            StructField("data1", StructType([
                StructField("nested", StructType([
                    StructField("field1", IntegerType(), True),
                    StructField("field2", IntegerType(), True),
                ]), True)
            ]), True)
        ])
    )

    df2 = spark.createDataFrame(
        data=[
            (1, {"nested": {"field1": 5, "field3": 50, "field4": None, 'field5': 70}}),
            (2, {"nested": {"field1": 6, "field3": 70, "field4": 80, 'field5': 90}})
        ],
        schema=StructType([
            StructField("id", IntegerType(), True),
            StructField("data2", StructType([
                StructField("nested", StructType([
                    StructField("field1", IntegerType(), True),
                    StructField("field3", IntegerType(), True),
                    StructField("field4", IntegerType(), True),
                    StructField("field5", IntegerType(), True),
                ]), True)
            ]), True)
        ])
    )
    # `feature`.`df2@data2`.`nested`.`field4`
    # `feature`.`df2@data2`.`nested`.`field4`
    # Define Variables
    variables = {
        "df1": {
            "data1": {
                "nested": ["field1", "field2"]
            }
        },
        "df2": {
            "data2": {
                "nested": ["field3", "field4"]
            }
        }
    }

    # Define Index Columns
    index_cols = ["id"]
    fill_values = {'`feature.df2@data2.nested.field4`': 0}

    # Define Expected Output
    expected_df = spark.createDataFrame(data=[
        (1, 10, 20, 50, 0),
        (2, 30, 40, 70, 80)
    ], schema=StructType([
        StructField("id", IntegerType(), True),
        StructField("feature.example@data.nested.field1", IntegerType(), True),
        StructField("feature.example@data.nested.field1", IntegerType(), True),
        StructField("feature.example@data.nested.field1", IntegerType(), True),
        StructField("feature.example@data.nested.field1", IntegerType(), True),
    ]))

    # Apply Function
    dataframes = {"df1": df1, "df2": df2}
    result_df = sp_funcs.join_featuresets(dataframes, variables, index_cols,
                                          fill_values=fill_values, ignore_missing_columns=False)

    # Assert Equality
    assert result_df.collect() == expected_df.collect(), f"Expected:\n{expected_df}\nGot:\n{result_df}"
