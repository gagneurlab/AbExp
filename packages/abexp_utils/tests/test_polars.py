import abexp_utils.polars_functions as pl_funcs
import polars as pl


def test_zip_explode_cols():
    # Create Sample DataFrame
    df = pl.DataFrame({
        "A": [1, 2],
        "B": [[4, 5], [6, 7]],
        "C": [[7, 8], [9, 10]]
    })

    # Define Expected Output
    expected_df = pl.DataFrame({
        "A": [1, 1, 2, 2],
        "new_name": [
            {"B": 4, "C": 7},
            {"B": 5, "C": 8},
            {"B": 6, "C": 9},
            {"B": 7, "C": 10}
        ]
    })

    # Apply Function
    result_df = pl_funcs.zip_explode_cols(df, cols=["B", "C"], result_name="new_name")

    # Assert Equality
    assert result_df.frame_equal(expected_df), f"Expected:\n{expected_df}\nGot:\n{result_df}"


def test_select_nested_fields():
    # Create Sample DataFrame
    df = pl.DataFrame({
        "vep": [
            {"any": {"transcript_ablation": {"max": 1}, "stop_gained.max": 2}},
            {"any": {"transcript_ablation": {"max": 3}, "stop_gained.max": 4}}
        ]
    })

    # Define Expected Output
    expected_output = [
        ('vep.any.transcript_ablation.max',
         pl.col('vep').struct.field('any').struct.field('transcript_ablation.max').alias(
             'vep.any.transcript_ablation.max')),
        ('vep.any.stop_gained.max',
         pl.col('vep').struct.field('any').struct.field('stop_gained.max').alias(
             'vep.any.stop_gained.max'))
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
    result = list(pl_funcs.select_nested_fields(fields))

    for i, x in enumerate(expected_output):
        assert result[i][0] == x[0], f"Expected: {x[0]}, Got: {result[i][0]}"
        assert str(result[i][1]) == str(x[1]), f"Expected: {x[1].name}, Got: {result[i][1].name}"


def test_transform_featureset():
    # Create Sample DataFrame
    df = pl.DataFrame({
        "id": [1, 2],
        "data": [
            {"nested": {"field1": 10, "field2": 20, 'field3': 30}},
            {"nested": {"field1": 30, "field2": 40, 'field3': 50}}
        ]
    })

    # Define Variables
    variables = {
        "data": {
            "nested": ["field1", "field3"]
        }
    }

    # Define Index Columns
    index_cols = ["id"]

    # Step 5: Define Expected Output
    expected_df = pl.DataFrame({
        "id": [1, 2],
        "feature.example@data.nested.field1": [10, 30],
        "feature.example@data.nested.field3": [30, 50]
    })

    # Apply Function
    result_df = pl_funcs.transform_featureset(df, fset_name="example", variables=variables, index_cols=index_cols)

    # Assert Equality
    assert result_df.frame_equal(expected_df), f"Expected:\n{expected_df}\nGot:\n{result_df}"


def test_join_featuresets():
    # Create Sample DataFrames
    df1 = pl.DataFrame({
        "id": [1, 2],
        "random": [3, 3],
        "data1": [
            {"nested": {"field1": 10, "field2": 20}},
            {"nested": {"field1": 30, "field2": 40}}
        ]
    })

    df2 = pl.DataFrame({
        "id": [1, 2],
        "data2": [
            {"nested": {"field1": 5, "field3": 50, "field4": None, 'field5': 70}},
            {"nested": {"field1": 6, "field3": 70, "field4": 80, 'field5': 90}}
        ]
    })

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
    fill_values = {'feature.df2@data2.nested.field4': 0}
    broadcast_columns = {
        "broadcast": pl.DataFrame({"test": [1, 2]})
    }

    # Define Expected Output
    expected_df = pl.DataFrame({
        "id": [1, 2],
        "feature.df1@data1.nested.field1": [10, 30],
        "feature.df1@data1.nested.field2": [20, 40],
        "feature.df2@data2.nested.field3": [50, 70],
        "feature.df2@data2.nested.field4": [0, 80]
    })

    # Apply Function
    dataframes = {"df1": df1, "df2": df2}
    result_df = pl_funcs.join_featuresets(dataframes, variables, index_cols,
                                          fill_values=fill_values, broadcast_columns=broadcast_columns,
                                          ignore_missing_columns=False)

    # Assert Equality
    assert result_df.frame_equal(expected_df), f"Expected:\n{expected_df}\nGot:\n{result_df}"
