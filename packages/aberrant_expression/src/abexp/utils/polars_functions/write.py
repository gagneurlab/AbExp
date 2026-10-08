def df_batch_writer(df_iter, output):
    """Write the polars DataFrames of `df_iter` to one CSV file.

    Raises StopIteration if `df_iter` is empty.
    """
    df = next(df_iter)
    with open(output, 'wb') as f:
        df.write_csv(f)
        for df in df_iter:
            df.write_csv(f, include_header=False)
