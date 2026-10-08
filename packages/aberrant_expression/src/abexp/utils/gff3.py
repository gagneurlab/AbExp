"""Read GFF3 files, such as the GENCODE annotation, with polars-bio. Needs the extra `enformer` or `pangolin`."""
import pathlib
import urllib.parse

import polars as pl
import polars_bio as pb


def read_gff3(path: str | pathlib.Path, attributes: tuple[str, ...]) -> pl.DataFrame:
    """
    Read a GFF3 file, such as a GENCODE annotation, into a polars DataFrame.

    The rows keep the order of the file. Start is 0-based, End is 1-based (0-based, half-open).
    An attribute with several values, such as `tag`, keeps them joined with ",", as in the file.
    Escaped characters, such as "%3B" for ";", are decoded.

    :param path: Path to the GFF3 file, may be gzipped
    :param attributes: GFF3 attributes to read as columns
    :return: DataFrame with the columns of polars-bio: chrom, start and end (Int64), type, source, score, strand,
        phase, and one column per attribute
    """
    if pathlib.Path(path).name.removesuffix('.gz').endswith('.gtf'):
        raise ValueError(f'{path} looks like a GTF file, but abexp reads the genome annotation from GFF3')
    df = pb.read_gff(str(path), attr_fields=list(attributes), use_zero_based=True)
    # polars-bio decodes only some escapes, e.g. "%3B" but not "%25", so decode the rest
    for name in attributes:
        escaped = df[name].drop_nulls().unique()
        escaped = escaped.filter(escaped.str.contains('%', literal=True))
        if len(escaped) > 0:
            df = df.with_columns(pl.col(name).replace(escaped, [urllib.parse.unquote(x) for x in escaped]))
    return df.with_columns(pl.col('start', 'end').cast(pl.Int64))
