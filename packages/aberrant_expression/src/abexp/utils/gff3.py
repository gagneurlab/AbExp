"""Read GFF3 files, such as the GENCODE annotation, with polars-bio. Needs the extra `gff3`."""
import pathlib
import urllib.parse

import polars as pl
import polars_bio as pb

# GENCODE marks the chrY PAR copies only in `ID`: with the suffix "_PAR_Y", e.g. ENST00000381192.10_PAR_Y or
# exon:ENST00000381192.10_PAR_Y:1, or in some lift37 entries with the prefix "ENSTR" or "ENSGR", e.g.
# ENSTR0000302805.2 or exon:ENSTR0000302805.2:1
_PAR_Y_ID = r'_PAR_Y|(^|:)ENS[GT]R'


def read_gff3(path: str | pathlib.Path, attributes: tuple[str, ...], *, par_y_suffix: bool = False) -> pl.DataFrame:
    """
    Read a GFF3 file, such as a GENCODE annotation, into a polars DataFrame.

    The rows keep the order of the file. Start is 0-based, End is 1-based (0-based, half-open).
    An attribute with several values, such as `tag`, keeps them joined with ",", as in the file.
    Escaped characters, such as "%3B" for ";", are decoded.

    :param path: Path to the GFF3 file, may be gzipped
    :param attributes: GFF3 attributes to read as columns
    :param par_y_suffix: Add the suffix "_PAR_Y" to `gene_id` and `transcript_id` of the chrY PAR copies, as in the
        GENCODE GTF, so that the IDs stay unique. GENCODE marks these copies only in `ID`.
    :return: DataFrame with the columns of polars-bio: chrom, start and end (Int64), type, source, score, strand,
        phase, and one column per attribute
    """
    if pathlib.Path(path).name.removesuffix('.gz').endswith('.gtf'):
        raise ValueError(f'{path} looks like a GTF file, but abexp reads the genome annotation from GFF3')
    read_id = par_y_suffix and 'ID' not in attributes
    columns = ('ID', *attributes) if read_id else attributes
    df = pb.read_gff(str(path), attr_fields=list(columns), use_zero_based=True)
    # polars-bio decodes only some escapes, e.g. "%3B" but not "%25", so decode the rest
    for name in columns:
        escaped = df[name].drop_nulls().unique()
        escaped = escaped.filter(escaped.str.contains('%', literal=True))
        if len(escaped) > 0:
            df = df.with_columns(pl.col(name).replace(escaped, [urllib.parse.unquote(x) for x in escaped]))
    if par_y_suffix:
        is_par_y = pl.col('ID').str.contains(_PAR_Y_ID)
        df = df.with_columns(
            pl.when(is_par_y & ~pl.col(c).str.ends_with('_PAR_Y'))
            .then(pl.col(c) + '_PAR_Y')
            .otherwise(pl.col(c))
            .alias(c)
            for c in ('gene_id', 'transcript_id') if c in df.columns
        )
    if read_id:
        df = df.drop('ID')
    return df.with_columns(pl.col('start', 'end').cast(pl.Int64))
