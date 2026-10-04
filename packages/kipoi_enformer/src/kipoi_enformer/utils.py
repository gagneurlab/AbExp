import pathlib
import urllib.parse

import polars as pl
import polars_bio as pb
import tensorflow as tf

# GFF3 attributes that kipoi_enformer uses
GFF3_ATTRIBUTES = ('gene_id', 'transcript_id', 'gene_type', 'tag')


def read_gff3(path: str | pathlib.Path, attributes: tuple[str, ...] = GFF3_ATTRIBUTES) -> pl.DataFrame:
    """
    Read a GFF3 file, such as a GENCODE annotation, into a polars DataFrame with pyranges-style column names.

    Start is 0-based, End is 1-based (0-based, half-open).
    An attribute with several values, such as `tag`, keeps them joined with ",", as in the file.
    Escaped characters, such as "%3B" for ";", are decoded.
    GENCODE marks the chrY PAR copies only in `ID`. Their `gene_id` and `transcript_id` get the suffix "_PAR_Y",
    as in the GENCODE GTF, so that the IDs stay unique.

    :param path: Path to the GFF3 file
    :param attributes: GFF3 attributes to read as columns
    :return: DataFrame with the columns Chromosome, Source, Feature, Start, End, Score, Strand, Frame
        and one column per attribute
    """
    if pathlib.Path(path).name.removesuffix('.gz').endswith('.gtf'):
        raise ValueError(f'{path} looks like a GTF file, but kipoi_enformer reads the genome annotation from GFF3')
    df = pb.read_gff(str(path), attr_fields=['ID', *attributes], use_zero_based=True)
    # polars-bio decodes only some escapes, e.g. "%3B" but not "%25", so decode the rest
    for name in ['ID', *attributes]:
        escaped = df[name].drop_nulls().unique()
        escaped = escaped.filter(escaped.str.contains('%', literal=True))
        if len(escaped) > 0:
            df = df.with_columns(pl.col(name).replace(escaped, [urllib.parse.unquote(x) for x in escaped]))
    # ID marks a chrY PAR copy with the suffix "_PAR_Y", e.g. ENST00000381192.10_PAR_Y or
    # exon:ENST00000381192.10_PAR_Y:1, or in some lift37 entries with the prefix "ENSTR" or "ENSGR"
    is_par_y = pl.col('ID').str.contains(r'_PAR_Y|(^|:)ENS[GT]R\d')
    df = df.with_columns(
        pl.when(is_par_y & ~pl.col(c).str.ends_with('_PAR_Y')).then(pl.col(c) + '_PAR_Y').otherwise(pl.col(c)).alias(c)
        for c in ('gene_id', 'transcript_id') if c in df.columns
    )
    if 'ID' not in attributes:
        df = df.drop('ID')
    return genome_annotation_to_polars(df.rename({
        'chrom': 'Chromosome', 'source': 'Source', 'type': 'Feature', 'start': 'Start', 'end': 'End',
        'score': 'Score', 'strand': 'Strand', 'phase': 'Frame',
    }))


def genome_annotation_to_polars(genome_annotation) -> pl.DataFrame:
    """
    Get the genome annotation as a polars DataFrame.

    :param genome_annotation: Path to a GFF3 file, or a polars or pandas DataFrame with pyranges-style column names
        (Chromosome, Start, End, Strand, Feature and the GFF3 attributes), such as the output of `read_gff3`
        or of `pyranges.read_gff3(..., as_df=True)`. Start is 0-based, End is 1-based.
    :return: polars DataFrame with Chromosome, Strand and Feature as strings and Start and End as Int64
    """
    if isinstance(genome_annotation, (str, pathlib.Path)):
        return read_gff3(genome_annotation)
    df = genome_annotation
    if not isinstance(df, pl.DataFrame):
        # e.g. a pandas DataFrame
        df = pl.from_pandas(df)
    return df.with_columns(
        *[pl.col(c).cast(pl.String) for c in ['Chromosome', 'Strand', 'Feature'] if c in df.columns],
        *[pl.col(c).cast(pl.Int64) for c in ['Start', 'End'] if c in df.columns],
    )


class RandomModel(tf.keras.Model):
    """
    A random model for testing purposes.
    """

    def __init__(self, lamda=10):
        super().__init__()
        self.lamda = lamda

    def predict_on_batch(self, input_tensor):
        # tf.random.set_seed(42)
        # Uniform integers from 0 to 2 * lamda have about the mean lamda, like poisson values, and are much faster to
        # draw. They also have few distinct values, which compress well in parquet.
        # Enformer reads only the human tracks.
        tracks = tf.random.uniform((input_tensor.shape[0], 896, 5313,), maxval=int(2 * self.lamda) + 1, dtype=tf.int32)
        return {
            'human': tf.cast(tracks, tf.float32),
        }
