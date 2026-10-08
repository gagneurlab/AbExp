import functools
import pathlib
import warnings

import polars as pl
import tensorflow as tf

from abexp.utils import gff3

# GFF3 attributes that abexp.enformer uses
GFF3_ATTRIBUTES = ('gene_id', 'transcript_id', 'gene_type', 'tag')


def renamed_parameter(old: str, new: str):
    """
    Accept the keyword argument `old` as a deprecated alias of the parameter `new`.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if old in kwargs:
                if new in kwargs:
                    raise TypeError(f'{func.__qualname__}() got both {new!r} and its deprecated alias {old!r}')
                warnings.warn(f'The parameter {old!r} of {func.__qualname__}() is deprecated, use {new!r}.',
                              DeprecationWarning, stacklevel=2)
                kwargs[new] = kwargs.pop(old)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def read_gff3(path: str | pathlib.Path, attributes: tuple[str, ...] = GFF3_ATTRIBUTES) -> pl.DataFrame:
    """
    Read a GFF3 file, such as a GENCODE annotation, into a polars DataFrame with pyranges-style column names.

    This is `abexp.utils.gff3.read_gff3` with `par_y_suffix=True` and other column names.
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
    df = gff3.read_gff3(path, attributes, par_y_suffix=True)
    return genome_annotation_to_polars(df.rename({
        'chrom': 'Chromosome', 'source': 'Source', 'type': 'Feature', 'start': 'Start', 'end': 'End',
        'score': 'Score', 'strand': 'Strand', 'phase': 'Frame',
    }))


@renamed_parameter('gtf', 'genome_annotation')
def genome_annotation_to_polars(genome_annotation) -> pl.DataFrame:
    """
    Get the genome annotation as a polars DataFrame.

    :param genome_annotation: Path to a GFF3 file, or a polars or pandas DataFrame with pyranges-style column names
        (Chromosome, Start, End, Strand, Feature and the GFF3 attributes), such as the output of `read_gff3`
        or of `pyranges.read_gff3(..., as_df=True)`. Start is 0-based, End is 1-based.
        The deprecated alias `gtf` still works.
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
