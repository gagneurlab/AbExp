import polars as pl
import polars_bio as pb
import tensorflow as tf
import pathlib

# GTF attributes that kipoi_enformer uses
GTF_ATTRIBUTES = ('gene_id', 'transcript_id', 'gene_type', 'tag')


def read_gtf(gtf: str | pathlib.Path, attributes: tuple[str, ...] = GTF_ATTRIBUTES) -> pl.DataFrame:
    """
    Read a GTF file into a polars DataFrame with pyranges-style column names.

    Start is 0-based, End is 1-based (0-based, half-open).
    Repeated attributes, such as `tag`, are joined with ",".

    :param gtf: Path to GTF file
    :param attributes: GTF attributes to read as columns
    :return: DataFrame with the columns Chromosome, Source, Feature, Start, End, Score, Strand, Frame
        and one column per attribute
    """
    df = pb.read_gtf(str(gtf), attr_fields=list(attributes), use_zero_based=True)
    return genome_annotation_to_polars(df.rename({
        'chrom': 'Chromosome', 'source': 'Source', 'type': 'Feature', 'start': 'Start', 'end': 'End',
        'score': 'Score', 'strand': 'Strand', 'phase': 'Frame',
    }))


def genome_annotation_to_polars(gtf) -> pl.DataFrame:
    """
    Get the genome annotation as a polars DataFrame.

    :param gtf: Path to a GTF file, or a polars or pandas DataFrame with pyranges-style column names
        (Chromosome, Start, End, Strand, Feature and the GTF attributes), such as the output of `read_gtf`
        or of `pyranges.read_gtf(..., as_df=True)`. Start is 0-based, End is 1-based.
    :return: polars DataFrame with Chromosome, Strand and Feature as strings and Start and End as Int64
    """
    if isinstance(gtf, (str, pathlib.Path)):
        return read_gtf(gtf)
    if not isinstance(gtf, pl.DataFrame):
        # e.g. a pandas DataFrame
        gtf = pl.from_pandas(gtf)
    return gtf.with_columns(
        *[pl.col(c).cast(pl.String) for c in ['Chromosome', 'Strand', 'Feature'] if c in gtf.columns],
        *[pl.col(c).cast(pl.Int64) for c in ['Start', 'End'] if c in gtf.columns],
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
        return {
            'human': tf.abs(tf.random.poisson((input_tensor.shape[0], 896, 5313,), lam=self.lamda)),
            'mouse': tf.abs(tf.random.poisson((input_tensor.shape[0], 896, 1643), lam=self.lamda)),
        }
