import pathlib

import polars as pl
import tensorflow as tf

from abexp.utils import gff3

# GFF3 attributes that abexp.promoterai uses
GFF3_ATTRIBUTES = ('gene_id', 'transcript_id', 'gene_type', 'tag')


def read_gff3(path: str | pathlib.Path, attributes: tuple[str, ...] = GFF3_ATTRIBUTES) -> pl.DataFrame:
    """
    Read a GFF3 file, such as a GENCODE annotation, into a polars DataFrame with pyranges-style column names.

    This is `abexp.utils.gff3.read_gff3` with `par_y_suffix=True` and other column names, the same as
    `abexp.enformer.utils.read_gff3`.
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


def genome_annotation_to_polars(genome_annotation) -> pl.DataFrame:
    """
    Get the genome annotation as a polars DataFrame.

    :param genome_annotation: Path to a GFF3 file, or a polars or pandas DataFrame with pyranges-style column names
        (Chromosome, Start, End, Strand, Feature and the GFF3 attributes), such as the output of `read_gff3`
        or of the gtf_transcripts module of the workflow. Start is 0-based, End is 1-based.
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
    A model with fixed random weights, for testing purposes.

    It maps each base of the sequence to `num_tracks` outputs with a weight matrix that depends only on `seed`, so
    its outputs are a deterministic function of the input: the same sequence always gives the same outputs, and two
    equal sequences give a variant effect of 0. `output_type` selects the output format of the published models:
    a dictionary with the key "human", a tuple of the human and the non-human head, or the bare human output.
    """

    def __init__(self, num_tracks: int = 5, seed: int = 0, output_type: str = 'dict'):
        super().__init__()
        if output_type not in ('dict', 'tuple', 'tensor'):
            raise ValueError(f'Unknown output_type: {output_type}')
        self.num_tracks = num_tracks
        self.seed = seed
        self.output_type = output_type
        # `weights` is a property of tf.keras.Model, so the weight matrix needs another name
        self.track_weights = tf.random.stateless_uniform((4, num_tracks), seed=[seed, 0],
                                                         minval=-1, maxval=1)

    def predict_on_batch(self, input_tensor):
        # (batch, length, 4) x (4, num_tracks) -> (batch, length, num_tracks)
        human = tf.matmul(tf.cast(input_tensor, tf.float32), self.track_weights)
        if self.output_type == 'dict':
            return {'human': human, 'nonhuman': tf.zeros_like(human)}
        if self.output_type == 'tuple':
            return human, tf.zeros_like(human)
        return human

