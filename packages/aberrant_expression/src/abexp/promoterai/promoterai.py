import math
import pathlib

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
import tensorflow as tf
from tqdm.autonotebook import tqdm

from abexp.promoterai.constants import MAX_POSITIONS_PER_CALL, SCORE_DECIMALS
from abexp.promoterai.dataloader import VCFTSSDataloader
from abexp.promoterai.logger import logger
from abexp.promoterai.utils import RandomModel, genome_annotation_to_polars

__all__ = ['PromoterAI', 'PromoterAIVeff']

# the aggregation modes of PromoterAIVeff.run
AGGREGATION_MODES = ('canonical', 'median')


def load_model(model_path: str | pathlib.Path):
    """
    Load a PromoterAI model from the folder that Illumina publishes under its license agreement.

    The published models are Keras 2 models in the TensorFlow SavedModel format. `tf.saved_model.load` restores
    them without Keras, so that they also load with Keras 3, which no longer reads this format. A model in the
    Keras format (`.keras`) is loaded with `tf.keras.models.load_model` instead.

    :param model_path: Path to the model folder, e.g. promoterAI_v1_hg38_finetune
    :return: the model, callable on a batch of sequences
    """
    model_path = pathlib.Path(model_path)
    if (model_path / 'saved_model.pb').exists():
        return tf.saved_model.load(str(model_path))
    return tf.keras.models.load_model(str(model_path))


def human_output(outputs):
    """
    The output of the human head of a PromoterAI model.

    The published models return a dict with the key "human", the legacy models a tuple or list with the human head
    first, see `promoterai.architecture._get_human_output` of PromoterAI.

    :param outputs: the outputs of the model
    :return: the output of the human head as a tensor
    """
    if isinstance(outputs, (tuple, list)):
        return outputs[0]
    if isinstance(outputs, dict):
        if 'human' in outputs:
            return outputs['human']
        # the signature of a SavedModel names its outputs, e.g. "output0" for the human head of PromoterAI
        return outputs[sorted(outputs.keys())[0]]
    return outputs


class PromoterAI:
    def __init__(self, model_path: str | pathlib.Path | None = None, is_random: bool = False, **random_kwargs):
        """
        :param model_path: Path to the model folder, see `load_model`. Required unless `is_random` is True.
        :param is_random: If True, load a model with fixed random weights for testing purposes.
        """
        if is_random:
            self._model = RandomModel(**random_kwargs)
        elif model_path is None:
            raise ValueError('model_path is required: PromoterAI publishes its models under a license agreement, '
                             'see https://github.com/Illumina/PromoterAI')
        else:
            logger.debug(f'Loading model from {model_path}')
            self._model = load_model(model_path)

    def predict(self, dataloader: VCFTSSDataloader, batch_size: int, filepath: str | pathlib.Path):
        """
        Predict on a dataloader and save the results in a parquet file.

        The file holds the metadata of the dataloader and the column `score`: the PromoterAI score of the variant
        for the strand of the transcript, between -1 and 1, where a negative score means under- and a positive
        score over-expression. The score is NaN for a variant that the dataloader could not score, see
        `abexp.promoterai.dataloader.dataloader.extract_variant_sequences`.

        :param dataloader: the dataloader with the variants
        :param batch_size: the number of records per batch
        :param filepath: the parquet file to write
        :return:
        """
        logger.debug('Predicting on dataloader')
        assert batch_size > 0

        # Hint: order matters
        schema = dataloader.pyarrow_metadata_schema
        schema = schema.insert(0, pa.field('score', pa.float32()))

        batch_counter = 0
        total_batches = math.ceil(len(dataloader) / batch_size)
        if total_batches == 0:
            logger.info('The dataloader is empty. No predictions to make.')
            writer = pq.ParquetWriter(filepath, schema)
            writer.close()
            return

        with pq.ParquetWriter(filepath, schema) as writer:
            for batch in tqdm(dataloader.batch_iter(batch_size=batch_size), total=total_batches):
                batch_counter += 1
                batch = self._to_pyarrow(self._process_batch(batch))
                writer.write_batch(batch)

        # sanity check for the dataloader
        assert batch_counter == total_batches

    def _process_batch(self, batch):
        """
        Process a batch of data. Run the model on the reference and the alternative sequences and prepare the
        results dict.

        The score is the tanh of the mean difference of the outputs, as `promoterai.score` of PromoterAI
        calculates it: the twin model of `promoterai.architecture.twin_wrap` subtracts the outputs of the human
        head and averages over the positions and the tracks, and the scoring script rounds the result to
        `SCORE_DECIMALS` decimals before the tanh.

        :param batch: batch of the dataloader
        :return: Results dict. Structure: {'metadata': {field: [values]}, 'score': [values]}
        """
        ref_predictions = self._predict_on_batch(batch['sequences']['ref'])
        alt_predictions = self._predict_on_batch(batch['sequences']['alt'])
        assert ref_predictions.shape == alt_predictions.shape, \
            f'the outputs of the sequences differ in shape: {ref_predictions.shape}, {alt_predictions.shape}'
        assert ref_predictions.shape[0] == len(batch['scorable'])

        # average the difference over the positions and the tracks of the human head
        axis = tuple(range(1, ref_predictions.ndim))
        score = np.tanh(np.round((alt_predictions - ref_predictions).mean(axis=axis), SCORE_DECIMALS))
        score = np.where(batch['scorable'], score, np.nan)

        return {
            'metadata': batch['metadata'],
            'score': score.astype(np.float32),
        }

    def _predict_on_batch(self, sequences: np.ndarray) -> np.ndarray:
        """
        Run the model on a batch of one-hot encoded sequences.

        A batch with `MAX_POSITIONS_PER_CALL` positions or more is split into several calls of the model, see
        `abexp.promoterai.constants.MAX_POSITIONS_PER_CALL`.

        :param sequences: float32 array of shape (batch, length, 4)
        :return: the output of the human head as a float32 array of shape (batch, positions, tracks)
        """
        call_size = max(1, (MAX_POSITIONS_PER_CALL - 1) // sequences.shape[1])
        if len(sequences) > call_size:
            return np.concatenate([self._call_model(sequences[i:i + call_size])
                                   for i in range(0, len(sequences), call_size)])
        return self._call_model(sequences)

    def _call_model(self, sequences: np.ndarray) -> np.ndarray:
        """
        Run the model once on a batch of one-hot encoded sequences, see `_predict_on_batch`.
        """
        input_tensor = tf.convert_to_tensor(sequences, dtype=tf.float32)
        model = self._model
        if hasattr(model, 'predict_on_batch'):
            # a Keras model, or RandomModel
            outputs = model.predict_on_batch(input_tensor)
        elif callable(model):
            # a SavedModel restored by tf.saved_model.load
            outputs = model(input_tensor)
        else:
            outputs = model.signatures['serving_default'](input_tensor)
        outputs = human_output(outputs)
        return np.asarray(outputs, dtype=np.float32)

    @staticmethod
    def _to_pyarrow(results: dict):
        """
        Convert the results dict of the _process_batch method to a pyarrow record batch.
        :param results: the results dict
        :return: pyarrow.RecordBatch object
        """
        logger.debug('Converting results to pyarrow')

        metadata = {k: pa.array(v.tolist()) for k, v in results['metadata'].items()}
        formatted_results = {
            'score': pa.array(results['score'], type=pa.float32()),
            **metadata
        }

        logger.debug('Constructing pyarrow record batch')
        return pa.RecordBatch.from_arrays(list(formatted_results.values()), names=list(formatted_results.keys()))


class PromoterAIVeff:
    # the columns of the output of `run`
    OUTPUT_SCHEMA = {
        'chrom': pl.String,
        'strand': pl.String,
        'gene_id': pl.String,
        'variant_start': pl.Int64,
        'variant_end': pl.Int64,
        'ref': pl.String,
        'alt': pl.String,
        'veff_score': pl.Float32,
    }

    def __init__(self, genome_annotation=None):
        """
        :param genome_annotation: The path to a GFF3 file or a polars or pandas DataFrame containing the genome
            annotation, see `abexp.promoterai.utils.genome_annotation_to_polars`. Required for the aggregation
            mode "canonical".
        """
        # if a genome annotation is given, then extract the canonical transcripts for the canonical aggregation mode
        self.canonical_transcripts = None
        if genome_annotation is not None:
            annotation = genome_annotation_to_polars(genome_annotation)
            # only keep protein_coding transcripts
            annotation = annotation.filter(pl.col('gene_type') == 'protein_coding')
            # check if Ensembl_canonical is in the set of tags
            annotation = annotation.filter(
                pl.col('tag').str.split(',').list.contains('Ensembl_canonical').fill_null(False))
            # the transcript IDs without their version, as `run` writes them
            self.canonical_transcripts = annotation['transcript_id'].str.replace(r'\..+$', '').unique()

    def run(self, score_path: str | pathlib.Path, output_path: str | pathlib.Path, aggregation_mode: str,
            upstream_tss: int | None = None, downstream_tss: int | None = None):
        """
        Given a file with the PromoterAI scores of the variant-transcript pairs, as `PromoterAI.predict` writes it,
        aggregate the scores by gene and variant. Save the results in a parquet file.

        A PromoterAI score depends only on the variant and on the strand of the transcript, so all transcripts of a
        gene have the same score and the aggregation modes differ only in which transcripts they keep.

        :param score_path: The parquet file that contains the scores
        :param output_path: The parquet file that will contain the variant effect scores
        :param aggregation_mode: One of ['canonical', 'median'].
        :param upstream_tss: Variants outside of the interval [-upstream_tss, downstream_tss] will be dropped.
        :param downstream_tss: Variants outside of the interval [-upstream_tss, downstream_tss] will be dropped.
        :return:
        """

        if aggregation_mode not in AGGREGATION_MODES:
            raise ValueError(f'Unknown mode: {aggregation_mode}')
        elif aggregation_mode == 'canonical':
            assert self.canonical_transcripts is not None, 'Canonical transcripts are required for this mode.'

        logger.debug(f'Calculating the variant effect for {score_path}')

        veff_ldf = pl.scan_parquet(score_path)

        # check if the scores are empty and write an empty file if that's the case
        if veff_ldf.select(pl.len()).collect().item() == 0:
            logger.warning('The scores are empty. No variant effect to calculate.')
            pl.DataFrame(schema=self.OUTPUT_SCHEMA).write_parquet(output_path)
            return

        veff_ldf = veff_ldf.select(['tss', 'chrom', 'strand',
                                    'gene_id', 'transcript_id', 'transcript_start', 'transcript_end',
                                    'variant_start', 'variant_end', 'ref', 'alt', 'score'])

        # drop the variants that PromoterAI could not score, see PromoterAI.predict
        veff_ldf = veff_ldf.filter(pl.col('score').is_not_nan() & pl.col('score').is_not_null())
        # filter out Y chromosome equivalent transcripts
        veff_ldf = veff_ldf.filter(~pl.col('transcript_id').str.contains('_PAR_Y'))
        # remove gene and transcript versions
        veff_ldf = veff_ldf.with_columns(
            pl.col('gene_id').str.replace(r'([^\.]+)\..+$', "${1}").alias('gene_id'),
            pl.col('transcript_id').str.replace(r'([^\.]+)\..+$', "${1}").alias('transcript_id')
        )

        # calculate variant position relative to the tss
        if downstream_tss is not None or upstream_tss is not None:
            veff_ldf = veff_ldf.with_columns(
                (pl.when(pl.col('strand') == '+').then(
                    pl.col('variant_start') - pl.col('transcript_start')
                ).otherwise(
                    pl.col('transcript_end') - pl.col('variant_start'))
                ).alias('relative_pos'))

            if downstream_tss is not None:
                veff_ldf = veff_ldf.filter(pl.col('relative_pos') <= downstream_tss)
            if upstream_tss is not None:
                veff_ldf = veff_ldf.filter(pl.col('relative_pos') >= -upstream_tss)

        veff_df = self._aggregate(veff_ldf, aggregation_mode)
        logger.debug(f'Writing the variant effect to {output_path}')
        veff_df = veff_df.rename({'score': 'veff_score'})
        veff_df.write_parquet(output_path)

    def _aggregate(self, veff_ldf, aggregation_mode: str):
        """
        Given a dataframe containing the scores of the variant-transcript pairs, aggregate the scores by gene and
        variant.

        :param veff_ldf: A polars lazy frame containing the scores.
        :param aggregation_mode: One of ['canonical', 'median'].
        :return: A polars DataFrame containing the aggregated scores.
        """

        group_keys = ['chrom', 'strand', 'gene_id', 'variant_start', 'variant_end', 'ref', 'alt']

        if aggregation_mode == 'canonical':
            # Keep only the canonical transcripts
            veff_ldf = veff_ldf.filter(pl.col('transcript_id').is_in(self.canonical_transcripts.to_list()))
            veff_ldf = veff_ldf.group_by(group_keys).agg(
                pl.col('score').first(),
                pl.len().alias('num_transcripts')
            )
            veff_df = veff_ldf.collect()
            # Verify that there is only one canonical transcript per gene
            max_transcripts_per_gene = veff_df['num_transcripts'].max()
            if max_transcripts_per_gene is not None and max_transcripts_per_gene > 1:
                logger.error('Multiple canonical transcripts found for a gene.')
                logger.error(veff_df.filter(pl.col('num_transcripts') > 1))
                raise ValueError('Multiple canonical transcripts found for a gene.')

            # Remove the num_transcripts column
            veff_df.drop_in_place('num_transcripts')
        elif aggregation_mode == 'median':
            veff_ldf = veff_ldf.group_by(group_keys).agg(pl.col('score').median().cast(pl.Float32))
            veff_df = veff_ldf.collect()
        else:
            raise ValueError(f'Unknown mode: {aggregation_mode}')

        # group_by returns the groups in a random order
        veff_df = veff_df.sort(group_keys)
        logger.debug(f'Aggregated table size: {len(veff_df)}')
        return veff_df
