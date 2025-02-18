import pathlib
import numpy as np
import tensorflow_hub as hub
import tensorflow as tf
from kipoi_enformer.dataloader import TSSDataloader
from kipoi_enformer.utils import RandomModel, gtf_to_pandas
from kipoi_enformer.logger import logger
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm.autonotebook import tqdm
import math
import yaml
import pickle
import polars as pl
from scipy.special import logsumexp
import xarray as xr
from sklearn import linear_model, pipeline, preprocessing
import pandas as pd
import sklearn as sk

__all__ = ['Enformer', 'EnformerAggregator', 'EnformerTissueMapper', 'EnformerVeff']

# Enformer model URI
MODEL_PATH = 'https://tfhub.dev/deepmind/enformer/1'


class Enformer:
    NUM_HUMAN_TRACKS = 5313
    # length of the bins in the enformer model
    BIN_SIZE = 128
    NUM_PREDICTION_BINS = 896
    NUM_SEEN_BINS = 1536
    # length of central sequence for which enformer gives predictions (896 bins)
    # ─────┆─────┆════════════════════════┆─────┆─────
    PRED_SEQUENCE_LENGTH = NUM_PREDICTION_BINS * BIN_SIZE
    # length of central sequence which enformer actually sees (1536 bins)
    # ─────┆═════┆════════════════════════┆═════┆─────
    SEEN_SEQUENCE_LENGTH = NUM_SEEN_BINS * BIN_SIZE

    def __init__(self, is_random: bool = False, **random_kwargs):
        """
        :param is_random: If True, load a random model for testing purposes.
        """
        if not is_random:
            logger.debug(f'Loading model from {MODEL_PATH}')
            self._model = hub.load(MODEL_PATH).model
        else:
            self._model = RandomModel(**random_kwargs)

    def predict(self, dataloader: TSSDataloader, batch_size: int, filepath: str | pathlib.Path,
                num_output_bins=NUM_PREDICTION_BINS):
        """
        Predict on a dataloader and save the results in a parquet file
        :param num_output_bins: The number of bins to extract from enformer's output
        :param filepath:
        :param dataloader:
        :param batch_size:
        :return: filepath to the parquet dataset
        """
        logger.debug('Predicting on dataloader')
        assert batch_size > 0

        # Hint: order matters
        schema = dataloader.pyarrow_metadata_schema
        schema = schema.insert(0, pa.field(f'tracks', pa.list_(pa.list_(pa.list_(pa.float32())))))

        shifts = [int(x) for x in schema.metadata[b'shifts'].split(b';')]
        max_abs_shift = max([abs(shift) for shift in shifts])
        assert math.ceil(max_abs_shift / self.BIN_SIZE) < num_output_bins <= self.NUM_PREDICTION_BINS, \
            f'num_output_bins must be fit the maximum shift and be at most {self.NUM_PREDICTION_BINS}'

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
                batch = self._to_pyarrow(self._process_batch(batch, num_output_bins=num_output_bins))
                writer.write_batch(batch)

        # sanity check for the dataloader
        assert batch_counter == total_batches

    def _process_batch(self, batch, num_output_bins=11):
        """
        Process a batch of data. Run the model and prepare the results dict.
        :param batch: list of data dicts
        :return: Results dict. Structure: {'metadata': {field: [values]}, 'predictions': {sequence_key: [values]}}
        """
        batch_size = batch['sequences'].shape[0]
        seqs_per_record = batch['sequences'].shape[1]
        sequences = np.reshape(batch['sequences'], (batch_size * seqs_per_record, 393_216, 4))

        # create input tensor
        input_tensor = tf.convert_to_tensor(sequences)
        assert input_tensor.shape == (batch_size * seqs_per_record, 393_216, 4)

        # run model
        predictions = self._model.predict_on_batch(input_tensor)['human'].numpy()
        assert predictions.shape == (batch_size * seqs_per_record, self.NUM_PREDICTION_BINS, self.NUM_HUMAN_TRACKS)
        predictions = predictions.reshape(batch_size, seqs_per_record, self.NUM_PREDICTION_BINS, self.NUM_HUMAN_TRACKS)

        # extract central bins if the number of output bins is different from the number of prediction bins
        if num_output_bins != self.NUM_PREDICTION_BINS:
            # calculate TSS bin
            tss_bin = (self.PRED_SEQUENCE_LENGTH // 2 + 1) // self.BIN_SIZE
            bins = [tss_bin + i for i in range(-math.floor(num_output_bins / 2), math.ceil(num_output_bins / 2))]
            assert len(bins) == num_output_bins
            predictions = predictions[:, :, bins, :]

        results = {
            'metadata': batch['metadata'],
            'tracks': predictions
        }
        return results

    @staticmethod
    def _to_pyarrow(results: dict):
        """
        Convert the results dict from the _process_batch method to a pyarrow table and write it in a parquet file.
        :param results: pyarrow.RecordBatch object
        """
        logger.debug('Converting results to pyarrow')

        # format predictions
        metadata = {}
        for k, v in results['metadata'].items():
            v = pa.array(v.tolist())
            if isinstance(v, np.ndarray):
                v = v.tolist()
            metadata[k] = v

        formatted_results = {
            'tracks': pa.array(results['tracks'].tolist(), type=pa.list_(pa.list_(pa.list_(pa.float32())))),
            **metadata
        }

        logger.debug('Constructing pyarrow record batch')
        # construct RecordBatch
        return pa.RecordBatch.from_arrays(list(formatted_results.values()), names=list(formatted_results.keys()))


class EnformerAggregator:
    def aggregate(self, enformer_scores_path: str | pathlib.Path, output_path: str | pathlib.Path, num_bins: int = 3):
        """
        Aggregate enformer predictions over the bins centered at the tss bin, and the shifts.
        :param enformer_scores_path:
        :param output_path:
        :param num_bins:
        :return:
        """

        enformer_file = pq.ParquetFile(enformer_scores_path)
        enformer_schema = enformer_file.schema.to_arrow_schema()
        metadata = enformer_schema.metadata
        metadata['nbins'] = str(num_bins)
        output_schema = enformer_schema.with_metadata(metadata)
        output_schema = output_schema.remove(0). \
            insert(0, pa.field('tracks', pa.list_(pa.float32(), list_size=Enformer.NUM_HUMAN_TRACKS)))
        # fix polars string issue when transforming to pyarrow
        for idx, x in enumerate(output_schema):
            if x.type == pa.string():
                output_schema = output_schema.remove(idx).insert(idx, pa.field(x.name, pa.large_string()))

        shifts = [int(x) for x in metadata[b'shifts'].split(b';')]

        logger.info(f'Iterating over the parquet files in {enformer_scores_path}')
        with pq.ParquetWriter(output_path, output_schema) as writer:
            for i in tqdm(range(enformer_file.num_row_groups)):
                df = self._aggregate_batch(pl.from_arrow(enformer_file.read_row_group(i)),
                                           Enformer.BIN_SIZE, num_bins, shifts)
                logger.debug('Writing to file')
                writer.write(df.to_arrow())

    @staticmethod
    def _aggregate_batch(frame, bin_size, num_bins, shifts):
        """
        Aggregate the predictions for a record.
        :param frame:
        :param bin_size:
        :param num_bins:
        :return:
        """
        logger.debug('Aggregating the predictions for a batch')

        pred = np.stack(frame['tracks'].to_list(), dtype=np.float32)
        pred_seq_length = bin_size * pred.shape[2]
        agg_pred = []
        for shift_i, shift in enumerate(shifts):
            # estimate the tss bin
            # todo verify this calculation
            tss_bin = (pred_seq_length // 2 + 1 - shift) // bin_size
            # get num_bins - 1 neighboring bins centered at tss bin
            bins = [tss_bin + i for i in range(-math.floor(num_bins / 2), math.ceil(num_bins / 2))]
            assert len(bins) == num_bins
            agg_pred.append(pred[:, shift_i, bins, :])

        agg_pred = np.stack(agg_pred).swapaxes(0, 1)

        assert agg_pred.shape == (len(frame), len(shifts), num_bins, Enformer.NUM_HUMAN_TRACKS)
        # average over shifts and bins
        agg_pred = agg_pred.mean(axis=(1, 2))

        assert agg_pred.shape == (len(frame), Enformer.NUM_HUMAN_TRACKS)

        return frame.with_columns(pl.Series(values=agg_pred.tolist(), name='tracks',
                                            dtype=pl.Array(pl.Float32, width=Enformer.NUM_HUMAN_TRACKS)))


class EnformerTissueMapper:
    def __init__(self, tracks_path: str | pathlib.Path, tissue_mapper_path: str | pathlib.Path | None = None):
        """
        :param tracks_path: A yaml file mapping the name of the tracks to the index in the predictions.
        Only the tracks in the file are considered for the mapping.
        :param tissue_mapper_path: A pickle containing a dictionary of linear models for each GTEx tissue.
        """
        self.tissue_mapper_lm_dict = None
        # If tissue_mapper_path is not None, load the linear models
        if tissue_mapper_path is not None:
            with open(tissue_mapper_path, 'rb') as f:
                self.tissue_mapper_lm_dict = pickle.load(f)

        with open(tracks_path, 'rb') as f:
            self.tracks_dict = yaml.safe_load(f)

    def train(self, agg_enformer_paths: list[str] | list[pathlib.Path], expression_path: str | pathlib.Path
              , output_path: str | pathlib.Path, model=linear_model.ElasticNetCV(cv=5)):
        """
        Load the predictions from the parquet file lazily.
        For each record, calculate the average predictions over the bins centered at the tss bin.
        Collect the average predictions and train a linear model for each tissue.
        Save the linear models in a pickle file.

        :param agg_enformer_paths: The parquet files that contain the aggregated enformer predictions.
        :param expression_path: The zarr file that contains the expression scores. (ground truth)
        :param output_path: The pickle file that will contain the linear models.
        :param model: The model to use for training the tissue mapper.
        :return:
        """

        logger.info(f'Loading the expression scores from {expression_path}')
        expression_xr = xr.open_zarr(expression_path)['tpm']

        logger.info('Calculating the average expression scores...')
        expression_xr = expression_xr.groupby('subtissue').mean('sample')
        # filter out Y chromosome equivalent transcripts
        expression_xr = expression_xr.sel(transcript=~expression_xr.transcript.str.endswith('_PAR_Y'))
        transcripts = [x.split('.')[0] for x in expression_xr.transcript.values]
        expression_xr = expression_xr.assign_coords(dict(transcript=transcripts))
        tracks = list(self.tracks_dict.values())

        logger.info(f'Loading the enformer scores from {agg_enformer_paths}')
        enformer_df = pl.concat([
            pl.scan_parquet(path).select(['transcript_id', 'tracks']) for path in agg_enformer_paths
        ]).collect()
        scores = enformer_df['tracks'].to_numpy()[:, tracks]
        transcripts = enformer_df['transcript_id'].to_list()
        enformer_xr = xr.DataArray(data=scores, dims=['transcript', 'tracks'],
                                   coords=dict(transcript=transcripts, tracks=tracks), name='enformer')
        # filter out Y chromosome equivalent transcripts
        enformer_xr = enformer_xr.sel(transcript=~enformer_xr.transcript.str.endswith('_PAR_Y'))
        transcripts = [x.split('.')[0] for x in enformer_xr.transcript.values]
        enformer_xr = enformer_xr.assign_coords(dict(transcript=transcripts))

        logger.info('Merging datasets...')
        # merge the two xr data arrays
        xrds = xr.merge([expression_xr, enformer_xr], join='inner')

        logger.info('Stared training.')
        # train the linear models
        model_dict = {}
        for subtissue, subtissue_xrds in xrds.groupby('subtissue'):
            logger.info(f'Training the model for {subtissue}')
            X = subtissue_xrds['enformer'].values
            X = np.log10(1 + X)
            y = subtissue_xrds['tpm'].values
            y = np.log10(1 + y)
            lm_pipe = pipeline.Pipeline([('scaler', preprocessing.StandardScaler()),
                                         ('model', sk.clone(model))])
            lm_pipe = lm_pipe.fit(X, y)
            logger.info('Training score: %f' % lm_pipe.score(X, y))
            model_dict[subtissue] = lm_pipe

        logger.info('Saving the models...')
        self.tissue_mapper_lm_dict = model_dict
        with open(output_path, 'wb') as f:
            pickle.dump(model_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    def predict(self, agg_enformer_path: str | pathlib.Path, output_path: str | pathlib.Path):
        """
        For each tissue in the tissue_mapper_lm_dict, predict a tissue-specific expression score.
        Save the expression scores in a new parquet file.

        :param agg_enformer_path: The parquet file that contains the aggregated enformer predictions.
        :param output_path: The parquet file that will contain the tissue-specific expression scores.

        The average predictions will be calculated at the tss bin of each record.
        """
        if self.tissue_mapper_lm_dict is None:
            raise ValueError('The tissue_mapper_lm_dict is not provided. Please train the linear models first.')

        tracks = list(self.tracks_dict.values())
        logger.debug(f'Iterating over the parquet files in {agg_enformer_path}')
        enformer_df = pl.read_parquet(agg_enformer_path, hive_partitioning=False)
        scores = enformer_df['tracks'].to_numpy()[:, tracks]
        scores = np.log10(scores + 1)
        dfs = []
        for tissue, lm in self.tissue_mapper_lm_dict.items():
            tissue_df = enformer_df.select(pl.exclude('tracks'))
            if len(scores) == 0:
                res = []
            else:
                res = lm.predict(scores).flatten()
            tissue_df = tissue_df.with_columns(pl.Series(name='score', values=res, dtype=pl.Float32),
                                               pl.lit(tissue).alias('tissue'))
            dfs.append(tissue_df)
        enformer_df = pl.concat(dfs)
        enformer_df.write_parquet(output_path)


class EnformerVeff:

    def __init__(self, isoforms_path: str | pathlib.Path | None = None,
                 gtf: pd.DataFrame | str | pathlib.Path | None = None):
        """

        :param isoforms_path: The path to the file containing the isoform proportions.
        :param gtf: The path to the GTF file or a pandas DataFrame containing the genome annotation.
        """

        self.isoform_proportion_ldf = None
        if isoforms_path is not None:
            self.isoform_proportion_ldf = (pl.scan_csv(isoforms_path, separator='\t').
                                           select(['gene', 'transcript', 'tissue', 'median_transcript_proportions']).
                                           rename({'median_transcript_proportions': 'isoform_proportion',
                                                   'gene': 'gene_id', 'transcript': 'transcript_id'}).
                                           filter(~pl.col('isoform_proportion').is_null()))

        # if GTF file is given, then extract the canonical transcripts for the canonical aggregation mode
        if gtf is not None:
            if isinstance(gtf, str) or isinstance(gtf, pathlib.Path):
                gtf = gtf_to_pandas(gtf)
            elif not isinstance(gtf, pd.DataFrame):
                raise ValueError('gtf must be a path or a pandas DataFrame')

            # only keep protein_coding transcripts
            gtf = gtf.query("`gene_type` == 'protein_coding'")
            # check if Ensembl_canonical is in the set of tags
            gtf = gtf[gtf['tag'].apply(lambda x: False if pd.isna(x) else ('Ensembl_canonical' in x.split(',')))]
            self.canonical_transcripts = gtf['transcript_id'].str.extract(r'([^\.]+)\..+$')[0].unique()

    def run(self, ref_paths: list[str] | list[pathlib.Path], alt_path: str | pathlib.Path,
            output_path: str | pathlib.Path, aggregation_mode: str, upstream_tss: int | None = None,
            downstream_tss: int | None = None):
        """
        Given a file containing enformer scores for alternative alleles, calculate the variant effect.
        Then aggregate the scores by gene, variant and tissue. Save the results in a parquet file.

        :param ref_paths: The parquet files that contains the reference scores
        :param alt_path: The parquet file that contains the alternate scores
        :param output_path: The parquet file that will contain the variant effect scores
        :param aggregation_mode: One of ['logsumexp', 'weighted_sum', 'median', 'canonical'].
        :param upstream_tss: Variant effects outside of the interval [-upstream_tss, downstream_tss] will be set to 0.
        :param downstream_tss: Variant effects outside of the interval [-upstream_tss, downstream_tss] will be set to 0.
        :return:
        """

        if aggregation_mode not in ['logsumexp', 'weighted_sum', 'median', 'canonical']:
            raise ValueError(f'Unknown mode: {aggregation_mode}')
        elif aggregation_mode in ['logsumexp', 'weighted_sum']:
            assert self.isoform_proportion_ldf is not None, 'Isoform proportions are required for this mode.'
        elif aggregation_mode == 'canonical':
            assert self.canonical_transcripts is not None, 'Canonical transcripts are required for this mode.'

        logger.debug(f'Calculating the variant effect for {alt_path}')

        ref_ldf = pl.concat(
            [pl.scan_parquet(path, hive_partitioning=True).rename({'score': 'ref_score'}) for path in ref_paths])
        alt_ldf = pl.scan_parquet(alt_path).rename({'score': 'alt_score'})

        # check if alt_ldf is empty and write empty file if that's the case
        if alt_ldf.select(pl.len()).collect().item() == 0:
            logger.warning('The alternate scores are empty. No variant effect to calculate.')
            pl.DataFrame(schema={
                "chrom": pl.String,
                "strand": pl.String,
                "gene_id": pl.String,
                "variant_start": pl.Int64,
                "variant_end": pl.Int64,
                "ref": pl.String,
                "alt": pl.String,
                "tissue": pl.String,
                "veff_score": pl.Float32
            }).write_parquet(output_path)
            return

        # refactored dataloader
        seq_start_col = 'seq_start'
        seq_end_col = 'seq_end'
        if seq_start_col not in set(alt_ldf.columns):
            seq_start_col = 'enformer_start'
            seq_end_col = 'enformer_end'

        on = ['tss', 'chrom', 'strand', 'gene_id', 'transcript_id', 'transcript_start', 'transcript_end',
              seq_start_col, seq_end_col, 'tissue']

        veff_ldf = alt_ldf.join(ref_ldf, how='left', on=on)
        veff_ldf = veff_ldf.select([seq_start_col, seq_end_col, 'tss', 'chrom', 'strand',
                                    'gene_id', 'transcript_id', 'transcript_start', 'transcript_end',
                                    'variant_start', 'variant_end', 'ref', 'alt', 'tissue',
                                    'ref_score', 'alt_score'])

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
        veff_df = veff_df.rename({'log2fc': 'veff_score'})
        veff_df.write_parquet(output_path)

    def _aggregate(self, veff_ldf, aggregation_mode: str):
        """
        Given a dataframe containing variant effect scores, aggregate the scores by gene, variant and tissue.

        :param veff_ldf: A polars dataframe containing the variant effect scores.
        :param aggregation_mode: One of ['logsumexp', 'weighted_sum', 'median', 'canonical'].
        :return: A pandas dataframe containing the aggregated scores.
        """

        def logsumexp_udf(score, weight):
            return logsumexp(score / math.log10(math.e), b=weight) / math.log(10)

        if aggregation_mode in ['logsumexp', 'weighted_sum']:
            veff_ldf = veff_ldf.join(self.isoform_proportion_ldf, on=['gene_id', 'tissue', 'transcript_id'],
                                     how='inner')

            if aggregation_mode == 'logsumexp':
                veff_ldf = veff_ldf. \
                    group_by(['chrom', 'strand', 'gene_id', 'variant_start', 'variant_end', 'ref', 'alt', 'tissue']). \
                    agg(
                    pl.struct(['ref_score', 'isoform_proportion']).
                    map_elements(
                        lambda x: logsumexp_udf(x.struct.field('ref_score'), x.struct.field('isoform_proportion')),
                        return_dtype=pl.Float64()).
                    alias('ref_score'),
                    pl.struct(['alt_score', 'isoform_proportion']).
                    map_elements(
                        lambda x: logsumexp_udf(x.struct.field('alt_score'), x.struct.field('isoform_proportion')),
                        return_dtype=pl.Float64()).
                    alias('alt_score'),
                )
                veff_ldf = veff_ldf.with_columns(
                    ((pl.col("alt_score") - pl.col("ref_score")) / np.log10(2)).alias('log2fc').fill_nan(0))

            elif aggregation_mode == 'weighted_sum':
                veff_ldf = veff_ldf.with_columns(
                    (pl.col('isoform_proportion') * (pl.col("alt_score") - pl.col("ref_score")) / np.log10(2)).alias(
                        'log2fc'))
                veff_ldf = veff_ldf.group_by(['chrom', 'strand', 'gene_id', 'variant_start',
                                              'variant_end', 'ref', 'alt', 'tissue']).agg(pl.col('log2fc').sum())
            veff_df = veff_ldf.collect()
        elif aggregation_mode == 'canonical':
            # Keep only the canonical transcripts
            veff_ldf = veff_ldf.filter(pl.col('transcript_id').is_in(self.canonical_transcripts))
            veff_ldf = veff_ldf.with_columns(
                ((pl.col("alt_score") - pl.col("ref_score")) / np.log10(2)).alias('log2fc'))
            veff_ldf = veff_ldf.group_by(['chrom', 'strand', 'gene_id', 'variant_start',
                                          'variant_end', 'ref', 'alt', 'tissue', ]).agg(
                pl.col(['tss', 'transcript_id', 'transcript_start', 'transcript_end',
                        'ref_score', 'alt_score', 'log2fc']).first(),
                pl.len().alias('num_transcripts')
            )
            veff_df = veff_ldf.collect()
            # Verify that there is only one canonical transcript per gene
            max_transcripts_per_gene = veff_df['num_transcripts'].max()
            if max_transcripts_per_gene is not None and max_transcripts_per_gene > 1:
                logger.error('Multiple canonical transcripts found for a gene.')
                logger.error(veff_df[veff_df['num_transcripts'] > 1])
                raise ValueError('Multiple canonical transcripts found for a gene.')

            # Remove the num_transcripts column
            veff_df.drop_in_place('num_transcripts')
        elif aggregation_mode == 'median':
            veff_ldf = veff_ldf.with_columns(
                ((pl.col("alt_score") - pl.col("ref_score")) / np.log10(2)).alias('log2fc'))
            veff_ldf = veff_ldf.group_by(['chrom', 'strand', 'gene_id', 'variant_start',
                                          'variant_end', 'ref', 'alt', 'tissue', ]).agg(pl.col('log2fc').median())
            veff_df = veff_ldf.collect()
        else:
            raise ValueError(f'Unknown mode: {aggregation_mode}')

        logger.debug(f'Aggregated table size: {len(veff_df)}')
        return veff_df
