# %%
from abexp.absplice import SplicingOutlierResult
import polars as pl

# %%
import os
import sys
import shutil

import json
import yaml

from pprint import pprint

from tqdm.auto import tqdm

import textwrap

# %%
import pyarrow as pa
import pyarrow.parquet as pq

# %%
import gc

# %%
snakefile_path = os.getcwd() + "/../../../Snakefile"

# %%
# del snakemake

# %%
try:
    snakemake
except NameError:
    from snakemk_util import load_rule_args
    
    snakemake = load_rule_args(
        snakefile = snakefile_path,
        rule_name = 'absplice_dna',
        default_wildcards={
            "vcf_file": "chrom=chr5/140871637-140899832.vcf.gz",
        }
    )

# %%
print(f'''loading df_mmsplice from "{snakemake.input['mmsplice_splicemap']}"...''', flush=True)
# the columns in the order of the file; the order breaks ties of the strongest junction
df_mmsplice = pl.read_csv(
    snakemake.input['mmsplice_splicemap'],
    columns=[
        'variant',
        'tissue',
        'junction',
        'event_type',
        'splice_site',
        'ref_psi',
        'median_n',
        'gene_id',
        'gene_name',
        'delta_logit_psi',
        'delta_psi',
    ],
    schema_overrides={
        'variant': pl.String,
        'gene_id': pl.String,
        'gene_name': pl.String,
        'event_type': pl.String,
        'splice_site': pl.String,
        'tissue': pl.String,
        'junction': pl.String,
        'delta_logit_psi': pl.Float64,
        'delta_psi': pl.Float64,
        'ref_psi': pl.Float64,
        'median_n': pl.Float64,
    }
)
print(f'''loading df_spliceai from "{snakemake.input['spliceai']}"...''', flush=True)
df_spliceai = pl.read_csv(
    snakemake.input['spliceai'],
    columns=[
        'variant',
        'gene_name',
        'delta_score',
        'acceptor_gain',
        'acceptor_loss',
        'donor_gain',
        'donor_loss',
        'acceptor_gain_position',
        'acceptor_loss_position',
        'donor_gain_position',
        'donor_loss_position',
    ],
    schema_overrides={
        'variant': pl.String,
        'gene_name': pl.String,
        'delta_score': pl.Float32,
        'acceptor_gain': pl.Float32,
        'acceptor_loss': pl.Float32,
        'donor_gain': pl.Float32,
        'donor_loss': pl.Float32,
        'acceptor_gain_position': pl.Int32,
        'acceptor_loss_position': pl.Int32,
        'donor_gain_position': pl.Int32,
        'donor_loss_position': pl.Int32,
    }
)

# %%
fake_chrom = pl.read_csv(snakemake.input["chrom_alias"], separator="\t")["chrom"][0] + "__FAKE__"
fake_variant = f"{fake_chrom}:1-2:A>B"
print(f"Fake variant: '{fake_variant}'", flush=True)

tissues = pl.read_csv(snakemake.input['tissue_mapping'])["tissue"]

fake_df_mmsplice = pl.DataFrame({
    'variant': fake_variant,
    'tissue': tissues,
    'junction': "",
    'event_type': "",
    'splice_site': "",
    'ref_psi': 0,
    'median_n': 0,
    'gene_id': "",
    'gene_name': "",
    'delta_logit_psi': 0,
    'delta_psi': 0,
}).cast(df_mmsplice.schema)

# %%
if df_mmsplice.is_empty():
    print("MMSplice dataframe empty. Faking mmsplice output...")
    df_mmsplice = fake_df_mmsplice

# %%
all_variants = pl.concat([df_spliceai["variant"], df_mmsplice["variant"]]).unique().sort()
all_variants

# %%
output_schema = pa.schema([
    pa.field('chrom', pa.string()),
    pa.field('start', pa.int64()),
    pa.field('end', pa.int64()),
    pa.field('ref', pa.string()),
    pa.field('alt', pa.string()),
    pa.field('pos', pa.int64()),
    pa.field('gene_id', pa.string()),
    pa.field('tissue', pa.string()),
    pa.field('delta_logit_psi', pa.float64()),
    pa.field('delta_psi', pa.float64()),
    pa.field('delta_score', pa.float64()),
    pa.field('splice_site_is_expressed', pa.int64()),
    pa.field('AbSplice_DNA', pa.float32()),
    pa.field('junction', pa.string()),
    pa.field('event_type', pa.string()),
    pa.field('splice_site', pa.string()),
    pa.field('ref_psi', pa.float64()),
    pa.field('median_n', pa.float64()),
    pa.field('acceptor_gain', pa.float32()),
    pa.field('acceptor_loss', pa.float32()),
    pa.field('donor_gain', pa.float32()),
    pa.field('donor_loss', pa.float32()),
    pa.field('acceptor_gain_position', pa.int32()),
    pa.field('acceptor_loss_position', pa.int32()),
    pa.field('donor_gain_position', pa.int32()),
    pa.field('donor_loss_position', pa.int32()),
])
output_schema

# %%
print("computing result...", flush=True)

batch_size = snakemake.params["variants_per_batch"]

with pq.ParquetWriter(snakemake.output['absplice_dna'], output_schema) as pqwriter:
    for i in tqdm(range(0, len(all_variants), batch_size)):
        batch = all_variants[i:i+batch_size]

        batch_df_mmsplice = df_mmsplice.filter(pl.col("variant").is_in(batch.implode()))
        if batch_df_mmsplice.is_empty():
            batch_df_mmsplice = fake_df_mmsplice

        batch_df_spliceai = df_spliceai.filter(pl.col("variant").is_in(batch.implode()))

        splicing_result = SplicingOutlierResult(
            df_mmsplice=batch_df_mmsplice, 
            df_spliceai=batch_df_spliceai,
        )
        df = splicing_result.predict_absplice_dna()
        # result_dfs.append(df)

        if fake_variant is not None:
            # print("Dropping fake variants...")
            df = df.filter(pl.col("variant") != fake_variant)

        # split variant[str] into chrom, pos, ref, alt columns, and add start and end columns
        out_df = df \
            .with_columns(pl.col("variant").str.extract_groups(
                r"^(?P<chrom>[^:]+):(?P<pos>\d+):(?P<ref>[^>]+)>(?P<alt>.+)$")) \
            .unnest("variant") \
            .with_columns(pl.col("pos").cast(pl.Int64)) \
            .with_columns(start=pl.col("pos") - 1, end=pl.col("pos") - 1 + pl.col("ref").str.len_chars())

        # write to parquet
        table = out_df.select(output_schema.names).to_arrow().cast(output_schema)
        pqwriter.write_table(table)

        # cleanup memory explicitly
        gc.collect()

# %%
snakemake.output['absplice_dna']

# %%
print('All done!', flush=True)

# %%
