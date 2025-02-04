# %%
from absplice import SplicingOutlierResult
import pandas as pd

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
snakefile_path = os.getcwd() + "/../../Snakefile"

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
df_mmsplice=pd.read_csv(
    snakemake.input['mmsplice_splicemap'],
    usecols=[
      'variant',
      'tissue',
      'junction',
      'event_type',
      'splice_site',
      'gene_id',
      'gene_name',
      'delta_logit_psi',
      'ref_psi',
      'median_n',
      'delta_psi'
    ],
    dtype={
        'variant': pd.StringDtype("pyarrow"),
        'gene_id': pd.StringDtype("pyarrow"),
        'gene_name': pd.StringDtype("pyarrow"),
        'event_type': pd.StringDtype("pyarrow"),
        'splice_site': pd.StringDtype("pyarrow"),
        'tissue': pd.StringDtype("pyarrow"),
        'junction': pd.StringDtype("pyarrow"),
        'delta_logit_psi': 'float64',
        'delta_psi': 'float64',
        'ref_psi': 'float64',
        'median_n': 'float64',
    }
)
print(f'''loading df_spliceai from "{snakemake.input['spliceai']}"...''', flush=True)
df_spliceai=pd.read_csv(
    snakemake.input['spliceai'],
    usecols=[
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
    dtype={
        'variant': pd.StringDtype("pyarrow"),
        'gene_name': pd.StringDtype("pyarrow"),
        'delta_score': 'float32',
        'acceptor_gain': 'float32',
        'acceptor_loss': 'float32',
        'donor_gain': 'float32',
        'donor_loss': 'float32',
        'acceptor_gain_position': 'int32',
        'acceptor_loss_position': 'int32',
        'donor_gain_position': 'int32',
        'donor_loss_position': 'int32',
    }
)

# %%
fake_chrom = pd.read_csv(snakemake.input["chrom_alias"], sep="\t")["chrom"].unique()[0] + "__FAKE__"
fake_variant = f"{fake_chrom}:1-2:A>B"
print(f"Fake variant: '{fake_variant}'", flush=True)

tissues = pd.read_csv(snakemake.input['tissue_mapping'])["tissue"]

fake_df_mmsplice = pd.DataFrame.from_dict({
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
}).astype(df_mmsplice.dtypes)

# %%
if df_mmsplice.empty:
    print("MMSplice dataframe empty. Faking mmsplice output...")
    df_mmsplice = fake_df_mmsplice

# %%
df_spliceai = df_spliceai.set_index("variant")
df_spliceai

# %%
df_mmsplice = df_mmsplice.set_index("variant")
df_mmsplice

# %%
all_variants = df_spliceai.index.unique().union(df_mmsplice.index.unique()).sort_values()
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

        batch_df_mmsplice = df_mmsplice.iloc[df_mmsplice.index.isin(batch)].reset_index()
        if batch_df_mmsplice.empty:
            batch_df_mmsplice = fake_df_mmsplice

        batch_df_spliceai = df_spliceai.iloc[df_spliceai.index.isin(batch)].reset_index()

        splicing_result = SplicingOutlierResult(
            df_mmsplice=batch_df_mmsplice, 
            df_spliceai=batch_df_spliceai,
        )
        df = splicing_result.predict_absplice_dna()
        # result_dfs.append(df)

        if fake_variant is not None:
            # print("Dropping fake variants...")
            df = df.query(f"`variant` != '{fake_variant}'")

        # make variant[str] a normal column
        out_df = df.reset_index()
        # split variant[str] into chrom, pos, ref, alt columns
        variant_idx = pd.DataFrame.from_records(
            out_df["variant"].str.split(":|>"),
            columns=["chrom", "pos", "ref", "alt"]
        )
        variant_idx = variant_idx.astype({"pos": "int"})

        # assign to output df and add start and end columns
        out_df = out_df.assign(**variant_idx.to_dict(orient='series'))
        out_df = out_df.assign(**{
            "start": out_df["pos"] - 1,
            "end": out_df["pos"] - 1 + out_df["ref"].str.len(),
        })
        out_df = out_df.drop(columns=["variant"])

        # write to parquet
        table = pa.Table.from_pandas(out_df, schema=output_schema, preserve_index=False)
        pqwriter.write_table(table)

        # cleanup memory explicitly
        gc.collect()

# %%
snakemake.output['absplice_dna']

# %%
print('All done!', flush=True)

# %%
