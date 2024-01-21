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

import textwrap

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
            "vcf_file": "chrom=chr3/178395969-178420969.vcf.gz",
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
fake_variant = None

# %%
if df_mmsplice.empty:
    print("Faking mmsplice output...")
    fake_chrom = pd.read_csv(snakemake.input["chrom_alias"], sep="\t")["chrom"].unique()[0] + "__FAKE__"
    fake_variant = f"{fake_chrom}:1-2:A>B"
    print(f"Fake variant: '{fake_variant}'", flush=True)
    
    tissues = pd.read_csv(snakemake.input['tissue_mapping'])["tissue"]
    
    df_mmsplice = pd.DataFrame.from_dict({
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
print("computing result...", flush=True)
splicing_result = SplicingOutlierResult(
    df_mmsplice=df_mmsplice, 
    df_spliceai=df_spliceai,
)
df = splicing_result.predict_absplice_dna()

# %%
if fake_variant is not None:
    print("Dropping fake variants...")
    df = df.query(f"`variant` != '{fake_variant}'")

# %%
print("cleaning up result...", flush=True)
# make variant[str] a normal column
out_df = df.reset_index()

# %%
# split variant[str] into chrom, pos, ref, alt columns
variant_idx = pd.DataFrame.from_records(
    out_df["variant"].str.split(":|>"),
    columns=["chrom", "pos", "ref", "alt"]
)
variant_idx = variant_idx.astype({"pos": "int"})

# %%
# assign to output df and add start and end columns
out_df = out_df.assign(**variant_idx.to_dict())
out_df = out_df.assign(**{
    "start": out_df["pos"] - 1,
    "end": out_df["pos"] - 1 + out_df["ref"].str.len(),
})

# %%
print(f'''writing result to "{snakemake.output['absplice_dna']}"...''', flush=True)
# drop variant column and save
(
    out_df
    .drop(columns=["variant"])
    .to_parquet(snakemake.output['absplice_dna'], index=False)
)

# %%
print('All done!', flush=True)

# %%
