from absplice import SplicingOutlierResult
import pandas as pd

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

print("computing result...", flush=True)
splicing_result = SplicingOutlierResult(
    df_mmsplice=df_mmsplice, 
    df_spliceai=df_spliceai,
)
df = splicing_result.predict_absplice_dna()

print("cleaning up result...", flush=True)
# make variant[str] a normal column
out_df = df.reset_index()

# split variant[str] into chrom, pos, ref, alt columns
variant_idx = pd.DataFrame.from_records(
    out_df["variant"].str.split(":|>"),
    columns=["chrom", "pos", "ref", "alt"]
)
variant_idx = variant_idx.astype({"pos": "int"})

# assign to output df and add start and end columns
out_df = out_df.assign(**variant_idx.to_dict())
out_df = out_df.assign(**{
    "start": out_df["pos"] - 1,
    "end": out_df["pos"] - 1 + out_df["ref"].str.len(),
})

print(f'''writing result to "{snakemake.output['absplice_dna']}"...''', flush=True)
# drop variant column and save
(
    out_df
    .drop(columns=["variant"])
    .to_parquet(snakemake.output['absplice_dna'], index=False)
)

print('All done!', flush=True)