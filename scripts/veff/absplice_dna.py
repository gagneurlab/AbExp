from absplice import SplicingOutlierResult
import pandas as pd

splicing_result = SplicingOutlierResult(
        df_mmsplice=snakemake.input['mmsplice_splicemap'], 
        df_spliceai=snakemake.input['spliceai'],
    )
df = splicing_result.predict_absplice_dna()

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

# drop variant column and save
(
  out_df
  .drop(columns=["variant"])
  .to_parquet(snakemake.output['absplice_dna'], index=False)
)