from absplice import SplicingOutlierResult
import pandas as pd

splicing_result = SplicingOutlierResult(
        df_mmsplice=snakemake.input['mmsplice_splicemap'], 
        df_spliceai=snakemake.input['spliceai'],
    )
df = splicing_result.predict_absplice_dna()


out_df = df.reset_index()

variant_idx = pd.DataFrame.from_records(
  out_df["variant"].str.split(":|>"),
  columns=["chrom", "pos", "ref", "alt"]
)
variant_idx = variant_idx.astype({"pos": "int"})

out_df = out_df.assign(**variant_idx.to_dict())
out_df = out_df.assign(**{
  "start": out_df["pos"] - 1,
  "end": out_df["pos"] - 1 + out_df["ref"].str.len(),
})

df.reset_index().to_parquet(snakemake.output['absplice_dna'], index=False)