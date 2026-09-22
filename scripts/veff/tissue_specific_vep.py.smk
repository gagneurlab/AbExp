
# consequence table from the configured annotator (vep.smk or mehari.smk)
VEP_PQ_INPUT_PATTERN=f"{VEFF_BASEDIR}/{VEFF_ANNOTATOR}/veff.parquet/{{vcf_file}}.parquet"

OUTPUT_BASEDIR=f"{VEFF_BASEDIR}/tissue_specific_vep.py"
VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"


rule veff__tissue_specific_vep:
    threads: lambda wildcards, attempt: 16 * attempt,
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (4000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        vep_pq=VEP_PQ_INPUT_PATTERN,
        isoform_proportions_pq=config["system"]["isoform_proportions_pq"],
        gtf_transcripts=f"{RESULTS_DIR}/gtf_transcripts.parquet",
        chrom_alias=ancient(CHROM_ALIAS_TSV),
    wildcard_constraints:
        ds_dir="[^/]+",
    conda:
        "../../envs/abexp-veff-py.yaml"
    script:
        "tissue_specific_vep.py.py"


del VEP_PQ_INPUT_PATTERN

del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN
