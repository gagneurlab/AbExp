import os

# mehari (https://github.com/varfish-org/mehari) as an alternative to VEP for the
# transcript consequence annotation. Selected with `veff.annotator: "mehari"` in the
# config. The annotation step calls the mehari Python package and writes mehari's
# per-transcript annotation with mehari's own consequence terms; `tissue_specific_vep.py.smk`
# aggregates it like the VEP table. The resulting features differ from the VEP features, so
# the shipped AbExp models need VEP; the mehari route is for training new models.
#
# The mehari transcript database is built by `veff__mehari_transcripts_db` from the
# GENCODE GFF3 and transcript FASTA given in the config (`veff.mehari_gencode_gff3`,
# `veff.mehari_gencode_transcripts_fasta`), of the same GENCODE release as `gtf_file`.
#
# LoF and NMD calls are not part of this step; LOFTEE and NMD-Scanner get their own steps.

OUTPUT_BASEDIR=f"{VEFF_BASEDIR}/mehari"

VEFF_VCF_PQ_PATTERN=f"{OUTPUT_BASEDIR}/veff.parquet/{{vcf_file}}.parquet"

MEHARI_TRANSCRIPTS_DB=config["system"]["mehari"]["transcripts_db"].format(
    human_genome_assembly=ASSEMBLY,
)
# conda environment of the mehari rules: `mehari_env.yaml` (built by `snakemake --sdm conda`),
# or the name of an existing environment with the mehari Python package (`mehari.conda_env`)
MEHARI_CONDA_ENV=config["system"]["mehari"]["conda_env"] or "mehari_env.yaml"

# the config schema requires both or none of them, and both for `veff.annotator: "mehari"`
MEHARI_GENCODE_GFF3 = config["veff"].get("mehari_gencode_gff3")
MEHARI_GENCODE_TRANSCRIPTS_FASTA = config["veff"].get("mehari_gencode_transcripts_fasta")


if MEHARI_GENCODE_GFF3 is not None:
    rule veff__mehari_transcripts_db:
        """
        Builds the mehari transcript database from a GENCODE GFF3 annotation and the
        matching transcript FASTA.
        """
        threads: 4
        resources:
            ntasks=1,
            mem_mb=lambda wildcards, attempt, threads: 8000 * attempt,
        output:
            transcripts_db=MEHARI_TRANSCRIPTS_DB,
            # transcript FASTA with the headers cut down to the transcript id
            transcripts_fasta=temp(f"{MEHARI_TRANSCRIPTS_DB}.transcripts.fa"),
            # transcript tags (basic, Ensembl_canonical, MANE_Select, MANE_Plus_Clinical)
            tags_tsv=f"{MEHARI_TRANSCRIPTS_DB}.tags.tsv",
        input:
            gff3=os.path.abspath(MEHARI_GENCODE_GFF3),
            transcripts_fasta=os.path.abspath(MEHARI_GENCODE_TRANSCRIPTS_FASTA),
        params:
            assembly=ASSEMBLY.lower(),
        conda:
            MEHARI_CONDA_ENV
        script:
            "mehari_transcripts_db.py.py"


rule veff__mehari_annotation:
    """
    Annotates the variants of one VCF with mehari and writes mehari's per-transcript table.
    """
    threads: 4
    resources:
        ntasks=1,
        mem_mb=lambda wildcards, attempt, threads: (2000 * threads) * attempt,
    output:
        veff_pq=VEFF_VCF_PQ_PATTERN,
    input:
        chrom_alias=ancient(CHROM_ALIAS_TSV),
        vcf=STRIPPED_VCF_FILE_PATTERN,
        transcripts_db=MEHARI_TRANSCRIPTS_DB,
        fasta=FASTA_FILE,
        fasta_index=FASTA_INDEX_FILE,
    wildcard_constraints:
        ds_dir="[^/]+",
        feature_set="[^/]+",
    conda:
        MEHARI_CONDA_ENV
    script:
        "mehari_annotation.py.py"


del OUTPUT_BASEDIR
del VEFF_VCF_PQ_PATTERN

del MEHARI_TRANSCRIPTS_DB
del MEHARI_CONDA_ENV
del MEHARI_GENCODE_GFF3
del MEHARI_GENCODE_TRANSCRIPTS_FASTA
