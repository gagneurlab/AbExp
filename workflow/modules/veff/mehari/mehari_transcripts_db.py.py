# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:anaconda-abexp-veff-py]
#     language: python
#     name: conda-env-anaconda-abexp-veff-py-py
# ---

# %% [markdown]
# # Build the mehari transcript database from GENCODE
#
# Input: a GENCODE GFF3 annotation and the matching transcript FASTA. The script
#
# 1. cuts the `|`-joined GENCODE FASTA headers down to the transcript id (mehari matches
#    sequences on the first whitespace-separated token of the header),
# 2. collects the transcript tags mehari knows (basic, Ensembl_canonical, MANE_Select,
#    MANE_Plus_Clinical) into the `mane_transcripts` TSV,
# 3. reads the GENCODE and Ensembl release from the GFF3 header,
# 4. calls `mehari.build_transcript_db`.

# %%
from IPython.display import display

# %% jupyter={"outputs_hidden": false} pycharm={"name": "#%%\n"}
import os
import sys
import gzip
import re
import inspect
import logging

import yaml

from pprint import pprint

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
        rule_name = 'veff__mehari_transcripts_db',
        default_wildcards={}
    )

# %%
try:
    from snakemk_util import pretty_print_snakemake
    print(pretty_print_snakemake(snakemake))
except ImportError:
    # not installed in the mehari conda environment
    pass

# %%
os.getcwd()

# %%
# mehari logs through Python's `logging`; show its progress messages
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s: %(message)s")

# %%
from tqdm.auto import tqdm
from mehari import build_transcript_db, SeqvarsAnnotator


# %%
def open_text(path):
    """Opens a plain or gzip-compressed text file for reading."""
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rt")
    return open(path, "rt")


# %% [markdown]
# # Transcript FASTA
#
# GENCODE headers look like `>ENST00000456328.2|ENSG00000290825.1|-|-|DDX11L2-202|DDX11L2|1657|lncRNA|`.
# mehari takes the first whitespace-separated token as the sequence id, so keep only the
# first `|` field. The output is uncompressed: mehari reads `.gz` FASTA files as BGZF only.

# %%
n_sequences = 0
with open_text(snakemake.input["transcripts_fasta"]) as fin, open(snakemake.output["transcripts_fasta"], "w") as fout:
    for line in fin:
        if line.startswith(">"):
            n_sequences += 1
            line = ">" + line[1:].rstrip("\n").split("|", 1)[0] + "\n"
        fout.write(line)
print(f"{n_sequences} transcript sequences")

# %% [markdown]
# # Transcript tags
#
# One row per transcript with the tags mehari understands, in the `mane_transcripts` TSV
# layout: transcript id without version, version, gene name, comma-separated labels.

# %%
KEPT_TAGS = {"basic", "Ensembl_canonical", "MANE_Select", "MANE_Plus_Clinical"}


# %%
def gff3_attributes(field):
    return dict(kv.split("=", 1) for kv in field.rstrip("\n").split(";") if "=" in kv)


# %%
gff3_path = snakemake.input["gff3"]
assert re.search(r"\.gff3?(\.gz)?$", gff3_path), \
    f"mehari recognizes GFF3 files by the extension .gff3, .gff, .gff3.gz or .gff.gz, got '{gff3_path}'"

# %%
header_lines = []
n_transcripts = 0
n_tagged = 0
with open_text(gff3_path) as fin, open(snakemake.output["tags_tsv"], "w") as fout:
    for line in fin:
        if line.startswith("#"):
            if len(header_lines) < 20:
                header_lines.append(line.rstrip("\n"))
            continue
        fields = line.split("\t")
        if len(fields) < 9 or fields[2] != "transcript":
            continue
        n_transcripts += 1
        attributes = gff3_attributes(fields[8])
        transcript_id = attributes.get("transcript_id")
        if transcript_id is None:
            continue
        tags = [tag for tag in attributes.get("tag", "").split(",") if tag in KEPT_TAGS]
        if not tags:
            continue
        n_tagged += 1
        transcript, _, version = transcript_id.partition(".")
        fout.write(f"{transcript}\t{version}\t{attributes.get('gene_name', '')}\t{','.join(tags)}\n")
print(f"{n_transcripts} transcripts in the GFF3, {n_tagged} with tags")

# %% [markdown]
# # Release information from the GFF3 header
#
# e.g. `#description: evidence-based annotation of the human genome (GRCh38), version 42 (Ensembl 108)`.
# mehari requires the Ensembl release for Ensembl transcripts. When the header does not name
# it, human GENCODE release N corresponds to Ensembl release N + 66.

# %%
header = "\n".join(header_lines)
print(header)

# %%
# the release is on the `#description:` line; `##gff-version 3` must not match
match = re.search(r"^#description:.*version (\d+)", header, re.MULTILINE)
gencode_release = match.group(1) if match else None
match = re.search(r"Ensembl (\d+)", header)
ensembl_release = match.group(1) if match else None
if ensembl_release is None and gencode_release is not None:
    ensembl_release = str(int(gencode_release) + 66)
match = re.search(r"GRCh3[78]\.p\d+", header)
assembly_version = match.group(0) if match else None

print(f"GENCODE release: {gencode_release}, Ensembl release: {ensembl_release}, assembly version: {assembly_version}")

# %%
assert ensembl_release is not None, "Could not read the Ensembl release from the GFF3 header"

# %% [markdown]
# # Build the database
#
# mehari compresses the output only when the file name ends in `.zst`. The database is
# written under a temporary name first and moved into place after a successful build.

# %%
output_db = snakemake.output["transcripts_db"]
assert output_db.endswith(".zst"), f"transcripts_db must end in .zst, got '{output_db}'"
building_db = output_db[:-len(".zst")] + ".building.zst"

# %%
build_kwargs = dict(
    assembly=snakemake.params["assembly"],
    annotation=[gff3_path],
    transcript_sequences=snakemake.output["transcripts_fasta"],
    transcript_source="Ensembl",
    transcript_source_version=ensembl_release,
    annotation_version=f"GENCODE {gencode_release}" if gencode_release is not None else None,
    assembly_version=assembly_version,
    mane_transcripts=snakemake.output["tags_tsv"],
    threads=snakemake.threads,
    output=building_db,
)
# progress bars need a mehari build that has the `progress` argument
if "progress" in inspect.signature(build_transcript_db).parameters:
    build_kwargs["progress"] = tqdm
pprint(build_kwargs)

# %%
build_transcript_db(**build_kwargs)

# %%
os.replace(building_db, output_db)
if os.path.exists(building_db + ".report.jsonl"):
    os.replace(building_db + ".report.jsonl", output_db + ".report.jsonl")
print(f"wrote {output_db} ({os.path.getsize(output_db) / 1e6:.1f} MB)")

# %% [markdown]
# # Check that mehari loads the database

# %%
SeqvarsAnnotator(transcript_db_paths=[output_db])
print("database loads")

# %%
