# GENCODE files for the mehari transcript database

tl;dr: run `download_gencode.sh` to fetch the GENCODE GFF3 annotation and transcript
FASTA of one release, then point `gff3_file` and `veff.mehari_gencode_transcripts_fasta` in
`config/config.yaml` at the two files. The pipeline rule `veff__mehari_transcripts_db`
(`workflow/modules/veff/mehari/mehari_transcripts_db.py.py`) builds the mehari transcript
database from them.

## Command

```
misc/mehari/download_gencode.sh <gencode_release> <GRCh38|GRCh37> <output_dir>
```

Example, GENCODE 42 (= Ensembl 108) for GRCh38:

```
misc/mehari/download_gencode.sh 42 GRCh38 data/gencode/release_42
```

Use the GENCODE release of your `gff3_file`. Human GENCODE release N corresponds to
Ensembl release N + 66.

## mehari version requirement

mehari 0.45.1 imports GENCODE GFF3 files incorrectly and pads incomplete CDS on
minus-strand transcripts wrong (fixed upstream in pull request #1062). Both are fixed
in the 0.46.0 release, which `workflow/modules/veff/envs/veff-py.yaml`
installs from bioconda.

## GRCh37

The GRCh37 route uses GENCODE's lift-over files and is untested end to end; treat it
as a starting point, not a verified path.
