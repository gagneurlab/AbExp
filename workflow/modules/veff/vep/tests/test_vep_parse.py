"""Tests of the rule veff__vep_parse

Each test runs `vep_parse.py.py` on a tiny VEP table of one variant in one transcript, and checks
the key columns and the columns that the parser derives from Consequence and IMPACT. The other
columns of the row (about 40 of the VEP output) are passed through or only cast, so a test of them
would only repeat the dtypes of the script. Run the tests in the environment of
`workflow/modules/veff/envs/veff-py.yaml`, with pytest added:

    pytest workflow/modules/veff/vep/tests
"""
import runpy
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import polars as pl

SCRIPT = Path(__file__).parent.parent / "vep_parse.py.py"

# the columns that the script reads, in the order of the VEP --tab output
HEADER = [
    "#Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "IMPACT",
    "CANONICAL",
    "GENE_PHENO",
    "SIFT",
    "PolyPhen",
    "Condel",
]

GENE = "ENSG00000000001"
TRANSCRIPT = "ENST00000000001"


@dataclass(frozen=True)
class VepRow:
    """A row of the VEP table: a deletion of 9 bases after the anchor base, which covers the transcript"""

    consequence: str
    impact: str
    gene: str = GENE
    transcript: str = TRANSCRIPT

    def to_line(self) -> str:
        columns = [
            "chr1:100:110:ACGTACGTAC>A",
            "chr1:102-110",
            "-",
            self.gene,
            self.transcript,
            "Transcript",
            self.consequence,
            self.impact,
            "YES",
            "-",
            "-",
            "-",
            "-",
        ]
        return "\t".join(columns)


def run_vep_parse(tmp_path: Path, rows: list[VepRow]) -> list[dict]:
    """Writes the input files, runs the script and returns the rows of its parquet file"""
    chrom_alias_path = tmp_path / "chrom_alias.tsv"
    chrom_alias_path.write_text("#alias\tchrom\nchr1\tchr1\n")

    # VEP writes ## lines first, then the column names with the prefix "#"
    veff_header_path = tmp_path / "veff.tsv.header"
    veff_header_path.write_text("## ENSEMBL VARIANT EFFECT PREDICTOR\n" + "\t".join(HEADER) + "\n")

    veff_tsv_path = tmp_path / "veff.tsv"
    veff_tsv_path.write_text("".join(row.to_line() + "\n" for row in rows))

    snakemake = SimpleNamespace(
        input={
            "chrom_alias": str(chrom_alias_path),
            "veff_header": str(veff_header_path),
            "veff_tsv": str(veff_tsv_path),
        },
        output={"veff_pq": str(tmp_path / "veff.parquet")},
        params={},
        config={},
    )

    runpy.run_path(str(SCRIPT), init_globals={"snakemake": snakemake})

    return pl.read_parquet(snakemake.output["veff_pq"]).to_dicts()


def true_terms(struct: dict) -> set[str]:
    """The names of the fields of a binarized struct that are True"""
    return {name for name, value in struct.items() if value}


# VEP 116 gives a variant that covers a whole protein-coding transcript the term
# coding_transcript_variant. tissue_specific_vep writes a feature for each field of the struct
# Consequence, so the parser must keep the term as a field.
def test_coding_transcript_variant_is_kept(tmp_path):
    rows = [VepRow(consequence="coding_transcript_variant", impact="MODIFIER")]

    parsed = run_vep_parse(tmp_path, rows)

    assert len(parsed) == 1
    row = parsed[0]
    assert "coding_transcript_variant" in row["Consequence"]
    assert {
        "chrom": row["chrom"],
        "start": row["start"],
        "end": row["end"],
        "ref": row["ref"],
        "alt": row["alt"],
        "pos": row["pos"],
        "len": row["len"],
        "gene": row["gene"],
        "transcript": row["transcript"],
        "CANONICAL": row["CANONICAL"],
        "Consequence": true_terms(row["Consequence"]),
        "IMPACT": true_terms(row["IMPACT"]),
    } == {
        "chrom": "chr1",
        "start": 100,
        "end": 110,
        "ref": "ACGTACGTAC",
        "alt": "A",
        "pos": 101,
        "len": 10,
        "gene": GENE,
        "transcript": TRANSCRIPT,
        "CANONICAL": True,
        "Consequence": {"coding_transcript_variant"},
        "IMPACT": {"MODIFIER"},
    }
