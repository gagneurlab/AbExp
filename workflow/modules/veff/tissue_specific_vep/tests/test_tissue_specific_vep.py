"""Tests of the rule veff__tissue_specific_vep

Each test runs `tissue_specific_vep.py.py` on small input tables of one variant in one gene, and
checks the full output row. Run them in the environment of
`workflow/modules/veff/envs/veff-py.yaml`, with pytest added:

    pytest workflow/modules/veff/tissue_specific_vep/tests
"""
import runpy
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import polars as pl
import pytest

SCRIPT = Path(__file__).parent.parent / "tissue_specific_vep.py.py"

TISSUE = "Whole Blood"

# A gene in a pseudoautosomal region (PAR), and a gene elsewhere. Each has two transcripts.
PAR_GENE = "ENSG00000000001"
NON_PAR_GENE = "ENSG00000000002"


@dataclass(frozen=True)
class GtfTranscript:
    """A row of the table of the gtf_transcripts module"""

    chromosome: str
    gene_id: str
    transcript_id: str


# As in GENCODE, the chrY copies of the PAR transcripts have the same IDs as the chrX copies, with
# the suffix "_PAR_Y".
GTF_TRANSCRIPTS = [
    GtfTranscript("chrX", "ENSG00000000001.1", "ENST00000000001.1"),
    GtfTranscript("chrX", "ENSG00000000001.1", "ENST00000000002.1"),
    GtfTranscript("chrY", "ENSG00000000001.1_PAR_Y", "ENST00000000001.1_PAR_Y"),
    GtfTranscript("chrY", "ENSG00000000001.1_PAR_Y", "ENST00000000002.1_PAR_Y"),
    GtfTranscript("chr22", "ENSG00000000002.1", "ENST00000000003.1"),
    GtfTranscript("chr22", "ENSG00000000002.1", "ENST00000000004.1"),
]


@dataclass(frozen=True)
class IsoformProportion:
    """A row of the isoform proportion table, in TISSUE"""

    gene: str
    transcript: str
    median_transcript_proportions: float


# The medians of each gene add up to 1 in the tissue.
ISOFORM_PROPORTIONS = [
    IsoformProportion(PAR_GENE, "ENST00000000001", 0.75),
    IsoformProportion(PAR_GENE, "ENST00000000002", 0.25),
    IsoformProportion(NON_PAR_GENE, "ENST00000000003", 0.75),
    IsoformProportion(NON_PAR_GENE, "ENST00000000004", 0.25),
]


@dataclass(frozen=True)
class Consequence:
    """A row of the consequence table of the vep or mehari module"""

    transcript: str
    stop_gained: bool
    missense_variant: bool


def run_tissue_specific_vep(tmp_path: Path, chrom: str, gene: str, consequences: list[Consequence]) -> list[dict]:
    """Writes the input tables for a variant at `chrom`:100 in `gene`, runs the script and returns
    its rows, with the struct `features` unnested."""
    vep_df = pl.DataFrame(
        [
            {
                "chrom": chrom,
                "start": 100,
                "end": 101,
                "ref": "C",
                "alt": "T",
                "gene": gene,
                "transcript": c.transcript,
                "Consequence": {"stop_gained": c.stop_gained, "missense_variant": c.missense_variant},
            }
            for c in consequences
        ],
        schema={
            "chrom": pl.String,
            "start": pl.Int64,
            "end": pl.Int64,
            "ref": pl.String,
            "alt": pl.String,
            "gene": pl.String,
            "transcript": pl.String,
            "Consequence": pl.Struct({"stop_gained": pl.Boolean, "missense_variant": pl.Boolean}),
        },
    )
    gtf_transcripts_df = pl.DataFrame(
        [
            {
                "Chromosome": t.chromosome,
                "gene_id": t.gene_id,
                "transcript_id": t.transcript_id,
                "transcript_biotype": "protein_coding",
            }
            for t in GTF_TRANSCRIPTS
        ],
        schema={
            "Chromosome": pl.String,
            "gene_id": pl.String,
            "transcript_id": pl.String,
            "transcript_biotype": pl.String,
        },
    )
    isoform_proportions_df = pl.DataFrame(
        [
            {
                "gene": p.gene,
                "tissue": TISSUE,
                "transcript": p.transcript,
                "mean_transcript_proportions": p.median_transcript_proportions,
                "median_transcript_proportions": p.median_transcript_proportions,
                "sd_transcript_proportions": 0.0,
            }
            for p in ISOFORM_PROPORTIONS
        ],
        schema={
            "gene": pl.String,
            "tissue": pl.String,
            "transcript": pl.String,
            "mean_transcript_proportions": pl.Float32,
            "median_transcript_proportions": pl.Float32,
            "sd_transcript_proportions": pl.Float32,
        },
    )

    snakemake = SimpleNamespace(
        input={
            "vep_pq": str(tmp_path / "veff.parquet"),
            "isoform_proportions_pq": str(tmp_path / "isoform_proportions.parquet"),
            "gtf_transcripts": str(tmp_path / "gtf_transcripts.parquet"),
        },
        output={"veff_pq": str(tmp_path / "tissue_specific_vep.parquet")},
        params={},
        config={},
    )
    vep_df.write_parquet(snakemake.input["vep_pq"])
    gtf_transcripts_df.write_parquet(snakemake.input["gtf_transcripts"])
    isoform_proportions_df.write_parquet(snakemake.input["isoform_proportions_pq"])

    runpy.run_path(str(SCRIPT), init_globals={"snakemake": snakemake})

    return pl.read_parquet(snakemake.output["veff_pq"]).unnest("features").to_dicts()


def features_row(chrom: str, gene: str) -> dict:
    """The output row of a variant at `chrom`:100 in `gene`, with a stop gained in the transcript
    with median 0.75 and a missense variant in the transcript with median 0.25"""
    return {
        "chrom": chrom,
        "start": 100,
        "end": 101,
        "ref": "C",
        "alt": "T",
        "gene": gene,
        "tissue": TISSUE,
        "stop_gained.max": True,
        "missense_variant.max": True,
        "stop_gained.sum": 1,
        "missense_variant.sum": 1,
        "stop_gained.proportion": 0.75,
        "missense_variant.proportion": 0.25,
        "num_transcripts": 2,
    }


def test_non_par_transcripts(tmp_path):
    consequences = [
        Consequence("ENST00000000003", stop_gained=True, missense_variant=False),
        Consequence("ENST00000000004", stop_gained=False, missense_variant=True),
    ]

    assert run_tissue_specific_vep(tmp_path, "chr22", NON_PAR_GENE, consequences) == [
        features_row("chr22", NON_PAR_GENE)
    ]


# VEP and mehari give PAR variants on chrX and on chrY the IDs of the chrX copy. Each transcript
# counts once, as in a gene outside the PAR.
@pytest.mark.parametrize("chrom", ["chrX", "chrY"])
def test_par_transcripts_count_once(tmp_path, chrom):
    consequences = [
        Consequence("ENST00000000001", stop_gained=True, missense_variant=False),
        Consequence("ENST00000000002", stop_gained=False, missense_variant=True),
    ]

    assert run_tissue_specific_vep(tmp_path, chrom, PAR_GENE, consequences) == [features_row(chrom, PAR_GENE)]
