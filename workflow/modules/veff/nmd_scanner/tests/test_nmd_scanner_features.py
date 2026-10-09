"""Tests of the rule veff__nmd_scanner_features

Each test runs `nmd_scanner_features.py.py` on small input tables of one variant in one gene, and
checks the full output row. Run them in the environment of
`workflow/modules/veff/envs/veff-py.yaml`, with pytest added:

    pytest workflow/modules/veff/nmd_scanner/tests
"""
import runpy
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import polars as pl
import pytest

SCRIPT = Path(__file__).parent.parent / "nmd_scanner_features.py.py"

VARIANT = {"chrom": "chr22", "start": 100, "end": 101, "ref": "C", "alt": "T"}
GENE = "ENSG00000000001"
TISSUE = "Whole Blood"

# the flags of the score table
FLAGS = [
    "start_loss",
    "stop_loss",
    "nmd_last_exon_rule",
    "nmd_50nt_penultimate_rule",
    "nmd_long_exon_rule",
    "nmd_start_proximal_rule",
    "nmd_single_exon_rule",
    "ptc_less_than_150nt_to_start",
]


@dataclass(frozen=True)
class Transcript:
    """A transcript of GENE: its row of NMD-Scanner's table for VARIANT, and its score if
    `nmd_model_status` is "ok"."""

    transcript: str
    nmd_model_status: str
    alt_has_ptc: bool | None
    nmd_pred_score: float | None = None
    nmd_escape: bool = False


@dataclass(frozen=True)
class IsoformProportion:
    """A row of the isoform proportion table"""

    transcript: str
    median_transcript_proportions: float | None
    gene: str = GENE
    tissue: str = TISSUE


def run_features(tmp_path: Path, transcripts: list[Transcript], proportions: list[IsoformProportion]) -> list[dict]:
    """Writes the input tables, runs the script and returns its rows, with the struct `features`
    unnested."""
    keys = {**VARIANT, "gene": GENE}
    key_schema = {
        "chrom": pl.String,
        "start": pl.Int64,
        "end": pl.Int64,
        "ref": pl.String,
        "alt": pl.String,
        "gene": pl.String,
        "transcript": pl.String,
    }

    nmd_scanner_df = pl.DataFrame(
        [
            {**keys, "transcript": t.transcript, "alt_has_ptc": t.alt_has_ptc, "nmd_model_status": t.nmd_model_status}
            for t in transcripts
        ],
        schema={**key_schema, "alt_has_ptc": pl.Boolean, "nmd_model_status": pl.Categorical},
    )
    score_df = pl.DataFrame(
        [
            {
                **keys,
                "transcript": t.transcript,
                "nmd_pred_score": t.nmd_pred_score,
                "alt_has_ptc": 1.0,
                "nmd_escape": float(t.nmd_escape),
                **{c: 0.0 for c in FLAGS},
            }
            for t in transcripts
            if t.nmd_model_status == "ok"
        ],
        schema={**key_schema, **{c: pl.Float64 for c in ["nmd_pred_score", "alt_has_ptc", "nmd_escape", *FLAGS]}},
    )
    isoform_proportions_df = pl.DataFrame(
        [
            {
                "gene": p.gene,
                "tissue": p.tissue,
                "transcript": p.transcript,
                "median_transcript_proportions": p.median_transcript_proportions,
            }
            for p in proportions
        ],
        schema={
            "gene": pl.String,
            "tissue": pl.String,
            "transcript": pl.String,
            "median_transcript_proportions": pl.Float32,
        },
    )

    snakemake = SimpleNamespace(
        input={
            "nmd_scanner_pq": str(tmp_path / "veff.parquet"),
            "nmd_score_pq": str(tmp_path / "score.parquet"),
            "isoform_proportions_pq": str(tmp_path / "isoform_proportions.parquet"),
        },
        output={"veff_pq": str(tmp_path / "features.parquet")},
    )
    nmd_scanner_df.write_parquet(snakemake.input["nmd_scanner_pq"])
    score_df.write_parquet(snakemake.input["nmd_score_pq"])
    isoform_proportions_df.write_parquet(snakemake.input["isoform_proportions_pq"])

    runpy.run_path(str(SCRIPT), init_globals={"snakemake": snakemake})

    return pl.read_parquet(snakemake.output["veff_pq"]).unnest("features").to_dicts()


def features_row(tissue: str = TISSUE, **features) -> dict:
    """The output row of VARIANT, GENE and `tissue`"""
    return {**VARIANT, "gene": GENE, "tissue": tissue, **features}


def flag_features(proportion: float | None, maximum: float | None) -> dict:
    """`<flag>.proportion` and `<flag>` of each flag in FLAGS"""
    return {**{f"{c}.proportion": proportion for c in FLAGS}, **{c: maximum for c in FLAGS}}


# The medians of the gene add up to 1 in the tissue, so they are the weights.
PROPORTIONS = [
    IsoformProportion("ENST00000000001", 0.5),
    IsoformProportion("ENST00000000002", 0.25),
    IsoformProportion("ENST00000000003", 0.25),
]


def test_ref_ptc_transcript_is_no_ptc(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=2.0),
        # the reference has the PTC already
        Transcript("ENST00000000002", "ref_ptc", True),
    ]

    assert run_features(tmp_path, transcripts, PROPORTIONS) == [
        features_row(
            **{
                "alt_has_ptc.proportion": 0.5,
                "num_ptc": 1,
                "nmd_pred_score.weighted_sum": 1.0,
                "nmd_escape.proportion": 0.0,
                "nmd_pred_score.weighted_mean": 2.0,
                "nmd_pred_score": 2.0,
                "nmd_pred_score.median": 2.0,
                "nmd_pred_score.mean": 2.0,
                "nmd_pred_score.std": None,
                "nmd_pred_score.high_proportion_weighted_max": 0.0,
                "num_escape": 0,
                "nmd_pred_score.high_expr_max": 2.0,
            },
            **flag_features(0.0, 0.0),
        )
    ]


def test_null_alt_has_ptc_is_no_ptc(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=2.0),
        Transcript("ENST00000000002", "unknown_effect", None),
    ]

    assert run_features(tmp_path, transcripts, PROPORTIONS) == [
        features_row(
            **{
                "alt_has_ptc.proportion": 0.5,
                "num_ptc": 1,
                "nmd_pred_score.weighted_sum": 1.0,
                "nmd_escape.proportion": 0.0,
                "nmd_pred_score.weighted_mean": 2.0,
                "nmd_pred_score": 2.0,
                "nmd_pred_score.median": 2.0,
                "nmd_pred_score.mean": 2.0,
                "nmd_pred_score.std": None,
                "nmd_pred_score.high_proportion_weighted_max": 0.0,
                "num_escape": 0,
                "nmd_pred_score.high_expr_max": 2.0,
            },
            **flag_features(0.0, 0.0),
        )
    ]


def test_unscored_ptc_transcript_is_ptc(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=2.0),
        Transcript("ENST00000000002", "no_annotated_stop", True),
    ]

    # The score features come from ENST00000000001 alone.
    assert run_features(tmp_path, transcripts, PROPORTIONS) == [
        features_row(
            **{
                "alt_has_ptc.proportion": 0.75,
                "num_ptc": 2,
                "nmd_pred_score.weighted_sum": 1.0,
                "nmd_escape.proportion": 0.0,
                "nmd_pred_score.weighted_mean": 2.0,
                "nmd_pred_score": 2.0,
                "nmd_pred_score.median": 2.0,
                "nmd_pred_score.mean": 2.0,
                "nmd_pred_score.std": None,
                "nmd_pred_score.high_proportion_weighted_max": 0.0,
                "num_escape": 0,
                "nmd_pred_score.high_expr_max": 2.0,
            },
            **flag_features(0.0, 0.0),
        )
    ]


def test_no_scored_ptc_transcript(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "no_ptc", False),
        Transcript("ENST00000000002", "no_annotated_start", True),
    ]

    assert run_features(tmp_path, transcripts, PROPORTIONS) == [
        features_row(
            **{
                "alt_has_ptc.proportion": 0.25,
                "num_ptc": 1,
                "nmd_pred_score.weighted_sum": None,
                "nmd_escape.proportion": None,
                "nmd_pred_score.weighted_mean": None,
                "nmd_pred_score": None,
                "nmd_pred_score.median": None,
                "nmd_pred_score.mean": None,
                "nmd_pred_score.std": None,
                "nmd_pred_score.high_proportion_weighted_max": None,
                "num_escape": None,
                "nmd_pred_score.high_expr_max": None,
            },
            **flag_features(None, None),
        )
    ]


@pytest.mark.parametrize("median", [0.0, None])
def test_gene_without_weights(tmp_path, median):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=2.0),
        Transcript("ENST00000000002", "no_annotated_stop", True),
    ]
    # The medians of the gene add up to 0 in the tissue.
    proportions = [
        IsoformProportion("ENST00000000001", median),
        IsoformProportion("ENST00000000002", median),
        IsoformProportion("ENST00000000003", median),
    ]

    assert run_features(tmp_path, transcripts, proportions) == [
        features_row(
            **{
                "alt_has_ptc.proportion": None,
                "num_ptc": 2,
                "nmd_pred_score.weighted_sum": None,
                "nmd_escape.proportion": None,
                "nmd_pred_score.weighted_mean": None,
                "nmd_pred_score": 2.0,
                "nmd_pred_score.median": 2.0,
                "nmd_pred_score.mean": 2.0,
                "nmd_pred_score.std": None,
                "nmd_pred_score.high_proportion_weighted_max": None,
                "num_escape": 0,
                "nmd_pred_score.high_expr_max": None,
            },
            **flag_features(None, 0.0),
        )
    ]


def test_ptc_transcripts_without_weight_in_gene_with_weights(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=2.0),
        Transcript("ENST00000000002", "ok", True, nmd_pred_score=1.0, nmd_escape=True),
    ]
    # The medians add up to 1, but both PTC transcripts have the median 0.
    proportions = [
        IsoformProportion("ENST00000000001", 0.0),
        IsoformProportion("ENST00000000002", 0.0),
        IsoformProportion("ENST00000000003", 1.0),
    ]

    assert run_features(tmp_path, transcripts, proportions) == [
        features_row(
            **{
                "alt_has_ptc.proportion": 0.0,
                "num_ptc": 2,
                "nmd_pred_score.weighted_sum": 0.0,
                "nmd_escape.proportion": 0.0,
                "nmd_pred_score.weighted_mean": None,
                "nmd_pred_score": 2.0,
                "nmd_pred_score.median": 1.5,
                "nmd_pred_score.mean": 1.5,
                "nmd_pred_score.std": 0.5**0.5,
                "nmd_pred_score.high_proportion_weighted_max": 0.0,
                "num_escape": 1,
                "nmd_pred_score.high_expr_max": None,
            },
            **flag_features(0.0, 0.0),
        )
    ]


def test_medians_that_sum_to_more_than_1(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=1.0),
        Transcript("ENST00000000002", "ok", True, nmd_pred_score=2.0, nmd_escape=True),
    ]
    # The medians add up to 2, so the weights are 0.5, 0.125 and 0.375. ENST00000000003 has no PTC,
    # but counts in the sum.
    proportions = [
        IsoformProportion("ENST00000000001", 1.0),
        IsoformProportion("ENST00000000002", 0.25),
        IsoformProportion("ENST00000000003", 0.75),
    ]

    assert run_features(tmp_path, transcripts, proportions) == [
        features_row(
            **{
                "alt_has_ptc.proportion": 0.625,
                "num_ptc": 2,
                "nmd_pred_score.weighted_sum": 0.75,
                "nmd_escape.proportion": 0.125,
                "nmd_pred_score.weighted_mean": 0.75 / 0.625,
                "nmd_pred_score": 2.0,
                "nmd_pred_score.median": 1.5,
                "nmd_pred_score.mean": 1.5,
                "nmd_pred_score.std": 0.5**0.5,
                # no weight above 0.8, though the median of ENST00000000001 is 1
                "nmd_pred_score.high_proportion_weighted_max": 0.0,
                "num_escape": 1,
                # ENST00000000002 has a weight below 0.2, though its median is 0.25
                "nmd_pred_score.high_expr_max": 1.0,
            },
            **flag_features(0.0, 0.0),
        )
    ]


def test_transcript_of_another_gene_in_isoform_table(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=2.0),
        Transcript("ENST00000000002", "ok", True, nmd_pred_score=1.0),
    ]
    # The table puts ENST00000000002 into another gene, so it has no weight in GENE. It counts in
    # the unweighted features.
    proportions = [
        IsoformProportion("ENST00000000001", 0.5),
        IsoformProportion("ENST00000000003", 0.5),
        IsoformProportion("ENST00000000002", 1.0, gene="ENSG00000000002"),
    ]

    assert run_features(tmp_path, transcripts, proportions) == [
        features_row(
            **{
                "alt_has_ptc.proportion": 0.5,
                "num_ptc": 2,
                "nmd_pred_score.weighted_sum": 1.0,
                "nmd_escape.proportion": 0.0,
                "nmd_pred_score.weighted_mean": 2.0,
                "nmd_pred_score": 2.0,
                "nmd_pred_score.median": 1.5,
                "nmd_pred_score.mean": 1.5,
                "nmd_pred_score.std": 0.5**0.5,
                "nmd_pred_score.high_proportion_weighted_max": 0.0,
                "num_escape": 0,
                "nmd_pred_score.high_expr_max": 2.0,
            },
            **flag_features(0.0, 0.0),
        )
    ]


def test_ptc_transcript_missing_from_isoform_table(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=2.0),
        # e.g. a transcript of a later GENCODE release than GTEx
        Transcript("ENST00000000004", "ok", True, nmd_pred_score=3.0, nmd_escape=True),
    ]

    # ENST00000000004 has no weight. It counts in the unweighted features.
    assert run_features(tmp_path, transcripts, PROPORTIONS) == [
        features_row(
            **{
                "alt_has_ptc.proportion": 0.5,
                "num_ptc": 2,
                "nmd_pred_score.weighted_sum": 1.0,
                "nmd_escape.proportion": 0.0,
                "nmd_pred_score.weighted_mean": 2.0,
                "nmd_pred_score": 3.0,
                "nmd_pred_score.median": 2.5,
                "nmd_pred_score.mean": 2.5,
                "nmd_pred_score.std": 0.5**0.5,
                "nmd_pred_score.high_proportion_weighted_max": 0.0,
                "num_escape": 1,
                "nmd_pred_score.high_expr_max": 2.0,
            },
            **flag_features(0.0, 0.0),
        )
    ]


def test_gene_missing_from_isoform_table(tmp_path):
    transcripts = [
        Transcript("ENST00000000001", "ok", True, nmd_pred_score=2.0),
        Transcript("ENST00000000002", "no_annotated_stop", True),
    ]
    # The table has another gene in two tissues, and no row of GENE.
    proportions = [
        IsoformProportion("ENST00000000005", 1.0, gene="ENSG00000000002", tissue="Lung"),
        IsoformProportion("ENST00000000005", 1.0, gene="ENSG00000000002", tissue=TISSUE),
    ]
    # GENE gets a row in each tissue of the table, with null weighted features.
    features = {
        "alt_has_ptc.proportion": None,
        "num_ptc": 2,
        "nmd_pred_score.weighted_sum": None,
        "nmd_escape.proportion": None,
        "nmd_pred_score.weighted_mean": None,
        "nmd_pred_score": 2.0,
        "nmd_pred_score.median": 2.0,
        "nmd_pred_score.mean": 2.0,
        "nmd_pred_score.std": None,
        "nmd_pred_score.high_proportion_weighted_max": None,
        "num_escape": 0,
        "nmd_pred_score.high_expr_max": None,
        **flag_features(None, 0.0),
    }

    assert run_features(tmp_path, transcripts, proportions) == [
        features_row("Lung", **features),
        features_row(TISSUE, **features),
    ]
