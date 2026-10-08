"""AbSplice2-DNA per variant, gene and GTEx tissue, from Pangolin, MMSplice with SpliceMaps, and the SpliceMaps.

Our own code in polars for the steps of the AbSplice2 example workflow after MMSplice and Pangolin
(https://github.com/gagneurlab/absplice2, commit a30120f: `pangolin_postprocess.py`, `pangolin_splicemap.py` and
`absplice_dna.py` in `example/workflow/splicing_pred/DNA`). The AbSplice2 model is not part of this package; the
caller loads it and passes its predictions as a function. Needs the extra `absplice2`.
"""
from abexp.absplice2.dna import FEATURES, GAIN_SCORE_CAP, OUTPUT_COLUMNS, SLACK, absplice2_dna, read_mmsplice_splicemap

__all__ = [
    'FEATURES',
    'GAIN_SCORE_CAP',
    'OUTPUT_COLUMNS',
    'SLACK',
    'absplice2_dna',
    'read_mmsplice_splicemap',
]
