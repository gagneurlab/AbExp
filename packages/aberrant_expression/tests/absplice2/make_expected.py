"""Print the expected rows of test_absplice2.py, computed by the upstream AbSplice2 scripts on constellations.py.

Run this script only to change the expected rows:

    python make_expected.py <checkout of gagneurlab/absplice2 at a30120f>

It runs the three scripts of `example/workflow/splicing_pred/DNA` unchanged, with a stand-in `snakemake` object:
`pangolin_postprocess.py`, `pangolin_splicemap.py` and `absplice_dna.py`. They need the environment of the AbSplice2
example workflow: Python 3.8, pandas, numpy, pysam, pyranges 0.1.4, splicemap 0.0.2 and tqdm. It does not need
aberrant-expression or interpret.

The inputs of each test go into a temporary folder: the SpliceMaps and the MMSplice table from write_inputs, and
the VCF file that Pangolin would write for the rows of abexp.pangolin. Pangolin rounds the scores to 2 decimals in
float32 and keeps the version of the gene id. `absplice_dna.py` reads the model from a path relative to the working
directory. There it finds a pickle of `StubModel`, the stand-in of constellations.py.

The rows are printed in the order and form of abexp.absplice2.absplice2_dna: the variant split into chrom, start
(0-based), end, ref and alt, the tissue renamed with TISSUE_MAPPING, and the rows sorted by these key columns.

Of tied rows, the upstream scripts keep any one. For `tied_rows`, the script therefore runs them twice, with the
MMSplice rows in the order of the constellation and reversed, and prints both results.
"""
import argparse
import dataclasses
import math
import os
import pickle
import runpy
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import constellations as c  # noqa: E402

UPSTREAM_COMMIT = 'a30120f'
OUTPUT_COLUMNS = ['chrom', 'start', 'end', 'ref', 'alt', 'gene', 'tissue', 'gain_score', 'gain_pos', 'loss_score',
                  'loss_pos', 'pangolin_tissue_score', 'ref_psi_pangolin', 'median_n_pangolin', 'ref_psi', 'median_n',
                  'delta_logit_psi', 'delta_psi', 'junction', 'event_type', 'splice_site', 'AbSplice_DNA']


class StubModel:
    """The stand-in for the AbSplice2 model, see constellations.py."""

    def predict_proba(self, features):
        x = np.asarray(features, dtype=np.float64)
        z = c.BIAS + x @ np.asarray(c.WEIGHTS)
        p = 1 / (1 + np.exp(-z))
        return np.column_stack([1 - p, p])


def pangolin_text(score):
    """A score as Pangolin writes it: rounded to 2 decimals in float32."""
    return str(round(np.float32(score), 2))


def write_pangolin_vcf(path, rows):
    """The VCF file of Pangolin for the rows of abexp.pangolin: one record per variant, one entry per gene."""
    entries = {}
    for r in rows:
        entry = (f'{r.gene_id}|{r.gain_pos}:{pangolin_text(r.gain_score)}|'
                 f'{r.loss_pos}:{pangolin_text(r.loss_score)}|Warnings:')
        entries.setdefault(r.variant, []).append(entry)
    chroms = sorted({variant.split(':')[0] for variant in entries})
    with open(path, 'w') as fd:
        fd.write('##fileformat=VCFv4.2\n')
        for chrom in chroms:
            fd.write(f'##contig=<ID={chrom}>\n')
        fd.write('##INFO=<ID=Pangolin,Number=.,Type=String,Description="Pangolin splice scores. '
                 'Format: gene|pos:score_change|pos:score_change|warnings,...">\n')
        fd.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        for variant, gene_entries in entries.items():
            parts = variant.split(':')
            alleles = parts[2].split('>')
            fd.write(f'{parts[0]}\t{parts[1]}\t.\t{alleles[0]}\t{alleles[1]}\t.\t.\tPangolin={",".join(gene_entries)}\n')


def run_upstream(scripts_dir, constellation, tmp):
    """The output table of `absplice_dna.py` for `constellation`."""
    inputs = c.write_inputs(constellation, tmp)
    write_pangolin_vcf(tmp / 'pangolin.vcf', constellation.pangolin)
    # absplice_dna.py reads ../../absplice/precomputed/AbSplice_2_DNA.pkl
    workdir = tmp / 'example' / 'workflow'
    workdir.mkdir(parents=True)
    (tmp / 'absplice' / 'precomputed').mkdir(parents=True)
    with open(tmp / 'absplice' / 'precomputed' / 'AbSplice_2_DNA.pkl', 'wb') as fd:
        pickle.dump(StubModel(), fd)
    steps = [
        ('pangolin_postprocess.py', {'pangolin_raw': str(tmp / 'pangolin.vcf')},
         {'pangolin_csv': str(tmp / 'pangolin.csv')}),
        ('pangolin_splicemap.py', {'pangolin_csv': str(tmp / 'pangolin.csv'), 'splicemap_5': inputs.splicemap5,
                                   'splicemap_3': inputs.splicemap3},
         {'pangolin_splicemap': str(tmp / 'pangolin_splicemap.csv')}),
        ('absplice_dna.py', {'pangolin_splicemap': str(tmp / 'pangolin_splicemap.csv'),
                             'mmsplice_splicemap': inputs.mmsplice},
         {'absplice_dna': str(tmp / 'absplice_dna.csv')}),
    ]
    cwd = Path.cwd()
    try:
        os.chdir(workdir)
        for script, step_inputs, step_outputs in steps:
            snakemake = SimpleNamespace(input=step_inputs, output=step_outputs)
            runpy.run_path(str(scripts_dir / script), init_globals={'snakemake': snakemake})
    finally:
        os.chdir(cwd)
    return pd.read_csv(tmp / 'absplice_dna.csv')


def value(x):
    """A cell of the upstream table as a Python value, with None for NaN."""
    if isinstance(x, float) and math.isnan(x):
        return None
    if isinstance(x, (np.integer, np.floating)):
        return x.item()
    return x


def expected_rows(df):
    """The rows of the upstream table in the order and form of absplice2_dna."""
    rows = []
    for record in df.to_dict('records'):
        parts = record['variant'].split(':')
        alleles = parts[2].split('>')
        start = int(parts[1]) - 1
        row = {
            'chrom': parts[0],
            'start': start,
            'end': start + len(alleles[0]),
            'ref': alleles[0],
            'alt': alleles[1],
            'gene': record['gene_id'],
            'tissue': c.TISSUE_MAPPING.get(record['tissue'], record['tissue']),
        }
        for column in OUTPUT_COLUMNS[7:]:
            row[column] = value(record[column])
        for column in ['gain_pos', 'loss_pos']:
            if row[column] is not None:
                row[column] = int(row[column])
        rows.append(tuple(row[column] for column in OUTPUT_COLUMNS))
    return sorted(rows, key=lambda row: row[:7])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('absplice2', type=Path, help='checkout of gagneurlab/absplice2')
    args = parser.parse_args()
    commit = subprocess.run(['git', '-C', str(args.absplice2), 'rev-parse', 'HEAD'], check=True,
                            capture_output=True, text=True).stdout.strip()
    assert commit.startswith(UPSTREAM_COMMIT), f'{args.absplice2} is at {commit}, not at {UPSTREAM_COMMIT}'
    scripts_dir = args.absplice2.resolve() / 'example' / 'workflow' / 'splicing_pred' / 'DNA'

    cases = dict(c.CONSTELLATIONS)
    cases['tied_rows (MMSplice rows reversed)'] = dataclasses.replace(
        c.TIED_ROWS, mmsplice=tuple(reversed(c.TIED_ROWS.mmsplice)))
    for name, constellation in cases.items():
        with tempfile.TemporaryDirectory() as tmp:
            df = run_upstream(scripts_dir, constellation, Path(tmp))
        print(name)
        for row in expected_rows(df):
            print(f'    {row},')


if __name__ == '__main__':
    main()
