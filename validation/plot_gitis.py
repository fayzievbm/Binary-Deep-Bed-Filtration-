"""Redraw the conditional Gitis comparison from existing CSV records only.

No solver, fitting routine, source-image curve or experimental input is modified.
The vertical label deliberately identifies the published normalized response;
its run-specific measurement interpretation remains conditional in the article.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'gitis'
CURVES = DATA / 'figure6_plot_curves.csv'
POINTS = DATA / 'figure6_experimental_points.csv'
DEPTHS = (0.82, 1.02, 1.62)


def read_csv(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    hashes_before = {p.name: sha256(p) for p in (CURVES, POINTS)}
    curves = read_csv(CURVES)
    points = read_csv(POINTS)
    if len(points) != 27 or len(curves) != 900:
        raise ValueError('Expected the fixed 27-symbol and 900-curve-record inputs.')

    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 12,
        'axes.labelsize': 12, 'axes.titlesize': 11,
        'xtick.labelsize': 11, 'ytick.labelsize': 11,
        'legend.fontsize': 10, 'lines.linewidth': 1.9,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': 0.85, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'savefig.dpi': 300, 'savefig.facecolor': 'white',
    })
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.6), layout='constrained')
    for i, (ax, depth) in enumerate(zip(axes, DEPTHS)):
        c = [r for r in curves if abs(float(r['depth_m']) - depth) < 1e-8]
        d = [r for r in points if abs(float(r['depth_m']) - depth) < 1e-8]
        if len(c) != 300 or len(d) != 9:
            raise ValueError(f'Unexpected data count at depth {depth}.')
        c.sort(key=lambda r: float(r['published_time_h']))
        ax.plot([float(r['published_time_h']) for r in c],
                [float(r['linear_ratio']) for r in c],
                color='#2267A0', linestyle='--', zorder=2)
        ax.plot([float(r['published_time_h']) for r in c],
                [float(r['multistage_ratio']) for r in c],
                color='#D16A24', linestyle='-', zorder=3)
        ax.errorbar([float(r['time_h']) for r in d],
                    [float(r['concentration_ratio']) for r in d],
                    xerr=[float(r['time_digitization_envelope_h']) for r in d],
                    yerr=[float(r['ratio_digitization_envelope']) for r in d],
                    fmt='o', color='black', ecolor='black', ms=4, capsize=2,
                    elinewidth=0.8, capthick=0.8, zorder=5)
        ax.set_title(f'{depth:.2f} m\n' + ('Calibration' if i == 0 else 'Spatial holdout'))
        ax.set_xlabel('Published time (h)')
        ax.set_ylabel('Normalized response')
        ax.set_ylim(0, 0.22)
        ax.grid(alpha=0.15)
        ax.set_axisbelow(True)

    handles = [
        Line2D([0], [0], color='#2267A0', linestyle='--', label='Linear'),
        Line2D([0], [0], color='#D16A24', linestyle='-', label='Multistage'),
        Line2D([0], [0], color='black', marker='o', linestyle='none',
               markersize=4, label='Digitized symbols'),
    ]
    fig.legend(handles=handles, loc='outside lower center', ncol=3, fontsize=10)
    for extension in ('png', 'pdf'):
        fig.savefig(DATA / f'gitis_figure6_conditional_fit.{extension}')
    plt.close(fig)

    hashes_after = {p.name: sha256(p) for p in (CURVES, POINTS)}
    if hashes_before != hashes_after:
        raise RuntimeError('Input CSV changed during rendering.')
    print(json.dumps({'input_sha256_unchanged': hashes_after,
                      'panels': 3, 'observations': 27,
                      'solver_or_fit_invoked': False}, indent=2))


if __name__ == '__main__':
    main()
