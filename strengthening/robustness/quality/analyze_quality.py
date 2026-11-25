"""Continuous effluent exceedance from bundled, previously computed CSVs.

No PDE is solved. NumPy and Matplotlib are the only dependencies.
Run from any directory: python /path/to/quality/analyze_quality.py
"""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / 'inputs'
FEED = np.array([0.05, 0.01])
GRIDS = (80, 160, 320)
MODES = ('pressure', 'rate')
KINETICS = ('capped_linear', 'charging_linear', 'multistage')
THETA = np.array([0.8, 0.9, 0.95])
CURVE = np.linspace(0.5, 1.0, 101)


def metrics(pv, ratios, theta):
    """Integrate max of two linear interpolants, clipped at theta, exactly.

    Exactness refers only to the declared piecewise-linear interpolant of the
    stored outlet samples. Breakpoints include both threshold crossings and
    intersections of the two component ratio lines.
    """
    total = 0.0
    good = 0.0
    for a, b, ya, yb in zip(pv[:-1], pv[1:], ratios[:-1], ratios[1:]):
        cut = [0., 1.]
        d = yb - ya
        for i in range(2):
            if d[i] != 0:
                f = (theta - ya[i]) / d[i]
                if 0 < f < 1: cut.append(float(f))
        dd = d[0] - d[1]
        if dd != 0:
            f = (ya[1] - ya[0]) / dd
            if 0 < f < 1: cut.append(float(f))
        cut.sort()
        for lo, hi in zip(cut[:-1], cut[1:]):
            mid = (lo + hi) / 2
            rmid = ya + mid * d
            if max(rmid) <= theta:
                good += (b - a) * (hi - lo)
            else:
                i = int(np.argmax(rmid))
                vlo = ya[i] + lo * d[i] - theta
                vhi = ya[i] + hi * d[i] - theta
                total += (b - a) * (hi - lo) * (vlo + vhi) / 2
    return total, good


def absolute_integral(x, values):
    """Exact integral of the absolute value of a scalar linear interpolant."""
    out = 0.
    for dx, a, b in zip(np.diff(x), values[:-1], values[1:]):
        if a * b >= 0:
            out += dx * (abs(a) + abs(b)) / 2
        else:
            out += dx * (a*a + b*b) / (2 * (abs(a) + abs(b)))
    return float(out)


def filename(mode, kin, n):
    return f'control_{mode}_{kin}{"" if n == 160 else f"_n{n}"}.csv'


def load(mode, kin, n):
    path = INPUTS / filename(mode, kin, n)
    data = np.genfromtxt(path, delimiter=',', names=True)
    x = data['PV']
    ratios = np.column_stack([data['c1_out'], data['c2_out']]) / FEED
    # Endpoints are inserted by interpolation, making the integration domain
    # explicit even if a future bundled input uses different sample spacing.
    use = (x > 1 + 1e-10) & (x < 14 - 1e-10)
    xx = np.r_[1., x[use], 14.]
    rr = np.column_stack([np.interp(xx, x, ratios[:, i]) for i in range(2)])
    assert np.all(np.diff(xx) > 0) and x[0] <= 1 and x[-1] >= 14 - 1e-9
    return xx, rr, path


def main():
    # Analytical triangular areas exercise threshold and component crossings.
    assert abs(metrics(np.array([0.,1.]),np.array([[0.,0.],[1.,0.]]),.5)[0]-.125)<1e-14
    assert abs(metrics(np.array([0.,1.]),np.array([[1.,0.],[0.,1.]]),.6)[0]-.16)<1e-14
    rows = []
    cases = {}
    curves = {}
    trajectories = {}
    provenance = []
    for mode in MODES:
        for kin in KINETICS:
            for n in GRIDS:
                x, r, path = load(mode, kin, n)
                key = f'{mode}_{kin}_n{n}'
                trajectories[(mode, kin, n)] = (x, r)
                curve = np.array([metrics(x, r, t)[0] for t in CURVE])
                curves[(mode, kin, n)] = curve
                result = []
                for theta in THETA:
                    value, good = metrics(x, r, theta)
                    trap = float(np.trapezoid(np.maximum(r.max(axis=1) - theta, 0), x))
                    coarse = []
                    for stride in (2, 4):
                        ix = np.unique(np.r_[np.arange(0, len(x), stride), len(x)-1])
                        coarse.append(metrics(x[ix], r[ix], theta)[0])
                    item = dict(theta=float(theta), exceedance_PV=value,
                                sampled_interpolant_good_PV=good,
                                nodal_trapezoid_exceedance_PV=trap,
                                trapezoid_minus_exact_interpolant=trap-value,
                                decimated_004_difference_PV=coarse[0]-value,
                                decimated_008_difference_PV=coarse[1]-value)
                    result.append(item)
                    rows.append(dict(mode=mode, kinetics=kin, n=n, **item))
                cases[key] = dict(mode=mode, kinetics=kin, n=n,
                                  stored_PV_spacing_max=float(np.diff(x).max()),
                                  source_file=path.name, metrics=result)
                provenance.append(dict(file='inputs/' + path.name,
                                       sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    checks = []
    for mode in MODES:
        for kin in KINETICS:
            for n0, n1 in ((80, 160), (160, 320)):
                x0, r0 = trajectories[(mode, kin, n0)]
                x1, r1 = trajectories[(mode, kin, n1)]
                x = np.unique(np.r_[x0, x1])
                d = np.column_stack([np.interp(x,x1,r1[:,i])-np.interp(x,x0,r0[:,i]) for i in range(2)])
                # |max(r_f)-max(r_c)| <= sum_i |r_fi-r_ci|. This is
                # a bound between saved interpolants, not a PDE error bound.
                bound = sum(absolute_integral(x, d[:,i]) for i in range(2))
                diffs = [abs(metrics(x1,r1,t)[0]-metrics(x0,r0,t)[0]) for t in CURVE]
                assert max(diffs) <= bound + 1e-12
                checks.append(dict(mode=mode, kinetics=kin, grids=[n0,n1],
                                   maximum_exceedance_difference_over_threshold_curve=max(diffs),
                                   saved_interpolant_L1_difference_bound=bound))
    report = dict(
        scope='Postprocessing of existing capmatched PDE simulations; no new physical equations or simulations.',
        definition='X_theta = integral_1^14 max(max_i(c_i,out/c_i,in)-theta,0) dPV',
        feed=FEED.tolist(), domain_PV=[1,14], illustrative_thresholds=THETA.tolist(),
        lower_is_better='Lower X_theta means less integrated normalized exceedance above this illustrative threshold over the declared fixed horizon.',
        interpolation='Each component outlet ratio is piecewise linear between PV=.02 saved samples. Crossings are inserted analytically, so integration is exact for that interpolant.',
        accuracy_limits=[
            'The saved sampling does not determine an error bound relative to the unsaved continuous solver trajectory or the continuum PDE solution.',
            'Nodal trapezoid differences and decimation to .04/.08 PV quantify sampling sensitivity only; decimation is not an independent PDE convergence study.',
            'The L1 bounds apply between saved piecewise-linear interpolants; they demonstrate continuity of this metric, not a certified PDE error.',
            'The sampled_interpolant_good_PV field integrates the saved interpolant and need not equal the good-PV counter accumulated at internal solver steps.',
            'All cases use PV14 and omit the first injected pore volume. This is neither an optimized stopping policy nor a water-quality standard.'
        ], cases=cases, grid_comparisons=checks, source_provenance=provenance)
    (ROOT/'quality_results.json').write_text(json.dumps(report, indent=2)+'\n')
    with (ROOT/'quality_metrics.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    with (ROOT/'quality_curves.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['mode','kinetics','n','theta','exceedance_PV'])
        for (mode,kin,n),v in curves.items():
            writer.writerows([mode,kin,n,float(t),float(z)] for t,z in zip(CURVE,v))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':300})
    fig,axes=plt.subplots(1,2,figsize=(9,3.8),sharey=True,layout='constrained')
    colors={'capped_linear':'#2267A0','charging_linear':'#22885B','multistage':'#D16A24'}
    labels={'capped_linear':'Capped linear reference','charging_linear':'Charging + linear','multistage':'Multistage'}
    for ax,mode,label in zip(axes,MODES,['Fixed pressure','Fixed flow']):
        ax.axhline(0,color='#666666',lw=1,label=labels['capped_linear'])
        for kin in ('charging_linear','multistage'):
            stack=np.array([curves[(mode,kin,n)]-curves[(mode,'capped_linear',n)] for n in GRIDS])
            ax.fill_between(CURVE,stack.min(axis=0),stack.max(axis=0),color=colors[kin],alpha=.17,lw=0)
            ax.plot(CURVE,stack[-1],color=colors[kin],lw=1.8,label=labels[kin])
        ax.set(xlabel='Illustrative threshold $\\theta$',title=label,xlim=(.5,1))
        ax.grid(alpha=.18)
    axes[0].set_ylabel('Exceedance difference $\\Delta X_\\theta$ (PV)')
    axes[1].legend(fontsize=9,loc='upper left')
    fig.get_layout_engine().set(rect=(0,.07,1,.93))
    fig.text(.5,.015,'Positive: more exceedance. Lines: 320 cells; bands: 80–320 cells. Domain: PV 1–14.',ha='center',fontsize=9)
    for ext in ('png','pdf'):fig.savefig(ROOT/f'figure_continuous_quality.{ext}')
    plt.close(fig)
    print(json.dumps({'cases':len(cases),'metrics':len(rows),
                      'maximum_nodal_trapezoid_error':max(abs(v['trapezoid_minus_exact_interpolant']) for v in rows),
                      'maximum_decimation_004_difference':max(abs(v['decimated_004_difference_PV']) for v in rows),
                      'maximum_decimation_008_difference':max(abs(v['decimated_008_difference_PV']) for v in rows)},indent=2))


if __name__ == '__main__':main()
