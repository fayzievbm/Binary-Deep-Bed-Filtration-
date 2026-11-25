"""History-selected interior equilibria in the unchanged filtration framework.

All runs start with clean water and zero retained particles. They differ only
in a finite preconditioning inlet history, after which the inlet and imposed
pressure drop are identical. Numerical outputs are illustrative, not fitted.
"""
from __future__ import annotations

from dataclasses import asdict
import csv
import importlib.util
import json
from pathlib import Path
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'history'
OUT.mkdir(exist_ok=True)
SOLVER_PATH = ROOT / 'numerics' / 'solver.py'
FINAL_FEED = np.array([0.002, 0.0005])
HIGH_FEED = np.array([0.006, 0.0015])
PRE_END = 2.0
END = 10.0
GRIDS = (80, 160)
CASES = {
    'LL': {'high': (False, False), 'label': 'No high preconditioning', 'color': '#2267A0', 'ls': '-'},
    'HL': {'high': (True, False), 'label': 'Type 1 preconditioned', 'color': '#D16A24', 'ls': '--'},
    'LH': {'high': (False, True), 'label': 'Type 2 preconditioned', 'color': '#22885B', 'ls': '-.'},
    'HH': {'high': (True, True), 'label': 'Both types preconditioned', 'color': '#8A56A3', 'ls': '-'},
}


def load_solver():
    spec = importlib.util.spec_from_file_location('novelty_solver_history', SOLVER_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def inlet_for(case):
    initial = np.where(CASES[case]['high'], HIGH_FEED, FINAL_FEED)
    def inlet(t):
        return initial.copy() if t < PRE_END - 1e-12 else FINAL_FEED.copy()
    return inlet


def net_source(B, c, p):
    """Printed piecewise source, used only as a diagnostic away from thresholds."""
    lo = np.asarray(p.B1)[:, None]
    mid = np.asarray(p.B2)[:, None]
    cap = np.asarray(p.cap)[:, None]
    alpha = np.asarray(p.alpha_s)[:, None]
    beta = np.asarray(p.beta_s)[:, None]
    chi = np.asarray(p.charging)[:, None]
    q = np.where(B < lo, chi * alpha * c,
                 np.where(B <= mid, alpha * c - beta * B,
                          alpha * cap * c / np.maximum(B, 1e-300) - beta * B))
    return np.where(B >= cap, 0.0, q)


def exact_state(p, case):
    ratio = np.asarray(p.alpha_s) * FINAL_FEED / np.asarray(p.beta_s)
    high = np.sqrt(np.asarray(p.cap) * ratio)
    B = np.where(CASES[case]['high'], high, ratio)
    occupancy = float(np.sum(B / (1 - np.asarray(p.eps))))
    k = float((1 - np.sqrt(occupancy)) ** 3)
    mu = float(1 + 2.5 * FINAL_FEED.sum())
    u = p.forcing * k / mu
    return {'B': B, 'occupancy': occupancy, 'k_over_k0': k,
            'mu_over_mu0': mu, 'flow': u, 'pressure': p.forcing}


def save_run(case, n, r, p, solver):
    stem = f'history_{case}_n{n}'
    t = np.asarray(r['time'])
    feed_values = np.asarray([inlet_for(case)(v) for v in t])
    h = np.asarray(r['hist'])
    # The existing history columns are retained verbatim; only the declared
    # inlet values at the recorded times are added in this study's CSV.
    table = np.column_stack([h[:, :15], feed_values])
    header = ('t_star,c1_out,c2_out,u_star,delta_p_star,B1_mean,B2_mean,'
              'inventory1,inventory2,input1,input2,output1,output2,'
              'max_balance1,max_balance2,c1_in,c2_in')
    np.savetxt(OUT / (stem + '.csv'), table, delimiter=',', header=header, comments='')
    y = np.asarray(r['final'])
    D, phi, c, k, u, dp, alpha, beta = solver.fields(y, p)
    q = net_source(y[1], c, p)
    expected = exact_state(p, case)
    np.savez_compressed(OUT / (stem + '_fields.npz'),
                        final=y, x_star=(np.arange(n) + 0.5) / n,
                        c=c, q=q, occupancy=D, mobile_fraction=phi, k_over_k0=k)
    qmax = np.max(np.abs(q), axis=1)
    delta = np.max(np.abs(y[1] - expected['B'][:, None]), axis=1)
    branch = np.where(y[1] > np.asarray(p.B2)[:, None], 'aging', 'transition')
    metrics = {
        'case': case, 'n': n, 'label': CASES[case]['label'],
        'preconditioned_species': [i + 1 for i, flag in enumerate(CASES[case]['high']) if flag],
        'preconditioning_inlet': np.where(CASES[case]['high'], HIGH_FEED, FINAL_FEED).tolist(),
        'common_final_inlet': FINAL_FEED.tolist(),
        'preconditioning_end_t_star': PRE_END, 'final_t_star': END,
        'final_mean_retained': np.mean(y[1], axis=1).tolist(),
        'final_min_retained': np.min(y[1], axis=1).tolist(),
        'final_max_retained': np.max(y[1], axis=1).tolist(),
        'exact_interior_equilibrium_retained': expected['B'].tolist(),
        'maximum_deposit_error_to_exact_uniform_equilibrium': delta.tolist(),
        'maximum_absolute_net_deposition_profile': qmax.tolist(),
        'maximum_mobile_concentration_error_to_final_inlet': np.max(abs(c - FINAL_FEED[:, None]), axis=1).tolist(),
        'final_branch_by_species': [sorted(set(row)) for row in branch],
        'final_flow': float(u), 'exact_uniform_equilibrium_flow': expected['flow'],
        'absolute_final_flow_error': abs(float(u) - expected['flow']),
        'final_pressure': float(dp),
        'maximum_inventory_balance_error': np.asarray(r['maxerr']).tolist(),
        'minimum_saturation': float(r['minstate']),
        'minimum_carrier_saturation': float(r['min_carrier']),
        'maximum_occupancy': float(r['maxD']),
        'cap_events': int(r['cap_hits']), 'holding_events': int(r['holds']),
        'steps': int(r['steps']),
        'total_particle_input': np.asarray(r['incoming']).tolist(),
        'total_particle_output': np.asarray(r['outgoing']).tolist(),
    }
    return metrics


def phase_piece(B, i, p):
    c = FINAL_FEED[i]
    return np.where(B < p.B1[i], p.charging[i] * c,
                    np.where(B <= p.B2[i], c - B,
                             p.cap[i] * c / np.maximum(B, 1e-300) - B))


def plot_runs(runs, p):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                         'axes.labelsize': 12, 'axes.titlesize': 12,
                         'xtick.labelsize': 11, 'ytick.labelsize': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'savefig.dpi': 300,
                         'lines.linewidth': 1.8})
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.5), sharex=True)
    fig.subplots_adjust(left=.10, right=.98, top=.93, bottom=.17, hspace=.35, wspace=.30)
    x = [0, PRE_END, PRE_END, END]
    for i, color in enumerate(('#2267A0', '#D16A24')):
        axes[0, 0].plot(x, [HIGH_FEED[i], HIGH_FEED[i], FINAL_FEED[i], FINAL_FEED[i]],
                        color=color, label=f'Type {i + 1}, preconditioned')
        axes[0, 0].plot(x, [FINAL_FEED[i]] * 4, color=color, ls=':',
                        label=f'Type {i + 1}, baseline')
    axes[0, 0].set_ylabel('Inlet volume fraction')
    axes[0, 0].legend(fontsize=8, loc='upper right')
    axes[0, 0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    # Draw the double-high solid curve first so dashed coincident trajectories
    # remain visible instead of hiding the fact that they share a plateau.
    for case in ('HH', 'LL', 'LH', 'HL'):
        r = runs[(160, case)]
        style = CASES[case]
        axes[0, 1].plot(r['time'], r['retained'][:, 0], color=style['color'], ls=style['ls'])
        axes[1, 0].plot(r['time'], r['retained'][:, 1], color=style['color'], ls=style['ls'])
        axes[1, 1].plot(r['time'], r['u'], color=style['color'], ls=style['ls'])
    axes[0, 1].set_ylabel('Mean retained type 1')
    axes[1, 0].set_ylabel('Mean retained type 2')
    axes[1, 1].set_ylabel('Darcy flux $U^*$')
    for ax, title, letter in zip(axes.flat,
        ('Prescribed inlet histories', 'Type 1 loading', 'Type 2 loading', 'Common fixed pressure'),
        ('a', 'b', 'c', 'd')):
        ax.set_title(f'({letter}) {title}', loc='left', pad=10)
        ax.axvspan(0, PRE_END, color='#D8D8D8', alpha=.25, lw=0)
        ax.axvline(PRE_END, color='#777777', lw=.8, ls=':')
        ax.grid(axis='y', alpha=.18)
        ax.set_xlim(0, END)
    for ax in axes[1]:ax.set_xlabel('Time $t^*=U_0t/L$')
    handles = [Line2D([0], [0], color=v['color'], ls=v['ls'], label=v['label']) for v in CASES.values()]
    fig.legend(handles=handles, loc='lower center', ncol=2, fontsize=10, bbox_to_anchor=(.5, .015))
    for ext in ('png', 'pdf'):fig.savefig(OUT / f'figure_history_multistability.{ext}')
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.6), layout='constrained')
    low = FINAL_FEED
    high = np.sqrt(np.asarray(p.cap) * FINAL_FEED)
    for i, ax in enumerate(axes[:2]):
        right = min(p.cap[i] * .9, high[i] * 1.4)
        intervals = [(0, p.B1[i], 'Charging'), (p.B1[i], p.B2[i], 'Transition'), (p.B2[i], right, 'Aging')]
        for a, b, label in intervals:
            bvals = np.linspace(a, b, 160)
            if label == 'Charging':bvals[-1] = np.nextafter(p.B1[i], 0)
            if label == 'Aging':bvals[0] = np.nextafter(p.B2[i], np.inf)
            ax.plot(1e3 * bvals, 1e3 * phase_piece(bvals, i, p), color='#303840')
        ax.axhline(0, color='#777777', lw=.8)
        ax.axvline(1e3 * p.B2[i], color='#999999', ls=':', lw=.9)
        ax.scatter(1e3 * np.array([low[i], high[i]]), [0, 0], color=['#2267A0', '#D16A24'], s=28, zorder=5)
        for value, label in [(low[i], 'L'), (high[i], 'H')]:
            ax.annotate(label, xy=(1e3 * value, 0), xytext=(0, 8),
                        textcoords='offset points', ha='center', fontsize=10)
        ax.set_title(f'({chr(97 + i)}) Type {i + 1} local source', loc='left')
        ax.set_xlabel(f'Deposit $10^3B_{i + 1}$')
        ax.set_ylabel(f'Net source $10^3q_{i + 1}^*$')
        ax.grid(alpha=.15)
    exact_flows = [exact_state(p, case)['flow'] for case in CASES]
    axes[2].bar(range(4), exact_flows, color=[v['color'] for v in CASES.values()], width=.65)
    axes[2].scatter(range(4), [runs[(160, c)]['u'][-1] for c in CASES], color='black', marker='x', s=30,
                    label='Numerical final state')
    axes[2].set_xticks(range(4), CASES)
    axes[2].set_ylim(0, 10)
    axes[2].set_ylabel('Darcy flux $U^*$')
    axes[2].set_title('(c) Four interior combinations', loc='left')
    axes[2].grid(axis='y', alpha=.15)
    axes[2].legend(fontsize=8, loc='lower left')
    for ext in ('png', 'pdf'):fig.savefig(OUT / f'figure_history_phase_portraits.{ext}')
    plt.close(fig)


def main():
    solver = load_solver()
    p = solver.Params(m0=.35, eps=(.5, .3), alpha_s=(1., 1.), beta_s=(1., 1.),
                      alpha_d=(0., 0.), beta_d=(0., 0.), charging=(.5, .5),
                      B1=(.001, .0002), B2=(.004, .001), cap=(.02, .004),
                      mode='pressure', forcing=10., kinetics='multistage')
    runs = {}; metrics = []; start = time.perf_counter()
    for n in GRIDS:
        for case in CASES:
            runstart = time.perf_counter()
            r = solver.simulate(p, n=n, tend=END, inlet=inlet_for(case), record_dt=.025)
            runs[(n, case)] = r
            m = save_run(case, n, r, p, solver)
            m['wall_seconds'] = time.perf_counter() - runstart
            metrics.append(m)
            print(json.dumps({'n': n, 'case': case, 'seconds': m['wall_seconds'],
                              'B_final': m['final_mean_retained'], 'flow': m['final_flow'],
                              'q_max': m['maximum_absolute_net_deposition_profile']}), flush=True)
    refinement = []
    for case in CASES:
        coarse = runs[(80, case)]; fine = runs[(160, case)]
        projected = np.asarray(fine['final']).reshape(2, 2, 80, 2).mean(axis=-1)
        refinement.append({'case': case,
            'maximum_final_state_difference_80_vs_projected_160': float(abs(coarse['final'] - projected).max()),
            'absolute_final_flow_difference': float(abs(coarse['u'][-1] - fine['u'][-1])),
            'maximum_recorded_flow_difference': float(abs(np.interp(fine['time'], coarse['time'], coarse['u']) - fine['u']).max()),
            'maximum_recorded_retained_difference': float(max(abs(np.interp(fine['time'], coarse['time'], coarse['retained'][:, i]) - fine['retained'][:, i]).max() for i in range(2)))})
    output = {
        'scope': 'Illustrative finite-inlet-history selection of interior equilibria from the same clean initial state; not experimental validation or a proof of coupled PDE stability.',
        'parameters': asdict(p), 'time_axis': 'Original t*=U0*t/L, not injected pore volumes',
        'preconditioning_end_t_star': PRE_END, 'final_t_star': END,
        'common_final_inlet': FINAL_FEED.tolist(), 'high_preconditioning_inlet': HIGH_FEED.tolist(),
        'initial_condition': 'All independent mobile/retained particle saturations zero; clean-water-filled original pores.',
        'same_final_forcing': 'Dimensionless pressure drop 10 for every history.',
        'interior_local_equilibria': {'low': FINAL_FEED.tolist(), 'high': np.sqrt(np.asarray(p.cap) * FINAL_FEED).tolist()},
        'all_results': metrics, 'grid_check': refinement, 'wall_seconds': time.perf_counter() - start,
        'interpretation_limits': [
            'The four studied equilibria are interior low/high combinations; absorbing capacity states are additional possible boundary states.',
            'Preconditioning protocols introduce different particle inputs and select different retained inventories; no matched-inventory hydraulic advantage is claimed.',
            'Local source phase portraits use prescribed final concentration; observed full-column trajectories do not constitute a general PDE stability theorem.',
            'The source law, original seven-saturation framework, fixed deposit porosities and hydraulic relations are unchanged.',
        ],
    }
    (OUT / 'history_results.json').write_text(json.dumps(output, indent=2) + '\n')
    plot_runs(runs, p)
    print(json.dumps({'completed': True, 'seconds': output['wall_seconds']}), flush=True)


if __name__ == '__main__':
    main()
