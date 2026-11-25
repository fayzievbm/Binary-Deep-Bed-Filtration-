"""Necessary hydraulic compatibility under the original 2017 permeability law.

No kinetic, permeability or observation parameter is fitted or changed. A and
the multiplier boundaries below are mathematical consequences or deterministic
sensitivity coordinates, not new constitutive coefficients or confidence sets.
"""
from __future__ import annotations

import csv
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / 'inputs'


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(name, rows):
    with (ROOT / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def P(s):
    return 2 * s ** 3 - 8 * s ** 2 + 12 * s - 3


def H(s):
    """Stable expression for (((1-s)**-3)-1)/s**2, 0<s<1."""
    return (s * s - 3 * s + 3) / (s * (1 - s) ** 3)


def Kbound(Bmean, eps, constant):
    return 1 / (1 + constant * np.asarray(Bmean) / (1 - eps))


def decimal_constants():
    with localcontext() as context:
        context.prec = 80
        s = Decimal('0.31')
        for _ in range(40):
            f = 2 * s ** 3 - 8 * s ** 2 + 12 * s - 3
            fp = 6 * s ** 2 - 16 * s + 12
            s -= f / fp
        A = (s * s - 3 * s + 3) / (s * (1 - s) ** 3)
        return str(s), str(A)


def plot_results(Bmean, Knom, eps_values, A):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12,
                         'axes.labelsize': 12, 'axes.titlesize': 12,
                         'xtick.labelsize': 11, 'ytick.labelsize': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'savefig.dpi': 400,
                         'lines.linewidth': 1.9})
    colors = ['#2267A0', '#22885B', '#D16A24', '#8A56A3']
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.6))
    fig.subplots_adjust(left=.09, right=.98, bottom=.27, top=.88, wspace=.30)
    B = np.linspace(0, .01, 401)
    loading_multiplier = np.linspace(0, 1.25, 251)
    for eps, color in zip(eps_values, colors):
        label = rf'$\varepsilon={eps:g}$'
        axes[0].plot(1e3 * B, Kbound(B, eps, A), color=color, label=label)
        axes[1].plot(loading_multiplier,
                     Kbound(loading_multiplier * Bmean, eps, A) / Knom,
                     color=color)
    axes[0].scatter([1e3 * Bmean], [Knom], marker='D', s=42,
                    color='black', zorder=5)
    axes[0].axvline(1e3 * Bmean, color='#92989E', lw=.9, ls=':')
    axes[0].axhline(Knom, color='#92989E', lw=.9, ls=':')
    axes[0].annotate('Nominal mapped pair', xy=(1e3 * Bmean, Knom),
                     xytext=(6.0, 1.005), textcoords='data', fontsize=10,
                     arrowprops={'arrowstyle': '-', 'color': '#555555', 'lw': .8})
    axes[0].set_xlabel(r'Mean solid loading $10^3\overline{B}$')
    axes[0].set_ylabel(r'Upper bound on $k_{\mathrm{eff}}/k_0$')
    axes[0].set_ylim(.55, 1.035)
    axes[0].set_xlim(0, 10)
    axes[0].set_title('(a) Original-law loading bound', loc='left', pad=11)
    axes[1].scatter([1], [1], marker='D', s=42, color='black', zorder=5)
    axes[1].axvline(1, color='#92989E', lw=.9, ls=':')
    axes[1].axhline(1, color='#92989E', lw=.9, ls=':')
    axes[1].set_xlabel(r'Loading multiplier $\lambda_B$')
    axes[1].set_ylabel(r'Maximum ratio multiplier $\lambda_K$')
    axes[1].set_xlim(0, 1.25)
    axes[1].set_ylim(.65, 1.08)
    axes[1].set_title('(b) Deterministic normalization sensitivity', loc='left', pad=11)
    axes[1].text(.03, .68, 'Below a curve: not excluded by that bound',
                 fontsize=9, color='#444B53')
    for ax in axes:
        ax.grid(axis='y', alpha=.2)
        ax.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    handles.append(Line2D([0], [0], color='black', marker='D', linestyle='none',
                          markersize=5, label='Nominal pair / multipliers (1, 1)'))
    labels.append('Nominal pair / multipliers (1, 1)')
    fig.legend(handles, labels, loc='lower center', ncol=3, fontsize=10,
               bbox_to_anchor=(.5, .045))
    fig.text(.09, .01, 'Necessary bounds and deterministic scenarios; no experimental error bars or confidence regions.',
             fontsize=9, color='#444B53')
    for ext in ('png', 'pdf'):
        fig.savefig(ROOT / f'hydraulic_compatibility.{ext}')
    plt.close(fig)


def main():
    parameters = json.loads((INPUTS / 'parameters.json').read_text())
    source_files = [INPUTS / name for name in (
        'parameters.json', 'data_retention_table1.csv',
        'data_retention_profile_table1_flow.csv')]
    input_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    rows = read_csv(INPUTS / 'data_retention_table1.csv')
    flow = [r for r in rows if r['mode'] == 'constant_flow']
    if len(flow) != 1:raise ValueError('Expected one selected flow-retention endpoint.')
    flow = flow[0]
    cin = float(flow['feed_g_L']) / parameters['particle_material_density_kg_m3']
    PV = float(flow['end_pore_volumes'])
    retained_fraction = float(flow['retained_fraction_of_injected'])
    Bmean = cin * PV * retained_fraction
    profile = read_csv(INPUTS / 'data_retention_profile_table1_flow.csv')
    widths = np.array([float(r['x_right_m']) - float(r['x_left_m']) for r in profile])
    weights = widths / widths.sum()
    Bbins = np.array([float(r['B_mean_in_bin_per_original_pore_volume']) for r in profile])
    Bprofile = float(np.dot(weights, Bbins))
    if abs(Bprofile - Bmean) > 1e-14:raise ValueError('Profile and total-retention input mismatch.')
    root = brentq(P, 0., 1., xtol=1e-14, rtol=1e-14)
    A = float(H(root))
    Dstar = root * root
    precise_root, precise_A = decimal_constants()
    K = float(parameters['nominal_hydraulic_ratio'])
    eps_values = parameters['deposit_porosity_values']
    inversion = []
    for eps in eps_values:
        Bmax = (1 - eps) * (1 / K - 1) / A
        ceiling = float(Kbound(Bmean, eps, A))
        profile_K = float(1 / np.dot(weights, (1 - np.sqrt(Bbins / (1 - eps))) ** -3))
        inversion.append({'deposit_porosity_assumed': eps,
                          'nominal_mean_retained_saturation': Bmean,
                          'nominal_mapped_hydraulic_ratio': K,
                          'original_law_upper_bound_at_nominal_loading': ceiling,
                          'piecewise_constant_section_profile_ratio': profile_K,
                          'maximum_loading_not_excluded_by_bound_at_nominal_ratio': Bmax,
                          'maximum_fraction_of_nominal_loading_not_excluded': Bmax / Bmean,
                          'maximum_percent_of_nominal_loading_not_excluded': 100 * Bmax / Bmean,
                          'maximum_ratio_multiplier_not_excluded_at_nominal_loading': ceiling / K,
                          'nominal_point_satisfies_necessary_bound': bool(K <= ceiling)})
    write_csv('compatibility_inversion.csv', inversion)
    curves = []
    for eps in eps_values:
        for B in np.linspace(0, .01, 401):
            curves.append({'deposit_porosity_assumed': eps, 'mean_retained_saturation': B,
                           'original_law_upper_bound': float(Kbound(B, eps, A))})
    write_csv('compatibility_curves.csv', curves)
    boundaries = []
    for eps in eps_values:
        for lb in np.linspace(0, 1.25, 251):
            upper = float(Kbound(lb * Bmean, eps, A))
            boundaries.append({'deposit_porosity_assumed': eps, 'loading_multiplier_lambda_B': lb,
                               'mean_retained_saturation_used': lb * Bmean,
                               'hydraulic_ratio_upper_bound': upper,
                               'maximum_hydraulic_ratio_multiplier_lambda_K': upper / K,
                               'interpretation': 'deterministic necessary-bound boundary, not CI or fitted correction'})
    write_csv('normalization_boundaries.csv', boundaries)
    matrix = []
    for eps in eps_values:
        for lb in (0., .25, .5, .75, 1., 1.25):
            upper = float(Kbound(lb * Bmean, eps, A))
            for lk in (.8, .85, .9, .95, 1., 1.05):
                mapped = lk * K
                matrix.append({'deposit_porosity_assumed': eps, 'loading_multiplier_lambda_B': lb,
                               'hydraulic_ratio_multiplier_lambda_K': lk,
                               'mapped_hydraulic_ratio': mapped,
                               'original_law_upper_bound': upper,
                               'not_excluded_by_necessary_bound': bool(mapped <= upper),
                               'interpretation': 'deterministic sensitivity, not a confidence probability'})
    write_csv('normalization_matrix.csv', matrix)
    samples = np.unique(np.concatenate(([0., Dstar],
                        np.geomspace(1e-20, .999999999, 30000),
                        np.linspace(0, .999999999, 30000))))
    resistance = (1 - np.sqrt(samples)) ** -3
    slack = resistance - 1 - A * samples
    opt = minimize_scalar(H, bounds=(1e-8, 1 - 1e-8), method='bounded',
                          options={'xatol': 1e-14})
    checks = {'polynomial_residual_at_brent_root': abs(P(root)),
              'double_root_absolute_error_against_80_digit_decimal': abs(root - float(precise_root)),
              'double_A_absolute_error_against_80_digit_decimal': abs(A - float(precise_A)),
              'independent_minimizer_A_absolute_difference': abs(A - float(opt.fun)),
              'sampled_pointwise_support_minimum_slack': float(slack.min()),
              'sampled_pointwise_support_count': len(samples),
              'D_zero_slack_exactly': 0.,
              'D_star_slack_floating_point': float((1 - root) ** -3 - 1 - A * Dstar),
              'section_profile_vs_aggregate_mean_absolute_error': abs(Bprofile - Bmean),
              'maximum_inverse_forward_ratio_error': max(abs(float(Kbound(r['maximum_loading_not_excluded_by_bound_at_nominal_ratio'], r['deposit_porosity_assumed'], A)) - K) for r in inversion),
              'any_stronger_support_A_plus_1e_minus8_fails_at_D_star':
                  bool((1 - root) ** -3 - 1 - (A + 1e-8) * Dstar < 0),
              'input_hashes_unchanged': input_hashes == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}}
    if checks['sampled_pointwise_support_minimum_slack'] < -2e-12:raise AssertionError(checks)
    if not checks['input_hashes_unchanged']:raise AssertionError('Input mutation.')
    limits = [
        'These are necessary deterministic model implications, not statistical rejection or experimental confidence bounds.',
        'The comparison requires the same homogeneous-reference bed interval, loading endpoint, head normalization and relevant fluid properties.',
        'The nominal 0.94 ratio is rounded from the published low-concentration flow series; its exact match to the 83-PV recovered-mass endpoint is not independently verified.',
        'Reported initial conductivity variability is about 7%; it does not supply an uncertainty interval for the normalized ratio or justify choosing a correction multiplier.',
        'Hydraulic conductivity and intrinsic permeability ratios require density/viscosity mapping; the original dilute viscosity correction at the stated feed is below 0.03%.',
        'Deposit porosity is an assumed fixed constitutive input; it is not measured by these comparisons.',
        'A smaller loading or altered ratio multiplier that passes the bound is not a fitted or validated correction and does not establish full model compatibility.',
        'The affine support constant is optimal for the permeability function alone; capacity, profile and kinetic reachability restrictions can give stricter conditions.',
    ]
    result = {'original_model_doi': parameters['original_model_doi'],
              'experimental_source_doi': parameters['experimental_source_doi'],
              'input_sha256': input_hashes,
              'feed_volume_fraction': cin, 'pore_volumes': PV,
              'retained_fraction_of_injected': retained_fraction,
              'nominal_mean_retained_saturation': Bmean,
              'nominal_hydraulic_ratio': K,
              'nominal_hydraulic_ratio_role': parameters['nominal_hydraulic_ratio_role'],
              'support_constant_A': A, 'support_root_s': root,
              'support_contact_D_star': Dstar,
              'decimal_80_digit_root_s': precise_root, 'decimal_80_digit_A': precise_A,
              'reported_initial_conductivity_fractional_variability':
                  parameters['reported_initial_conductivity_plusminus_m_s'] / parameters['reported_initial_conductivity_m_s'],
              'inversion': inversion, 'numerical_precision_checks': checks,
              'limitations': limits,
              'deterministic_multiplier_definition': 'lambda_B multiplies nominal retained loading; lambda_K multiplies nominal mapped hydraulic ratio. Neither multiplier is estimated or assigned a probability.'}
    (ROOT / 'compatibility_results.json').write_text(json.dumps(result, indent=2) + '\n')
    plot_results(Bmean, K, eps_values, A)
    print(json.dumps({'A': A, 's_star': root, 'D_star': Dstar, 'Bmean': Bmean,
                      'inversion': inversion, 'precision_checks': checks}, indent=2))


if __name__ == '__main__':
    main()
