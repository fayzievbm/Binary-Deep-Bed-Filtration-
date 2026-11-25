#!/usr/bin/env python3
"""Uniform open-branch equilibria of the unchanged binary filtration model.

Run from any directory: python /path/to/run_equilibria.py
Dependencies: numpy, matplotlib, mpmath. No other project files are read.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
import mpmath as mp
import numpy as np


BRANCHES = ("LL", "LH", "HL", "HH")
EPS = np.array([0.5, 0.3])
W = 1.0 / (1.0 - EPS)
B1 = np.array([0.001, 0.0002])
B2 = np.array([0.004, 0.001])
CAP = np.array([0.02, 0.004])
GAMMA = np.array([0.002, 0.0005])
MU = 1.0 + 2.5 * GAMMA.sum()
BETA_D = np.array([0.05, 0.05])
ALPHA_D_EXAMPLE = np.array([0.1, 0.1])
D_MAX = float(np.dot(W, CAP))
KAPPA_MIN = float((1.0 - np.sqrt(D_MAX)) ** 3)
CASES = (
    {"name": "pure_dynamic", "mode": "head", "forcing": 10.0,
     "alpha_s": 0.0, "beta_s": 0.0, "title": "Purely dynamic rates"},
    {"name": "mixed_head", "mode": "head", "forcing": 10.0,
     "alpha_s": 0.1, "beta_s": 0.1, "title": r"Mixed rates: $P=10$"},
    {"name": "mixed_flow", "mode": "flow", "forcing": 10.0,
     "alpha_s": 0.1, "beta_s": 0.1, "title": r"Mixed rates: $u_f=10$"},
)


def rates_ratio(kappa, alpha_d, case):
    """r_i=alpha_i gamma_i/beta_i; kappa may have shape (...)."""
    kappa = np.asarray(kappa)
    if case["mode"] == "head":
        numerator = GAMMA * (case["alpha_s"] + alpha_d * case["forcing"] * kappa[..., None] / MU)
        return numerator / (case["beta_s"] + BETA_D * case["forcing"])
    numerator = GAMMA * (case["alpha_s"] + alpha_d * case["forcing"]) * kappa[..., None]
    denominator = case["beta_s"] * kappa[..., None] + BETA_D * case["forcing"] * MU
    return numerator / denominator


def deposits(kappa, alpha_d, branch, case):
    ratio = rates_ratio(kappa, alpha_d, case)
    return np.stack([ratio[..., i] if branch[i] == "L" else np.sqrt(CAP[i] * ratio[..., i])
                     for i in range(2)], axis=-1)


def scalar_residual(kappa, alpha_d, branch, case):
    return np.sum(W * deposits(kappa, alpha_d, branch, case), axis=-1) - (1 - np.cbrt(kappa)) ** 2


def candidates(alpha_d, branch, case):
    """Vector bisection; branch admissibility is checked separately."""
    shape = np.asarray(alpha_d).shape[:-1]
    lo = np.zeros(shape)
    hi = np.ones(shape)
    g_lo = scalar_residual(lo, alpha_d, branch, case)
    g_hi = scalar_residual(hi, alpha_d, branch, case)
    interior_root = (g_lo < 0) & (g_hi > 0)
    for _ in range(60):
        mid = lo + (hi - lo) / 2
        positive = scalar_residual(mid, alpha_d, branch, case) > 0
        hi = np.where(positive, mid, hi)
        lo = np.where(positive, lo, mid)
    root = lo + (hi - lo) / 2
    dep = deposits(root, alpha_d, branch, case)
    lower = np.array([B1[i] if branch[i] == "L" else B2[i] for i in range(2)])
    upper = np.array([B2[i] if branch[i] == "L" else CAP[i] for i in range(2)])
    margin = np.minimum(dep - lower, upper - dep)
    accepted = interior_root & np.all(margin > 0, axis=-1) & (root > KAPPA_MIN)
    return root, dep, accepted, margin, hi - lo


def boundary_xy(branch, species, endpoint, case, samples=401):
    """Exact branch-rectangle edge parameterization, evaluated in doubles."""
    other = 1 - species
    lower = B1[other] if branch[other] == "L" else B2[other]
    upper = B2[other] if branch[other] == "L" else CAP[other]
    dep = np.zeros((samples, 2))
    dep[:, species] = endpoint
    dep[:, other] = np.linspace(lower, upper, samples)
    kap = (1 - np.sqrt(dep @ W)) ** 3
    rho = np.stack([dep[:, i] if branch[i] == "L" else dep[:, i] ** 2 / CAP[i]
                    for i in range(2)], axis=-1)
    if case["mode"] == "head":
        d = (rho * (case["beta_s"] + BETA_D * case["forcing"]) - GAMMA * case["alpha_s"]) / (
            BETA_D * case["forcing"] * kap[:, None])
    else:
        d = rho / kap[:, None] + (rho * case["beta_s"] - GAMMA * case["alpha_s"]) / (
            BETA_D * MU * case["forcing"])
    return d[:, 0] / B2[0], d[:, 1] / B2[1]


def high_precision_enclosure(branch, case, displayed_kappa):
    """A high-precision bracket with outward interval signs at both endpoints.

    The input constants are their exact documented decimal values. Endpoints
    enclose a unique zero by the analytical monotonicity result. This certifies
    these example roots, not the complete parameter-map grid.
    """
    mp.mp.dps = 100
    mp.iv.dps = 80

    def eval_g(k, ctx):
        one = ctx.mpf(1)
        gamma = [ctx.mpf("0.002"), ctx.mpf("0.0005")]
        cap = [ctx.mpf("0.02"), ctx.mpf("0.004")]
        w = [ctx.mpf(2), one / ctx.mpf("0.7")]
        mu = one + ctx.mpf("2.5") * sum(gamma)
        ad = ctx.mpf("0.1")
        bd = ctx.mpf("0.05")
        a_s = ctx.mpf(str(case["alpha_s"]))
        b_s = ctx.mpf(str(case["beta_s"]))
        force = ctx.mpf(str(case["forcing"]))
        ans = ctx.mpf(0)
        for i in range(2):
            if case["mode"] == "head":
                r = gamma[i] * (a_s + ad * force * k / mu) / (b_s + bd * force)
            else:
                r = gamma[i] * (a_s + ad * force) * k / (b_s * k + bd * force * mu)
            b = r if branch[i] == "L" else ctx.sqrt(cap[i] * r)
            ans += w[i] * b
        return ans - (one - k ** (one / ctx.mpf(3))) ** 2

    lo, hi = mp.mpf("1e-90"), mp.mpf(1)
    assert eval_g(lo, mp.mp) < 0 < eval_g(hi, mp.mp)
    for _ in range(260):
        mid = (lo + hi) / 2
        if eval_g(mid, mp.mp) > 0:
            hi = mid
        else:
            lo = mid
    center = (lo + hi) / 2
    radius = mp.mpf("1e-45")
    lo_text = mp.nstr(center - radius, 70)
    hi_text = mp.nstr(center + radius, 70)
    g_lo_iv = eval_g(mp.iv.mpf(lo_text), mp.iv)
    g_hi_iv = eval_g(mp.iv.mpf(hi_text), mp.iv)
    lower_sign = bool(g_lo_iv < 0)
    upper_sign = bool(g_hi_iv > 0)
    assert lower_sign and upper_sign
    fp = mp.mpf(float(displayed_kappa))
    fp_bound = max(abs(fp - mp.mpf(lo_text)), abs(fp - mp.mpf(hi_text)))
    return {
        "kappa_high_precision": mp.nstr(center, 60),
        "kappa_bracket_lower": lo_text,
        "kappa_bracket_upper": hi_text,
        "bracket_width_upper_bound": "2.0000000000000000000000001e-45",
        "lower_endpoint_interval_G_strictly_negative": lower_sign,
        "upper_endpoint_interval_G_strictly_positive": upper_sign,
        "double_kappa_absolute_error_upper_bound": mp.nstr(fp_bound * (1 + mp.mpf("1e-20")), 35),
        "interval_arithmetic_decimal_precision": 80,
    }


def example_rows():
    rows = []
    certificates = {}
    for case in CASES:
        for branch in BRANCHES:
            kap, dep, accepted, margin, width = candidates(ALPHA_D_EXAMPLE, branch, case)
            kap = float(kap)
            assert bool(accepted), (case["name"], branch, dep)
            if case["mode"] == "head":
                flow = case["forcing"] * kap / MU
                pressure = case["forcing"]
            else:
                flow = case["forcing"]
                pressure = flow * MU / kap
            alpha = case["alpha_s"] + ALPHA_D_EXAMPLE * flow
            beta = case["beta_s"] + BETA_D * pressure
            q = np.array([alpha[i] * GAMMA[i] - beta[i] * dep[i] if branch[i] == "L"
                          else alpha[i] * CAP[i] * GAMMA[i] / dep[i] - beta[i] * dep[i]
                          for i in range(2)])
            phi = 1 - float(np.dot(W, dep))
            mobile = phi * GAMMA
            eta = W - 1
            saturations = np.array([phi - mobile.sum(), mobile[0], mobile[1],
                                    eta[0] * dep[0], eta[1] * dep[1], dep[0], dep[1]])
            row = {
                "case": case["name"], "branch": branch,
                "kappa": kap, "B_1": float(dep[0]), "B_2": float(dep[1]),
                "phi": phi, "u": float(flow), "pressure_drop": float(pressure),
                "pressure_drop_at_u10_same_state": float(10 * MU / kap),
                "u_at_P10_same_state": float(10 * kap / MU),
                "q_1": float(q[0]), "q_2": float(q[1]),
                "scalar_equilibrium_residual": float(scalar_residual(kap, ALPHA_D_EXAMPLE, branch, case)),
                "constitutive_kappa_residual": float(kap - (1 - np.sqrt(np.dot(W, dep))) ** 3),
                "darcy_residual": float(pressure - flow * MU / kap),
                "species_1_lower_branch_margin": float(dep[0] - (B1[0] if branch[0] == "L" else B2[0])),
                "species_1_upper_branch_margin": float((B2[0] if branch[0] == "L" else CAP[0]) - dep[0]),
                "species_2_lower_branch_margin": float(dep[1] - (B1[1] if branch[1] == "L" else B2[1])),
                "species_2_upper_branch_margin": float((B2[1] if branch[1] == "L" else CAP[1]) - dep[1]),
                "minimum_branch_margin": float(margin.min()),
                "minimum_saturation": float(saturations.min()),
                "saturation_sum_residual": float(saturations.sum() - 1),
                "derived_reserve_z": float(phi * (1 - np.dot(W, GAMMA))),
                "double_bisection_bracket_width": float(width),
            }
            certificate_key = case["name"] + ":" + branch
            certificates[certificate_key] = high_precision_enclosure(branch, case, kap)
            rows.append(row)
    return rows, certificates


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_maps(output_dir, n, axis_max):
    axis = np.linspace(0, axis_max, n)
    x, y = np.meshgrid(axis, axis)
    d = np.stack([x * B2[0], y * B2[1]], axis=-1)
    alpha_d = d * BETA_D * MU / GAMMA
    maps = []
    grid_summaries = {}
    for case in CASES:
        count = np.zeros(x.shape, dtype=np.uint8)
        mask = np.zeros(x.shape, dtype=np.uint8)
        root_width_max = 0.0
        min_boundary_margin = np.inf
        for index, branch in enumerate(BRANCHES):
            _, _, accepted, margin, width = candidates(alpha_d, branch, case)
            count += accepted.astype(np.uint8)
            mask += accepted.astype(np.uint8) * (1 << index)
            root_width_max = max(root_width_max, float(width.max()))
            min_boundary_margin = min(min_boundary_margin, float(np.min(np.abs(margin))))
        maps.append({"count": count, "mask": mask})
        unique, number = np.unique(count, return_counts=True)
        grid_summaries[case["name"]] = {
            "grid_point_counts_by_interior_multiplicity": {str(int(k)): int(v) for k, v in zip(unique, number)},
            "total_grid_points": int(x.size),
            "maximum_double_root_bracket_width": root_width_max,
            "minimum_absolute_branch_margin_over_grid_candidates": min_boundary_margin,
            "boundary_policy": "Strict open-branch inequalities; exact threshold/cap states excluded. Grid classification is floating-point, not interval-certified.",
        }
        path = output_dir / (case["name"] + "_regime_map.csv")
        with path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["d1_over_T1", "d2_over_T2", "interior_branch_count", "branch_bitmask_LL1_LH2_HL4_HH8"])
            for xv, yv, cv, mv in zip(x.ravel(), y.ravel(), count.ravel(), mask.ravel()):
                writer.writerow([f"{xv:.12g}", f"{yv:.12g}", int(cv), int(mv)])

    count_difference = maps[1]["count"] != maps[2]["count"]
    family_difference = maps[1]["mask"] != maps[2]["mask"]
    differing = np.argwhere(count_difference)
    representative = None
    if differing.size:
        # Choose a positive-axis point well inside the domain with maximal count
        # discrepancy; the selection rule is declared and uses no PDE trajectory.
        discrepancy = np.abs(maps[1]["count"].astype(int) - maps[2]["count"].astype(int))
        max_delta = int(discrepancy.max())
        chosen = np.argwhere((discrepancy == max_delta) & (x > 0.1) & (y > 0.1))
        if chosen.size == 0:
            chosen = np.argwhere(discrepancy == max_delta)
        center_distance = ((chosen[:, 0] - (n - 1) / 2) ** 2 + (chosen[:, 1] - (n - 1) / 2) ** 2)
        iy, ix = chosen[np.argmin(center_distance)]
        representative = {
            "d1_over_T1": float(x[iy, ix]), "d2_over_T2": float(y[iy, ix]),
            "head_branch_count": int(maps[1]["count"][iy, ix]),
            "flow_branch_count": int(maps[2]["count"][iy, ix]),
            "head_branches": [b for j, b in enumerate(BRANCHES) if int(maps[1]["mask"][iy, ix]) & (1 << j)],
            "flow_branches": [b for j, b in enumerate(BRANCHES) if int(maps[2]["mask"][iy, ix]) & (1 << j)],
        }
        rep_alpha_d = np.array([x[iy, ix], y[iy, ix]]) * B2 * BETA_D * MU / GAMMA
        rep_rows = []
        for case in CASES[1:]:
            for branch in BRANCHES:
                root, dep, accepted, margin, width = candidates(rep_alpha_d, branch, case)
                rep_rows.append({
                    "case": case["name"], "branch": branch,
                    "alpha_d_1": float(rep_alpha_d[0]), "alpha_d_2": float(rep_alpha_d[1]),
                    "kappa_candidate": float(root), "B_1_candidate": float(dep[0]), "B_2_candidate": float(dep[1]),
                    "accepted_interior_equilibrium": bool(accepted),
                    "species_1_signed_branch_margin": float(margin[0]),
                    "species_2_signed_branch_margin": float(margin[1]),
                    "double_bisection_bracket_width": float(width),
                })
        write_csv(output_dir / "mixed_mode_difference_example.csv", rep_rows)
        representative["candidate_details_csv"] = "mixed_mode_difference_example.csv"
        representative["minimum_absolute_candidate_branch_margin"] = min(
            min(abs(r["species_1_signed_branch_margin"]), abs(r["species_2_signed_branch_margin"])) for r in rep_rows)
    comparison = {
        "different_interior_counts_grid_points": int(count_difference.sum()),
        "different_admissible_branch_sets_grid_points": int(family_difference.sum()),
        "same_admissible_branch_sets_grid_points": int((~family_difference).sum()),
        "representative_count_difference": representative,
        "interpretation": "Parameter-grid counts only, not probabilities, solution selection, stability, or transient reachability.",
    }

    cmap = ListedColormap(["#f0f0f0", "#c6dbef", "#6baed6", "#756bb1", "#efbd53"])
    norm = BoundaryNorm(np.arange(-0.5, 5.5, 1), cmap.N)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelsize": 12, "axes.titlesize": 12,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, all_axes = plt.subplots(2, 2, figsize=(9, 8.7))
    axes = list(all_axes.ravel()[:3])
    legend_ax = all_axes.ravel()[3]
    legend_ax.axis("off")
    line_styles = {"LL": "-", "LH": "--", "HL": ":", "HH": "-."}
    marker = ALPHA_D_EXAMPLE * GAMMA / (BETA_D * MU * B2)
    for panel, (ax, case, regime) in enumerate(zip(axes, CASES, maps)):
        pic = ax.pcolormesh(x, y, regime["count"], cmap=cmap, norm=norm, shading="nearest", rasterized=True)
        for branch in BRANCHES:
            for species in range(2):
                endpoints = (B1[species], B2[species]) if branch[species] == "L" else (B2[species], CAP[species])
                for endpoint in endpoints:
                    bx, by = boundary_xy(branch, species, endpoint, case)
                    ax.plot(bx, by, color="#222222", lw=1.05, ls=line_styles[branch])
        ax.plot(marker[0], marker[1], "*", ms=13, mfc="white", mec="black", mew=0.9, zorder=5)
        ax.set(xlim=(0, axis_max), ylim=(0, axis_max), xlabel=r"$d_1/T_1$", ylabel=r"$d_2/T_2$",
               title="(" + "abc"[panel] + ") " + case["title"])
        ax.set_aspect("equal")
        ax.set_xticks(np.arange(0, axis_max + 0.01, 0.5))
        ax.set_yticks(np.arange(0, axis_max + 0.01, 0.5))
        ax.grid(alpha=0.16, linewidth=0.5)
    fig.subplots_adjust(left=0.10, right=0.98, top=0.90, bottom=0.08, wspace=0.25, hspace=0.33)
    pos = legend_ax.get_position()
    cax = fig.add_axes([pos.x0 + 0.03, pos.y0 + 0.04, pos.width - 0.06, 0.025])
    cb = fig.colorbar(pic, cax=cax, orientation="horizontal", ticks=range(5))
    cb.set_label("Interior branch combinations", fontsize=10)
    handles = [Line2D([0], [0], color="#222222", ls=line_styles[b], lw=1.3, label=b + " boundary") for b in BRANCHES]
    handles.append(Line2D([0], [0], marker="*", color="none", markerfacecolor="white", markeredgecolor="black", markersize=10, label="Example"))
    legend_ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 1.0), ncol=2, frameon=False, fontsize=10)
    legend_ax.text(0.08, 0.63, "Colors count admissible homogeneous\nopen-branch equilibria.\n\nZero excludes only this interior family.\nThreshold and absorbing-capacity\nequilibria are not counted.",
                   transform=legend_ax.transAxes, va="top", fontsize=10, linespacing=1.3)
    fig.suptitle("Coupled Darcy feedback and equilibrium multiplicity", y=0.98, fontsize=14)
    for ext in ("png", "pdf"):
        fig.savefig(output_dir / ("equilibrium_regime_map." + ext), dpi=220 if ext == "png" else None)
    plt.close(fig)
    return grid_summaries, comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--grid-size", type=int, default=321)
    parser.add_argument("--axis-max", type=float, default=2.0)
    args = parser.parse_args()
    if args.grid_size < 3 or args.axis_max <= 0:
        parser.error("grid-size must be at least 3 and axis-max positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows, certificates = example_rows()
    write_csv(args.output_dir / "example_equilibria.csv", rows)
    summaries, comparison = make_maps(args.output_dir, args.grid_size, args.axis_max)
    summary = {
        "model_scope": "Only deposition kinetics differ from the 2017 model. These are uniform steady states, not transient simulations or PDE stability results.",
        "parameter_provenance": {
            "original_model": "Makhmudov, Saidullaev and Khuzhayorov, Fluid Dynamics 52 (2017), 299–308, doi:10.1134/S0015462817020132",
            "porosities": "Existing manuscript baseline, epsilon=(0.5,0.3).",
            "thresholds_and_capacities": "Existing manuscript illustrative Table 1; not fitted experimental parameters.",
            "feed": "Existing below-capacity history-study final feed gamma=(0.002,0.0005).",
            "pure_dynamic_example": "Existing dynamic-rate coefficients alpha_di=0.1, beta_di=0.05, now combined with that lower final feed for steady-state analysis.",
            "mixed_example": "Original static coefficients alpha_si=beta_si=0.1, as in the existing static-rate sweep; this steady-state map is new.",
            "map_axes": "Fixed gamma, mu and beta_di; alpha_di varies so d_i/T_i=alpha_di*gamma_i/(beta_di*mu*T_i).",
        },
        "parameters": {
            "epsilon": EPS.tolist(), "w": W.tolist(), "B1": B1.tolist(), "T_B2": B2.tolist(),
            "C": CAP.tolist(), "gamma": GAMMA.tolist(), "mu_gamma": float(MU),
            "beta_d": BETA_D.tolist(), "example_alpha_d": ALPHA_D_EXAMPLE.tolist(),
            "D_max": D_MAX, "kappa_capacity_lower_bound": KAPPA_MIN,
            "weighted_inlet_concentration": float(np.dot(W, GAMMA)),
            "cases": CASES,
        },
        "grid": {"size_each_axis": args.grid_size, "axis_min": 0, "axis_max": args.axis_max,
                 "axis_spacing": args.axis_max / (args.grid_size - 1),
                 "membership_bits": {b: 1 << j for j, b in enumerate(BRANCHES)}},
        "example_states": rows,
        "example_root_certificates": certificates,
        "example_max_abs_source_residual": max(max(abs(r["q_1"]), abs(r["q_2"])) for r in rows),
        "example_max_abs_constitutive_residual": max(abs(r["constitutive_kappa_residual"]) for r in rows),
        "example_minimum_branch_margin": min(r["minimum_branch_margin"] for r in rows),
        "example_max_double_root_absolute_error_bound": max(float(c["double_kappa_absolute_error_upper_bound"]) for c in certificates.values()),
        "regime_maps": summaries,
        "mixed_mode_comparison": comparison,
        "limitations": [
            "No PDE stability, transient branch selection, clean-bed reachability, or experimental identification is inferred.",
            "Counts exclude charging-floor, repelling-interface and absorbing-capacity states and nonuniform equilibria.",
            "Mapped boundaries are exact parameterizations evaluated in floating point; strict grid classifications at boundaries are not certified.",
            "Interval root certificates apply to the twelve example states only, with documented decimal input constants.",
            "Grid-point fractions depend on the chosen parameter rectangle and are not probabilities.",
        ],
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir),
                      "example_max_abs_source_residual": summary["example_max_abs_source_residual"],
                      "example_minimum_branch_margin": summary["example_minimum_branch_margin"],
                      "example_max_double_root_absolute_error_bound": summary["example_max_double_root_absolute_error_bound"],
                      "regime_maps": summaries, "mixed_mode_comparison": comparison}, indent=2))


if __name__ == "__main__":
    main()
