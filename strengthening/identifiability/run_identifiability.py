"""Descriptive identifiability screens using the original calibration subsets.

Run from an unpacked article package as:
  python strengthening/identifiability/run_identifiability.py
The sibling validation/ directory contains the unchanged scalar solver and data.
For a source checkout, work_validation/ is also discovered. No held-out response
is used, no published fit is overwritten, and no experimental noise model or
confidence interval is inferred. Python 3, NumPy, SciPy, Matplotlib and g++ needed.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
import numpy as np
from scipy.optimize import least_squares, minimize, minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent


def discover_validation(explicit=None):
    if explicit:
        candidates = [Path(explicit).resolve()]
    else:
        candidates = [p / name for p in [HERE, *HERE.parents]
                      for name in ("validation", "work_validation")]
    for candidate in candidates:
        if (candidate / "validation_solver.py").is_file() and (
                candidate / "gitis/fit_gitis_results.json").is_file():
            return candidate
    raise FileNotFoundError("Cannot find validation/; use --validation-root PATH")


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(v) for v in value]
    return value


def load_solver(validation, output):
    """Keep compilation/cache writes entirely inside the diagnostic output."""
    spec = importlib.util.spec_from_file_location(
        "identifiability_scalar_solver", validation / "validation_solver.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    cache = output / "kernel_cache"
    cache.mkdir(exist_ok=True)
    source = validation / "validation_kernel.cpp"
    target = cache / source.name
    if not target.exists() or sha256(target) != sha256(source):
        shutil.copyfile(source, target)
    module.ROOT = cache
    return module


def summarize(prediction, observed, result=None):
    residual = prediction - observed
    record = dict(SSE=float(residual @ residual),
                  RMSE=float(np.sqrt(np.mean(residual**2))))
    if result is not None:
        record.update(maximum_B=float(result["B_profile"].max()),
                      minimum_carrier=float(result["minimum_carrier"].min()),
                      maximum_balance_error=float(result["max_balance_error"].max()),
                      cap_events=int(result["cap_events"][-1]),
                      hold_events=int(result["hold_events"][-1]))
    return record


def svd_record(jacobian, labels, step, normalization, scheme):
    _, singular, vt = np.linalg.svd(jacobian, full_matrices=False)
    return dict(parameter_coordinates=labels, difference_step=step,
                finite_difference_scheme=scheme,
                response_scale_RMS_observed=normalization,
                singular_values=singular,
                relative_singular_values=singular / singular[0],
                condition_number=float(singular[0] / singular[-1])
                if singular[-1] > 0 else None,
                right_singular_vectors=vt,
                weakest_direction=vt[-1],
                normalized_jacobian=jacobian)


def log_screen(predictor, point, labels, observed, steps=(0.02, 0.01, 0.005)):
    scale = float(np.sqrt(np.mean(observed**2)))
    records = []
    for step in steps:
        columns = []
        for j in range(len(point)):
            increment = np.zeros(len(point))
            increment[j] = step
            columns.append((predictor(point + increment) -
                            predictor(point - increment)) / (2 * step * scale))
        records.append(svd_record(np.array(columns).T, labels, step, scale,
                                  "central differences in natural-log coordinates"))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-root")
    parser.add_argument("--output", type=Path, default=HERE)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    validation = discover_validation(args.validation_root)
    sources = ["validation_solver.py", "validation_kernel.cpp",
               "gitis/figure6_experimental_points.csv", "gitis/fit_gitis_results.json",
               "alem/data_breakthrough.csv", "alem/results/fit_results.json"]
    original_hashes = {name: sha256(validation / name) for name in sources}
    solver = load_solver(validation, output)
    Params, simulate = solver.Params, solver.simulate
    started = time.perf_counter()

    # Read/filter response data once. No Figure 7, deeper Gitis response or Alem
    # head, permeability, or retained-mass response enters an objective or screen.
    grows = [r for r in read_csv(validation / sources[2])
             if abs(float(r["depth_m"]) - 0.82) < 1e-8]
    arows = [r for r in read_csv(validation / sources[4])
             if r["mode"] == "constant_flow"]
    assert len(grows) == 9 and len(arows) == 26
    gfit = json.loads((validation / sources[3]).read_text())
    afit = json.loads((validation / sources[5]).read_text())
    gy = np.array([float(r["concentration_ratio"]) for r in grows])
    gt = np.array([float(r["time_h"]) for r in grows])
    ay = np.array([float(r["value"]) for r in arows])
    ax = np.array([float(r["pore_volumes"]) for r in arows])
    gtheta = np.array(gfit["fitted_theta_log"])
    gin = 1e-6
    ain = afit["concentration_volume_fraction"]
    apoint = Params(**afit["fits"]["multistage"]["parameters"])
    calls = {"gitis": 0, "alem": 0}

    def gpredict(theta=gtheta, n=160, cin=gin, b2factor=1., capfactor=1.,
                 full=False):
        calls["gitis"] += 1
        alpha, chi, ratio = np.exp(theta)
        p = Params(m0=.44, eps=.5, alpha_s=float(alpha), beta_s=0.,
                   charging=float(chi), B1=float(ratio * cin),
                   B2=.01 * cin / gin * b2factor,
                   cap=.02 * cin / gin * capfactor,
                   mode="rate", forcing=1., kinetics="multistage")
        result = simulate(p, sample=10 / 1.62 * (gt + .1), sample_axis="time",
                          cin=cin, depths=(.82 / 1.62,), n=n)
        prediction = result["c_normalized"][:, 0]
        return (prediction, result, p) if full else prediction

    def apredict(p=apoint, n=120, full=False):
        calls["alem"] += 1
        result = simulate(p, sample=ax, cin=ain, n=n)
        prediction = result["c_normalized"][:, 0]
        return (prediction, result) if full else prediction

    gbase, gr, gp = gpredict(full=True)
    abase, ar = apredict(full=True)
    gbase_metrics = summarize(gbase, gy, gr)
    abase_metrics = summarize(abase, ay, ar)
    print("Saved-point baselines reproduced; screening sensitivities", flush=True)

    gsvd = log_screen(gpredict, gtheta,
                      ["ln(alpha_s)", "ln(chi)", "ln(B1/c0)"], gy)
    atheta = np.log([apoint.alpha_s, apoint.cap])
    asvd = log_screen(lambda z: apredict(replace(
        apoint, alpha_s=float(np.exp(z[0])), cap=float(np.exp(z[1])))),
        atheta, ["ln(alpha_s)", "ln(capacity)"], ay)
    # A log(beta) derivative at beta~0 can look small merely because beta is tiny.
    # Complement it with a physically scaled, one-sided beta coordinate based on
    # the deposition equation m0*dB/dt = ... - beta*B: beta*t_end/m0=beta*PV_end.
    beta_scale = 1 / float(ax.max())
    amixed = []
    response_scale = float(np.sqrt(np.mean(ay**2)))
    for step in (0.02, 0.01, 0.005):
        j = np.array(asvd[(0.02, 0.01, 0.005).index(step)]["normalized_jacobian"])
        dbeta = (apredict(replace(apoint, beta_s=apoint.beta_s +
                                 step * beta_scale)) - abase) / (step * response_scale)
        amixed.append(svd_record(np.column_stack([j, dbeta]),
                      ["ln(alpha_s)", "ln(capacity)", "beta_s*PV_end"],
                      step, response_scale,
                      "central alpha/cap; forward beta from saved near-boundary point"))

    inactive = []
    for parameter, b2factor, capfactor in [
            ("baseline", 1., 1.), ("B2_half", .5, 1.),
            ("B2_times_1.5", 1.5, 1.), ("capacity_times_0.75", 1., .75),
            ("capacity_double", 1., 2.)]:
        prediction, result, p = gpredict(b2factor=b2factor,
                                       capfactor=capfactor, full=True)
        inactive.append(dict(variation=parameter, B2=p.B2, capacity=p.cap,
                             max_absolute_prediction_change=float(abs(prediction - gbase).max()),
                             **summarize(prediction, gy, result)))
    concentration_scale = []
    for cin in (1e-7, 1e-6, 1e-5):
        prediction, result, p = gpredict(cin=cin, full=True)
        concentration_scale.append(dict(
            assumed_c0=cin, B1=p.B1, B2=p.B2, capacity=p.cap,
            max_absolute_prediction_change=float(abs(prediction - gbase).max()),
            **summarize(prediction, gy, result)))

    # Descriptive objective profile: fix B1/c0, reoptimize only alpha and chi.
    # All profiles retain the declared clock, c0, porosity and other restrictions.
    print("Profiling Gitis B1/c0 using only 9 calibration symbols", flush=True)
    gprofiles = []
    for factor in (.5, .7, .85, 1., 1.2, 1.5, 2.):
        fixed = gtheta[2] + np.log(factor)
        def residual(z):
            return gpredict(np.r_[z, fixed]) - gy
        attempts = [least_squares(residual, start,
                bounds=(np.log([.1, .01]), np.log([100., 1.])),
                diff_step=1e-3, ftol=1e-8, xtol=1e-7, gtol=1e-8, max_nfev=70)
                for start in (gtheta[:2], np.log([8., .45]))]
        opt = min(attempts, key=lambda candidate: float(candidate.fun @ candidate.fun))
        prediction, result, p = gpredict(np.r_[opt.x, fixed], full=True)
        row = dict(B1_over_c0_factor=factor, B1_over_c0=float(np.exp(fixed)),
                   alpha_s=p.alpha_s, chi=p.charging,
                   **summarize(prediction, gy, result),
                   objective_ratio_to_saved_point=float(np.sum(opt.fun**2) / gbase_metrics["SSE"]),
                   optimizer_success=bool(opt.success), optimizer_nfev=int(opt.nfev),
                   optimizer_starts=2,
                   optimizer_total_nfev=int(sum(candidate.nfev for candidate in attempts)),
                   saved_point_start_SSE=float(attempts[0].fun @ attempts[0].fun),
                   alternative_start_SSE=float(attempts[1].fun @ attempts[1].fun))
        gprofiles.append(row)
        print("Gitis", factor, row["RMSE"], flush=True)

    # Alem: keep the same fixed early thresholds, explore capacity and retain
    # original alpha/beta search limits. This is a conditional ridge screen,
    # not a fit claiming a newly established material capacity.
    print("Profiling Alem attachment-capacity ridge at N120", flush=True)
    aprofiles = []
    for cap in (.01, np.sqrt(.0002), .02, np.sqrt(.0008), .04):
        guess = np.log10([apoint.alpha_s * apoint.cap / cap, apoint.beta_s])
        def objective(z):
            p = replace(apoint, alpha_s=float(10**z[0]),
                        beta_s=float(10**z[1]), cap=float(cap))
            residual = apredict(p) - ay
            return float(np.mean(residual**2))
        opt = minimize(objective, guess, method="Nelder-Mead",
                       bounds=[(-3., .5), (-7., .3)],
                       options=dict(maxiter=75, xatol=3e-4, fatol=2e-8))
        p = replace(apoint, alpha_s=float(10**opt.x[0]),
                    beta_s=float(10**opt.x[1]), cap=float(cap))
        prediction, result = apredict(p, full=True)
        aprofiles.append(dict(capacity=cap, alpha_s=p.alpha_s, beta_s=p.beta_s,
            alpha_times_capacity=p.alpha_s * cap,
            **summarize(prediction, ay, result),
            objective_ratio_to_saved_point=float(np.mean((prediction-ay)**2) /
                                                (abase_metrics["SSE"] / len(ay))),
            optimizer_success=bool(opt.success), optimizer_nfev=int(opt.nfev)))
        print("Alem capacity", cap, aprofiles[-1]["RMSE"], flush=True)

    print("Profiling one-sided Alem detachment, including beta=0", flush=True)
    beta_profiles = []
    for beta in (0., 1e-7, 1e-5, 1e-3, .01, .03):
        def objective(logalpha):
            residual = apredict(replace(apoint,
                alpha_s=float(np.exp(logalpha)), beta_s=beta)) - ay
            return float(np.mean(residual**2))
        opt = minimize_scalar(objective, method="bounded",
                    bounds=(np.log(1e-3), np.log(10**.5)),
                    options=dict(xatol=3e-4, maxiter=32))
        p = replace(apoint, alpha_s=float(np.exp(opt.x)), beta_s=beta)
        prediction, result = apredict(p, full=True)
        beta_profiles.append(dict(beta_s=beta, alpha_s=p.alpha_s,
            **summarize(prediction, ay, result),
            objective_ratio_to_saved_point=float(np.sum((prediction-ay)**2) /
                                                abase_metrics["SSE"]),
            optimizer_success=bool(opt.success), optimizer_nfev=int(opt.nfev)))
        print("Alem beta", beta, beta_profiles[-1]["RMSE"], flush=True)

    # Fixed-parameter checks on the screen, not an additional fit on any mesh.
    refinement = []
    for n in (120, 240):
        prediction, result = apredict(n=n, full=True)
        refinement.append(dict(dataset="Alem saved point", n=n,
            **summarize(prediction, ay, result)))
    for n in (160, 320):
        prediction, result, _ = gpredict(n=n, full=True)
        refinement.append(dict(dataset="Gitis saved point", n=n,
            **summarize(prediction, gy, result)))

    after_hashes = {name: sha256(validation / name) for name in sources}
    assert after_hashes == original_hashes, "Original input changed during diagnostic"
    result = dict(
        purpose="Descriptive conditional identifiability screen; no statistical confidence intervals",
        scope="Original 9 Gitis shallow-depth symbols and 26 Alem constant-flow concentration symbols only",
        held_out_responses_used=False,
        original_fits_replaced=False,
        noise_model=None,
        normalization="Jacobian of normalized prediction divided by RMS calibration observations",
        provenance=dict(input_sha256=original_hashes, input_hashes_unchanged=True,
                        script_sha256=sha256(Path(__file__)),
                        relative_validation_directory=validation.name),
        gitis=dict(n=160, calibration_n=9, baseline_parameters=asdict(gp),
                   baseline_metrics=gbase_metrics,
                   saved_objective_SSE_N160=gfit["main_multistage"]["objective_SSE"],
                   saved_theta_log=gtheta, sensitivity=gsvd,
                   inactive_threshold_screen=inactive,
                   concentration_scale_screen=concentration_scale,
                   B1_objective_profile=gprofiles,
                   fixed_constraints="beta=0, dynamic rates=0, c0=1e-6, eps=.5, clock=.1h, inactive B2/cap fixed"),
        alem=dict(n=120, calibration_n=26, baseline_parameters=asdict(apoint),
                  baseline_metrics=abase_metrics,
                  saved_fit_mean_squared_error_N120=afit["fits"]["multistage"]["fit_mean_squared_error_N120"],
                  log_alpha_capacity_sensitivity=asvd,
                  mixed_beta_sensitivity=amixed,
                  beta_scale=beta_scale,
                  capacity_objective_profile=aprofiles,
                  one_sided_beta_objective_profile=beta_profiles,
                  fixed_constraints="chi=.5, B1/c0=.001, B2/c0=.1, eps=.5, dynamic rates=0"),
        fixed_parameter_refinement=refinement,
        elapsed_seconds=time.perf_counter()-started,
        simulation_calls=calls,
        interpretation_limits=[
            "Singular values and directions depend on coordinates, scales, finite-difference step, mesh and threshold crossings.",
            "The tested finite-difference steps do not demonstrate Jacobian step convergence; these are finite-step screens, not converged derivative estimates.",
            "These local screens and conditional objective profiles are not global structural-identifiability proofs.",
            "An inactive threshold has no response influence as long as the same active branches are retained.",
            "c0 scaling is approximate, because original pore-volume and hydraulic couplings depend on absolute deposited mass.",
            "The aging source contains alpha_s*capacity, but early branches can break this symmetry; no exact full-model symmetry is asserted.",
            "Profiles are bounded local optimizations, not certified global minima or probabilistic intervals.",
            "Zero beta is a diagnostic extension to the physical nonnegative boundary, not a replacement of the original bounded fit.",
            "The zero-versus-1e-7 beta profile RMSE difference is below fixed-parameter N120-to-N240 drift; no physical preference for exactly zero beta or ranking of nearby capacities is resolved.",
            "No terminal-capacity event or mixed-component experiment is supplied by these calibration subsets."])
    (output / "identifiability_results.json").write_text(
        json.dumps(jsonable(result), indent=2) + "\n")
    write_csv(output / "gitis_B1_profile.csv", gprofiles)
    write_csv(output / "alem_capacity_profile.csv", aprofiles)
    write_csv(output / "alem_beta_profile.csv", beta_profiles)
    write_csv(output / "gitis_inactive_thresholds.csv", inactive)
    write_csv(output / "gitis_concentration_scale.csv", concentration_scale)
    write_csv(output / "fixed_parameter_refinement.csv", refinement)
    sensitivity_rows = []
    for dataset, series in [("Gitis log", gsvd), ("Alem alpha/cap log", asvd),
                            ("Alem mixed beta", amixed)]:
        for entry in series:
            for j, singular in enumerate(entry["singular_values"]):
                sensitivity_rows.append(dict(dataset=dataset, step=entry["difference_step"],
                    singular_index=j+1, singular_value=singular,
                    relative_singular_value=entry["relative_singular_values"][j]))
    write_csv(output / "sensitivity_singular_values.csv", sensitivity_rows)
    predictions = [dict(dataset="Gitis", point_id=r["point_id"],
                        observed=float(y), prediction=float(p), residual=float(p-y))
                   for r, y, p in zip(grows, gy, gbase)]
    predictions += [dict(dataset="Alem", point_id=str(i+1),
                        observed=float(y), prediction=float(p), residual=float(p-y))
                    for i, (y,p) in enumerate(zip(ay,abase))]
    write_csv(output / "calibration_predictions.csv", predictions)
    plot_results(result, output)
    print("Completed", result["elapsed_seconds"], "seconds;", calls, flush=True)


def plot_results(result, output):
    plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10,
                         "axes.spines.top":False, "axes.spines.right":False,
                         "savefig.dpi":300})
    fig, axes = plt.subplots(2, 2, figsize=(9, 6.7), layout="constrained")
    ax = axes[0,0]
    for label, entries, color, marker in [
            ("Gitis: 3 log parameters", result["gitis"]["sensitivity"], "#2673a5", "o"),
            ("Alem: log attachment/capacity", result["alem"]["log_alpha_capacity_sensitivity"], "#c95e2c", "s")]:
        for i, entry in enumerate(entries):
            values=np.array(entry["relative_singular_values"])
            ax.semilogy(np.arange(1,len(values)+1), values, marker=marker,
                        color=color, alpha=(.35,.65,1)[i],
                        label=label if i==2 else None)
    ax.set(xticks=[1,2,3], xlabel="Singular-value index", ylabel="Singular value / largest",
           title="(a) Scaled local sensitivity")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, fontsize=9)
    ax.grid(alpha=.18)
    ax=axes[0,1]; p=result["gitis"]["B1_objective_profile"]
    ax.plot([r["B1_over_c0"] for r in p], [r["RMSE"] for r in p], "o-", color="#2673a5")
    ax.axvline(np.exp(result["gitis"]["saved_theta_log"][2]), color=".5",ls=":")
    ax.set(xlabel=r"Fixed charging threshold $B_1/c_0$", ylabel="Calibration RMSE",
           title="(b) Gitis threshold profile")
    ax.grid(alpha=.18)
    ax=axes[1,0]; p=result["alem"]["capacity_objective_profile"]
    ax.plot([r["capacity"] for r in p], [r["RMSE"] for r in p], "s-",color="#c95e2c")
    ax.set(xscale="log", xlabel="Fixed capacity", ylabel="Calibration RMSE",
           title="(c) Alem capacity profile")
    ax.grid(alpha=.18)
    ax=axes[1,1]; p=result["alem"]["one_sided_beta_objective_profile"]
    ax.plot([r["beta_s"] for r in p], [r["RMSE"] for r in p], "o-",color="#657a39")
    ax.set_xscale("symlog", linthresh=1e-6)
    ax.set(xlabel=r"Fixed detachment $\beta_s^*$", ylabel="Calibration RMSE",
           title="(d) Alem detachment profile")
    ax.grid(alpha=.18)
    fig.savefig(output / "identifiability_diagnostics.png")
    fig.savefig(output / "identifiability_diagnostics.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
