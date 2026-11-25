# Binary Deep Bed Filtration with Multistage Deposition Kinetics

Research code, digitized published observations, numerical results, and figures for a two-component deep-bed filtration study under fixed pressure or fixed flow.

The governing framework is the seven-saturation model of Makhmudov, Saidullaev, and Khuzhayorov (2017). **Only the net deposition kinetics are extended.** The original saturation balances, fixed deposit porosities, immobilized-liquid relations, concentration definitions, Darcy equation, permeability law, and suspension viscosity are retained. The studies examine operating-mode correspondence, stage-dependent retention, admissibility, equilibrium selection, parameter information, and limits of agreement with existing published observations.

No new laboratory experiments are reported. The published-data comparisons concern a conditional one-component specialization; they do not establish binary, complete four-stage, or hydraulic validation.

## Model and numerical conventions

The primary kinetic closure has four loading regimes for each particle species: charging, reversible transition, loading-dependent aging, and absorbing terminal capacity. It follows the phenomenological form in Fayziev (2020), Equation 8, with independent species thresholds and capacities. Aging changes a kinetic multiplier; deposit porosity does not evolve and no maturation state is added.

The implementation advances mobile particle saturation `a_i=S_2i` and deposited particle saturation `B_i=S_4i`. The remaining original saturations are reconstructed from the unchanged fixed-porosity identities. Time retains the original scale `t*=U0*t/L`, with media porosity `m0` explicit:

```text
m0 * d(a_i)/dt* + d(u* c_i)/dx* = -q_i
m0 * d(B_i)/dt*                  =  q_i
```

The original constitutive laws are `k/k0=(1-sqrt(D))^3`, `mu/mu0=1+2.5*(c1+c2)`, and `D=sum(B_i/(1-eps_i))`. Original velocity-dependent attachment and pressure-gradient-dependent detachment remain available. Positive flow is assumed in the supplied scenarios.

The sharp-law solver uses conservative first-order upwind transport and Lie splitting. Frozen-concentration source branches are integrated analytically with threshold events; each deposit increment has an equal and opposite mobile increment. The first charging threshold uses the declared attractive-interface holding completion. Capacity arrest sets the complete net source to zero, including detachment. Exact local branch integration does not make the coupled method exact or second order. The smoothing study is a separately declared kinetic sensitivity variant.

Injected pore volumes are integrated as `PV=integral(u*/m0 dt*)`; they are distinct from `t*`. Inventories are normalized by original pore volume. Numerical verification, analytical results, and comparison with published observations are separate forms of evidence.

## Repository contents

| Directory or file | Contents |
|---|---|
| [`numerics/`](numerics/README.md) | Baseline binary solver; seven illustrative scenarios; local-event, conservation, spatial/time refinement, and species-concentration overshoot checks; CSV histories, NPZ fields, JSON metrics, and plots. |
| [`analysis/`](analysis/README.md) | Theory and algebraic checks; a separate solver implementation for injected volume/work diagnostics and kinetic controls; 32 operating-mode/control calculations and eight feed-history calculations. |
| [`strengthening/equilibria/`](strengthening/equilibria/README.md) | Coupled uniform equilibria, three 321-by-321 regime maps, twelve example states, and interval checks for those example roots. |
| [`strengthening/robustness/`](strengthening/robustness/README.md) | Independent source checks and 58 history/refinement calculations testing finite smoothing of the aging transition. |
| [`strengthening/robustness/quality/`](strengthening/robustness/quality/README.md) | Continuous effluent-exceedance postprocessing of 18 supplied matched-capacity histories, including grid and sampling sensitivity. |
| [`strengthening/identifiability/`](strengthening/identifiability/README.md) | Calibration-subset sensitivity screens, conditional objective profiles, inactive-stage checks, and fixed-parameter refinement. |
| [`strengthening/compatibility/`](strengthening/compatibility/README.md) | Necessary retained-loading/permeability bound and nominal hydraulic compatibility diagrams using the unchanged permeability law. |
| [`validation/`](validation/ASSESSMENT.md) | Digitized published symbols, metadata, extraction provenance, restricted fits, holdout predictions, scalar Python/C++ solver, plots, and consistency checks. |
| `figures/` | Assembled article and supplemental figures in PNG and vector PDF formats. Study scripts regenerate figures within their own directories; assembled copies are not automatically synchronized. |
| `literature/` | Reference catalog, primary-source verification notes, and literature screening records. |
| [`verification/REPRODUCIBILITY_CHECKS.md`](verification/REPRODUCIBILITY_CHECKS.md) | Executed checks and independent preservation audit of the recorded scientific data. |
| `SHA256SUMS.json` | SHA-256 digests of distributed files, excluding this manifest itself. |
| `requirements.txt` | Pinned Python dependency versions from the recorded computational environment. |

The `numerics/` baseline and reference copies are kept separately from the extended analysis and smoothing implementations. Their purpose and numerical differences are documented in the linked READMEs. The source distribution includes C++ source rather than a machine-specific compiled library.

## Installation

The recorded environment is **Python 3.12.14**, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8, mpmath 1.3.0, and Pillow 13.0.0. A Linux environment with Python 3.12 is the closest match to the supplied execution. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Published-data fitting, the scalar benchmark, and parameter-information calculations also need **`g++` with C++17 support**, available on `PATH`. `validation/validation_solver.py` builds the supplied `validation_kernel.cpp` on first use with `-O3 -std=c++17 -shared -fPIC`. The identifiability study uses its own disposable kernel cache. The compiler is a system dependency and is not installed by `pip`.

Pillow is used only to recreate extraction overlays from source images. The supplied CSV observations suffice for fitting, plotting, and all scientific analyses. Pandoc and office software are not required to reproduce the code/results.

## Quick checks and saved-result reproduction

Run the following from the repository root. They avoid full filtration sweeps and optimizer searches:

```bash
# Verify the distributed record before regenerating outputs.
python tools/verify_checksums.py

# Algebraic identities and analytical states; Python standard library only.
python analysis/theory_checks.py

# Independent frozen-source integration checks; no full filtration histories.
python strengthening/robustness/source_checks.py

# Saved-history quality metrics and necessary hydraulic compatibility analysis.
python strengthening/robustness/quality/analyze_quality.py
python strengthening/compatibility/analyze_compatibility.py

# Published-data diagnostics and plots using supplied observations/predictions.
python validation/hydraulic_consistency.py
python validation/check_alem_overlap.py
python validation/plot_gitis.py
python validation/plot_alem.py

# Integrity and saved-record checks for the parameter-information study.
python strengthening/identifiability/verify_identifiability.py
```

The following short benchmark performs numerical integrations and compiles the scalar kernel on first use, but does not fit any observations:

```bash
python validation/benchmark_validation_solver.py
```

To redraw the operating-mode/control figures and reaggregate complete saved smoothing runs without rerunning those histories:

```bash
(cd analysis && python -c "import run_mode_controls as s; s.redraw_saved_figures()")
python strengthening/robustness/run_smoothing.py --reuse
```

The `--reuse` option reruns an individual case if its saved record is missing. Most commands regenerate outputs in place. Use a fresh checkout when comparing regenerated files with the distributed record. The saved-record integrity check should be run before refitting, because it checks hashes of the original calibration inputs/results.

## Full reproduction

All paths below are relative to the repository root. Each study can be run separately. Complete sweeps and restricted calibration/profile searches may take several minutes or longer depending on hardware. Recorded CSV, JSON, and NPZ outputs are already supplied.

| Study | Command | Main saved record |
|---|---|---|
| Baseline kinetics, profiles, sensitivities, and verification | `python numerics/solver.py` | `numerics/summary.json`, `case_metrics.csv`, per-case CSV/NPZ files |
| Operating-mode correspondence, static-rate perturbations, matched-capacity controls | `python analysis/run_mode_controls.py` | `analysis/mode_controls/results.json` |
| Additional control refinement | `python analysis/run_control_refinement.py` | `analysis/mode_controls/control_refinement.json` |
| Four feed-history-selected interior states on two meshes | `python analysis/run_history.py` | `analysis/history/history_results.json` |
| Coupled Darcy-feedback equilibrium maps | `python strengthening/equilibria/run_equilibria.py` | `strengthening/equilibria/results/summary.json` and regime-map CSVs |
| Aging-transition smoothing sensitivity and refinements | `python strengthening/robustness/run_smoothing.py` | `strengthening/robustness/results/smoothing_results.json` |
| Restricted Gitis calibration and held-out predictions | `python validation/fit_gitis.py` | `validation/gitis/fit_gitis_results.json` |
| Restricted Alem calibration and held-out predictions | `python validation/fit_alem.py` | `validation/alem/results/fit_results.json` |
| Calibration-subset parameter-information diagnostics | `python strengthening/identifiability/run_identifiability.py` | `strengthening/identifiability/identifiability_results.json` |

After regenerating fits or histories, rerun the corresponding diagnostic/plot commands from the quick-check section. The quality analysis has its own bundled input trajectories; its README records their provenance. Optional script arguments, parameter definitions, output schemas, and numerical precision limits are described in each study's README. Optimizer and floating-point differences across environments can change the last digits; the supplied files record the reported execution.

## Interpretation and applicability limits

- The binary simulations use illustrative parameters. Operating-mode correspondence has stated positive-flow and rate-coefficient assumptions. The invariant-region result is sufficient, and local source/equilibrium analysis does not prove general coupled-PDE stability or uniqueness.
- Matched-capacity controls share identical filled-state hydraulics. Transient effluent quality and work depend on operating mode and the illustrative concentration threshold; no universal multistage performance advantage is established. Very small front-related quality windows can be dominated by numerical spreading.
- Different preconditioning histories inject and retain different particle inventories. Their distinct plateau fluxes are not hydraulic benefits at matched retained mass. Smoothing sensitivity supports a closure-specific selection result, not experimental physical hysteresis.
- **Gitis et al. (2010):** nine shallow-depth symbols calibrate a restricted charging/transition model; eighteen deeper symbols are spatial holdouts in the same run. The observation mapping, absolute inlet scale, deposit porosity, and startup clock include explicit assumptions. Aging and capacity arrest are not visited. Transfer to Figure 7 is an unknown-condition diagnostic and does not establish controlled external validation or validated hydraulics.
- **Alem et al. (2015):** 26 constant-flow concentration symbols calibrate attachment and detachment; constant-head concentrations, hydraulic observations, and retained mass are withheld. The tested multistage specification fits concentrations worse than the linear baseline. Hydraulic discrepancies remain for both closures. Polydisperse aggregate observations do not validate two resolved particle species.
- Parameter-information screens use calibration observations only. Jacobian step convergence is not demonstrated; conditional objective profiles are not confidence intervals or global identifiability proofs. Boundary detachment fits do not establish physical irreversibility.
- The hydraulic compatibility diagrams are necessary deterministic bounds, conditional on a common measurement domain, endpoint, and normalization. The nominal published loading/ratio pair lacks independently verified endpoint matching; the comparison is not a statistical rejection test.

See [`validation/ASSESSMENT.md`](validation/ASSESSMENT.md) for the reported comparison metrics and [`validation/validation_protocol.md`](validation/validation_protocol.md) for dimensional mappings and calibration/holdout discipline.

## Published-data provenance

The distributed observations were digitized from experimental symbols or transcribed from published tables, not sampled from source-model curves. Metadata retain source links, image hashes, marker coordinates, axis calibrations, and extraction notes. Pixel-reading envelopes are not experimental standard deviations. Original laboratory spreadsheets, complete third-party papers, and source figure rasters are not included. The digitization scripts are supplied for audit; running them requires the identified original images at the paths documented by the scripts and metadata.

## Model and data sources

- Makhmudov, Zh.M.; Saidullaev, U.Zh.; Khuzhayorov, B.Kh. *Mathematical model of deep-bed filtration of a two-component suspension through a porous medium.* Fluid Dynamics **2017**, 52, 299–308. [DOI: 10.1134/S0015462817020132](https://doi.org/10.1134/S0015462817020132).
- Fayziev, B. *A phenomenological model of suspension filtration in porous medium.* International Journal of Applied Mathematics **2020**, 33, 511–521. [DOI: 10.12732/ijam.v33i3.10](https://doi.org/10.12732/ijam.v33i3.10).
- Gitis, V.; Rubinstein, I.; Livshits, M.; Ziskind, G. *Deep-bed filtration model with multistage deposition kinetics.* Chemical Engineering Journal **2010**, 163, 78–85. [DOI: 10.1016/j.cej.2010.07.044](https://doi.org/10.1016/j.cej.2010.07.044).
- Alem, A.; Ahfir, N.-D.; Elkawafi, A.; Wang, H.Q. *Hydraulic Operating Conditions and Particle Concentration Effects on Physical Clogging of a Porous Medium.* Transport in Porous Media **2015**, 106, 303–321. [DOI: 10.1007/s11242-014-0402-8](https://doi.org/10.1007/s11242-014-0402-8).

The literature records also discuss the close two-component multistage (2022) and charging (2025) formulations. The contribution is framed around derived properties and controlled comparisons; no first-ever multistage-filtration claim is made.
