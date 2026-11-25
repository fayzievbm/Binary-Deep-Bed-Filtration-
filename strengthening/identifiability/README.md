# Conditional parameter identifiability

This supplement diagnoses what the **existing calibration subsets** can inform.
It adds no experiment, measurement, physical mechanism or new fit used in the
article's main comparisons. The original solver, input files, parameter fits
and holdout predictions remain unchanged.

## Reproduction

From the unpacked reproducibility package, run:

```sh
python strengthening/identifiability/run_identifiability.py
```

The script discovers the sibling `validation/` directory. A source checkout's
`work_validation/` is also supported. An explicit location can be supplied with
`--validation-root PATH`; `--output PATH` selects the output directory. No
author-specific workspace path is embedded. Requirements: Python 3, NumPy,
SciPy, Matplotlib, and a C++17 compiler (`g++`). The unchanged scalar kernel is
compiled into this diagnostic's `kernel_cache/`, keeping writes away from the
original validation files. That cache may be deleted and regenerated.

## Data and protocol

- Gitis: only the nine Figure 6 symbols at 0.82 m; 160 cells, assumed
  `c0=1e-6`, deposit porosity 0.5, and the existing 0.1 h clock offset.
- Alem: only the 26 constant-flow concentration symbols; 120 cells, measured
  feed conversion and existing prescribed early thresholds/charging.
- No deeper Gitis observations, Gitis Figure 7 responses, Alem constant-head
  concentrations, permeability symbols or retained-mass values enter the
  sensitivities or profile objectives. There is no holdout-driven selection.
- The starting points are the archived multistage fits. These are not replaced
  by diagnostic reoptimizations. Comparisons with original baseline errors
  always use the **same** mesh as each screen.

The sensitivity matrix is the derivative of the normalized model response,
divided by the RMS calibration response, with respect to natural-log parameter
coordinates. Three central-difference steps (0.02, 0.01, 0.005) check local step
sensitivity. Gitis uses attachment, charging and the normalized charging
threshold; Alem first uses attachment and capacity with detachment fixed.
The tested differences **do not demonstrate step convergence** of the Jacobian:
Gitis's second singular value changes from 1.184 through 2.225 to 3.957, and
Alem's smaller attachment/capacity singular value from 0.03845 through 0.04139
to 0.11027. These are numerical finite-step screens, not converged derivative
estimates or evidence for a precise statistical rank.

For Alem, an additional mixed-coordinate screen includes the scaled additive
detachment coordinate `beta_s * PV_end`. Because the archived detachment lies
near its nonnegative boundary, a forward difference is used. This avoids
mistaking the smallness of a log-detachment derivative near zero for evidence
that any detachment value is admissible. Singular values depend on the declared
coordinates, scales, differencing and mesh; no universal rank threshold is used.

Selected bounded objective profiles fix Gitis `B1/c0` and reoptimize attachment
and charging from two starting points, fix Alem capacity and reoptimize attachment and detachment, or fix
Alem detachment (including zero) and reoptimize attachment. These are conditional
local searches, not certified global minima. Their objective is unweighted
squared residuals because experimental variances are unavailable. They are
**not** profile-likelihood confidence intervals or significance tests.

## Outputs

- `identifiability_results.json`: complete protocol, matrices, singular values,
  directions, objective profiles, stage events, baseline/refinement checks,
  input SHA-256 hashes, call counts and elapsed time.
- `sensitivity_singular_values.csv`: all step-dependent singular values.
- `gitis_B1_profile.csv`, `alem_capacity_profile.csv`, `alem_beta_profile.csv`:
  conditional profile fits and objective ratios relative to the saved fit.
- `gitis_inactive_thresholds.csv`: independently varied, unvisited aging and
  capacity thresholds and their response changes.
- `gitis_concentration_scale.csv`: absolute concentration and all deposit
  thresholds rescaled together, with active normalized parameters fixed.
- `calibration_predictions.csv`: the archived parameter predictions on the
  calibration subsets only.
- `fixed_parameter_refinement.csv`: fixed-parameter 160/320-cell Gitis and
  120/240-cell Alem checks; no additional calibration on finer meshes.
- `identifiability_diagnostics.png` and `.pdf`: four-panel diagnostic figure.
- `manuscript.md`: ready-to-insert methods/results text, caption and a compact
  stage-activation/estimability table.
- `stage_activation_estimability.csv`: the same two-row scientific table in
  reusable form.
- `verify_identifiability.py` and `verification_checks.json`: archived-fit
  reproduction, original-input integrity, inactive-stage, admissibility, profile
  convergence and portable-directory-discovery checks.

## Interpretation limits

An unvisited threshold is uninformative only over parameter perturbations that
keep it unvisited. Gitis cannot establish aging or capacity arrest. Rescaling
`c0` and deposit thresholds approximately preserves normalized concentration
responses, but the unchanged mobile-pore-volume and hydraulic relations depend
on absolute loading; no exact scale invariance is asserted.

In Alem's aging branch the source contains the product `alpha_s * capacity`.
This can produce a shallow objective ridge, while the early charging/transition
branches and threshold crossings can break the product symmetry. No exact
global nonidentifiability or unique material capacity is inferred. A fitted
detachment at a search boundary does not establish physical irreversibility.
The optimized RMSE difference between zero and `1e-7` detachment is about
`5.9e-7`, below the `5.26e-4` fixed-parameter drift from 120 to 240 cells. It
cannot resolve a physical preference for exactly zero versus small positive
detachment. Small objective differences between nearby capacities likewise do
not establish a ranking or a material capacity estimate.

The original conditional observation mapping, unknown inputs and hydraulic
mismatch remain limitations. None of these diagnostic calculations validates
binary interactions, all four kinetic stages or physical hysteresis.
