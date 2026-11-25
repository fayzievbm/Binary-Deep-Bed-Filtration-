# Analytical and controlled numerical supplement

This directory contains the added analysis of the existing multistage closure
within the unchanged 2017 seven-saturation framework. It introduces no additional
physical state, balance equation, evolving deposit porosity, dispersion,
aggregation or external cake. The original seven-case study is retained in
the package's separate `numerics/` directory.

## Files

- `theory.md`: propositions and proofs concerning operating-mode reparameterization,
  a sufficient admissible state region, and conditional interior branch selection.
- `theory_checks.py` and `theory_checks_results.json`: independent algebraic checks,
  including exact rational conservation identities and four uniform steady states.
- `prior_art_comparison.md/.json`: bounded comparison with the closest 2017,
  2022 and 2025 formulations. The screening status and
  manuscript reference numbers have been updated to the completed results; the
  record does not assert worldwide first-ever novelty.
- `numerics/solver.py`: extended computational implementation of the same
  physical framework, with PV tracking, restart/inlet options, safe source
  restrictions and the declared kinetic ablation controls.
- `numerics/reference_solver.py`: unchanged preceding solver used by the
  static-rate consistency check.
- `run_mode_controls.py`, `mode_controls/`: operating-mode equivalence tests,
  static-rate perturbations and capacity-matched controls.
- `run_control_refinement.py`: additional 80/320-cell checks of those controls.
- `run_history.py`, `history/`: eight clean-bed runs demonstrating selection of
  four distinct interior final states by specified inlet histories.

## Reproduction

Python 3, NumPy and Matplotlib are sufficient. `theory_checks.py` uses only the
Python standard library. From this directory run:

```bash
python theory_checks.py
python run_mode_controls.py
python run_control_refinement.py
python run_history.py
```

The mode script produces 20 independent time integrations and two primary
figures. The refinement script adds 12 cases. The history script produces eight
cases and two figures. Machine timings and minor floating-point differences may
vary. Recorded CSV, JSON and NPZ files document the execution reported in the
manuscript. To redraw the mode/control figures without repeating any simulation:

```bash
python -c "import run_mode_controls as s; s.redraw_saved_figures()"
```

The original-reference consistency check can be run separately:

```bash
python -c "import run_mode_controls as s; print(s.reference_check())"
```

All paths resolve inside this directory; no network access or files from the
authoring workspace are needed. The reference comparisons advance four short
static-rate cases rather than repeating the full research sweep.

## Interpretation

The operating-mode theorem concerns continuous solution correspondence under
positive flow, purely hydrodynamic rate coefficients and matching inlet
programs expressed in injected PV. It is not a global existence or uniqueness
theorem. The computations independently advance the original dimensionless time
equations and track injected PV from the transport flux. Their finite-grid
differences are therefore measured and refined rather than forced to zero.
The chosen static-rate perturbations break the equivalence during loading,
although final capacity-filled states coincide.

The capped-linear and charging-without-aging controls use the same absorbing
capacities as the multistage law. They reach identical terminal hydraulic
states. Differences in quality-admissible volume and pump work are transient
effects, and depend on the illustrative quality threshold and operating mode.
The threshold curves are not regulatory water-quality standards. Their volume
omits the first initial PV; small remaining front-related windows can still be
dominated by numerical spreading. Pump work is integrated over the same 14-PV
horizon, including off-specification effluent, rather than optimized stopping.
See `mode_controls/README.md` for definitions, refinement and energy scaling.

The four history-selected interior plateaus arise after different preconditioning
inputs, followed by identical final inlet concentration and pressure. They have
different retained inventories; their hydraulic differences are not benefits at
matched retained mass. Local source root stability is distinguished from full
coupled-PDE stability, and absorbing boundary states remain additional permitted
states. See `history/README.md` for the exact roots and numerical checks.

All calculations are illustrative. The algebraic tests, numerical verification,
physical interpretation and assessment against experimental observations are
separate kinds of evidence. No complete experimental validation is claimed.
