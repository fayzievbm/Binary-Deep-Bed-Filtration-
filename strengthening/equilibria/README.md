# Coupled Darcy-feedback equilibria

This self-contained addition uses the existing model and kinetic parameters to analyze **uniform steady equilibria**, including the original velocity-dependent attachment and pressure-gradient-dependent detachment. It changes no balance, liquid-exchange, permeability or viscosity equation. The computation is algebraic; it does not run a transient filtration solver or claim PDE stability or clean-start reachability.

## Reproduce

Requirements: Python 3.10 or later, NumPy, Matplotlib and mpmath. Tested with Python 3.12.14, NumPy 2.3.5, Matplotlib 3.10.8 and mpmath 1.3.0. Run from any working directory:

```bash
python run_equilibria.py
```

If calling from elsewhere, give the path to the script. By default it writes to the `results/` directory next to the script. Optional flags are `--output-dir`, `--grid-size` (default 321), and `--axis-max` (default 2.0). The default grid spacing is 0.00625. No other manuscript files, web downloads, experimental datasets, or absolute workspace paths are read.

## Files

- `theory.md`: complete assumptions, propositions and proofs; general static/dynamic rates, pure-dynamic specialization, branch-admissibility test, hydraulic ordering and analytic map-boundary formulas.
- `manuscript.md`: compact ready-to-insert theory with four unnumbered display equations, a proposed results paragraph, four-state table and figure caption. Suggested integration: add after the current interior-equilibrium subsection; move the old terminal-equivalence subsection one place later. Add the map/results after the current history-study results. Assign equation/table/figure numbers during manuscript integration.
- `run_equilibria.py`: portable computation and figure generation.
- `results/summary.json`: parameters, provenance, twelve example states and interval-root enclosures, diagnostics and map counts.
- `results/example_equilibria.csv`: all twelve strict example states, original-source residuals, full-saturation checks, hydraulic observables and branch margins.
- `results/pure_dynamic_regime_map.csv`, `mixed_head_regime_map.csv`, `mixed_flow_regime_map.csv`: grid coordinates, number of accepted interior assignments, and a branch bitmask (`LL=1`, `LH=2`, `HL=4`, `HH=8`). Bitmask 15 means all four.
- `results/mixed_mode_difference_example.csv`: individual candidate permeabilities, deposits and signed branch margins at the explicitly identified point where mixed-rate controls admit different branch families.
- `results/equilibrium_regime_map.png` and `.pdf`: publication figure. The PDF keeps lines and text as vectors and rasterizes the sampled color layer. Figure dimensions are 9 by 8.7 inches; a two-column-by-two-row arrangement keeps labels legible when scaled to manuscript width.

## Parameter provenance and interpretation

The retained physical structure and affine coefficient definitions are from Makhmudov, Zh.M.; Saidullaev, U.Zh.; Khuzhayorov, B.Kh. *Fluid Dynamics* **2017**, *52*, 299–308. DOI: [10.1134/S0015462817020132](https://doi.org/10.1134/S0015462817020132).

Porosities `(0.5,0.3)`, charging thresholds `(0.001,0.0002)`, aging thresholds `(0.004,0.001)` and capacities `(0.02,0.004)` are the existing manuscript's illustrative parameters. They are not fitted experimental values. The fixed feed `(0.002,0.0005)` is the existing history study's final feed. The example dynamic coefficients `alpha_di=0.1`, `beta_di=0.05` and forcing magnitude 10 already occur in the existing dynamic-rate calculations. Their combination with this lower feed for the coupled equilibrium calculation is a new analysis, and its transient reachability is not inferred from the earlier static-rate history runs. Mixed panels use original static coefficients `alpha_si=beta_si=0.1`, also present in the existing static-rate sweep.

The map axes `d_i/T_i=alpha_di*gamma_i/(beta_di*mu_gamma*T_i)` vary the original dynamic attachment coefficients at **fixed feed**, fixed viscosity `mu_gamma=1.00625` and fixed `beta_di=0.05`. The first panel has zero static coefficients. The other two use the same mixed coefficients under `P=10` and `u_f=10`. These are not matched final flow/resistance conditions. At this example's axis position, all four equilibria exist in every panel, but their states need not coincide across the mixed controls.

The number 0–4 counts only homogeneous equilibria strictly inside the L/H branches. Zero is **not** absence of every equilibrium: charging floors, thresholds, absorbing caps and nonuniform states are omitted. Branchwise uniqueness does not imply a unique full-model equilibrium. The maps do not identify which branch a time-dependent experiment selects.

Columns named `u_at_P10_same_state` and `pressure_drop_at_u10_same_state` are Darcy observables evaluated at the recorded state. For the purely dynamic case they also give corresponding equilibria under the other control. For mixed rates they do **not** assert that the same recorded state remains an equilibrium when the control changes.

## Checks and numerical precision

All twelve example states pass the strict branch inequalities. The smallest branch margin is `2.3809015277e-4`, the largest absolute original-source residual is `8.6736173799e-19`, and the maximum double-precision permeability error bound is `1.9456330413e-16`.

The map grid uses 60 ordinary floating-point bisection steps; actual final brackets have maximum width `1.1102230246e-16`. This should not be confused with proof of an outward-rounded enclosure for every grid root. Boundaries are drawn from analytic parameterizations independently of the grid.

For the twelve example states only, `high_precision_enclosure` separately constructs finite-decimal brackets of width approximately `2e-45` using 100-digit point arithmetic, then checks strict opposite endpoint signs using **80-digit outward-rounded mpmath interval arithmetic**. These interval evaluations, combined with the proved scalar monotonicity, enclose the unique candidate roots. The endpoint strings, interval sign booleans and upper bounds for the double-precision error are saved in the JSON. This is a software-assisted interval check, not an additional PDE theorem. It does not certify all grid classifications, experimental parameters, or transient dynamics.

On the default grid, the mixed control panels have different admissible branch sets at 11,110 of 103,041 points. These are reproducible counts over a declared rectangle, not probabilities. At `(d1/T1,d2/T2)=(1.4625,0.33125)`, the pressure-controlled case admits `LL,LH,HL,HH`, while the fixed-flow case admits only `HL`. The smallest absolute signed branch margin among the eight candidates there exceeds `6.5e-6`, so this example is not a threshold-equality convention.

## Pure-dynamic four-state summary

| Assignment | Relative permeability | Darcy flux at pressure 10 | Pressure drop at flux 10 |
|:--|--:|--:|--:|
| LL | 0.766671783772 | 7.619098472270 | 13.124912398 |
| LH | 0.749184528678 | 7.445312086239 | 13.431270421 |
| HL | 0.670565722868 | 6.664007183782 | 15.005986225 |
| HH | 0.659429027867 | 6.553331953959 | 15.259413181 |

These states differ in retained inventory. Their hydraulic differences are not a comparison at matched retained mass, an experimental hysteresis validation, or evidence of a universal performance advantage.
