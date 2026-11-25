# Additional analytical and numerical evidence

All four studies preserve the original 2017 balances, fixed deposit porosities,
liquid exchange, concentrations, Darcy law, permeability and viscosity. Only
kinetic closures are compared. The original sharp law remains the primary
model; aging-transition smoothing is an explicitly declared sensitivity variant.

## Studies

- `equilibria/`: general static/dynamic rate proposition, pure-dynamic strict
  admissibility test, hydraulic ordering, three 321-by-321 regime maps and twelve
  example roots independently checked by outward-rounded interval arithmetic.
  Interior branchwise uniqueness does not imply a unique or stable PDE state.
- `robustness/`: 56 history runs plus two refinements, conservative smoothing
  source map and independent source verification. Four plateaus persist at the
  tested relative half-widths through 0.4; all four tested histories approach one
  interior plateau at 0.6. Local merger widths 0.521859/0.529981 concern this static
  final-feed smoothing family. Absorbing caps remain additional states.
- `robustness/quality/`: continuous exceedance computed from 18 saved original
  control histories, with piecewise-linear crossing integration and sampling
  sensitivity. No new PDE solution or experimental standard is introduced.
- `identifiability/`: only the existing 9/26 calibration observations enter
  finite-step sensitivity screens and local residual-objective profiles.
  Original fits and holdouts stay unchanged. Jacobian step convergence is not
  established; no confidence interval, global identifiability or resolved
  physical preference for zero detachment is claimed.
- `compatibility/`: a global support inequality for the unchanged permeability
  yields necessary mass–hydraulic limits and deterministic normalization curves.
  The rounded hydraulic ratio and recovered mass are a nominal mapped pair,
  without verified endpoint pairing or a statistical rejection interpretation.

## Reproduce from the package root

```bash
python strengthening/equilibria/run_equilibria.py
python strengthening/robustness/source_checks.py
python strengthening/robustness/run_smoothing.py
python strengthening/robustness/quality/analyze_quality.py
python strengthening/identifiability/run_identifiability.py
python strengthening/compatibility/analyze_compatibility.py
```

Full history runs and diagnostic profile searches take several minutes; all
reported CSV, JSON and NPZ outputs are already supplied. To reuse the saved
smoothing runs and regenerate their summary/figures, add `--reuse` to that
script. Consult each directory README for options and precision limits.
Python3, NumPy, SciPy, Matplotlib and mpmath are required across the studies;
the identifiability scalar kernel compiles with a C++17 compiler on first use.
Compiled caches and third-party full-text/source-image files are excluded.

The new main-article sections are 3.13–3.15,4.4/4.6,5.6–5.7 and6.3–6.4;
Section 5.4 and Figure S4 add the continuous quality metric. The other recorded
numerical results, original source/kernel files and calibration holdouts were
preserved. All 58 display bodies of the preceding manuscript remain present.

These additions strengthen model analysis and quantify applicability limits.
They do not supply missing binary experimental observations or establish
physical hysteresis.
