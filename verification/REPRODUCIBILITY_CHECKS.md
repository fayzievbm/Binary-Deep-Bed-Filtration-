# Reproducibility checks

All **12 commands passed**. Commands below were executed from the repository root against the supplied scripts and saved observations/results. All **28 Python source files** also passed syntax parsing.

| Command | Check performed | Result | Elapsed (s) |
|---|---|---|---:|
| `python analysis/theory_checks.py` | algebraic assertions | PASS | 0.035 |
| `python validation/benchmark_validation_solver.py` | numerical reference consistency | PASS | 1.053 |
| `python strengthening/identifiability/verify_identifiability.py` | saved-fit and input-integrity checks | PASS | 0.588 |
| `python validation/plot_gitis.py` | saved-data figure regeneration | PASS | 0.626 |
| `python validation/plot_alem.py` | saved-data figure regeneration | PASS | 2.077 |
| `python validation/check_alem_overlap.py` | saved-prediction subset statistics | PASS | 0.022 |
| `python validation/hydraulic_consistency.py` | retained-loading/permeability calculation | PASS | 0.318 |
| `python strengthening/compatibility/analyze_compatibility.py` | analytical bound, precision checks and figure | PASS | 0.976 |
| `python strengthening/robustness/quality/analyze_quality.py` | saved-history continuous-quality calculation | PASS | 11.050 |
| `python strengthening/robustness/run_smoothing.py --reuse` | saved-run summary and figure regeneration | PASS | 1.449 |
| `python strengthening/identifiability/run_identifiability.py --help` | command-line parsing only | PASS | 0.520 |
| `python strengthening/equilibria/run_equilibria.py --help` | command-line parsing only | PASS | 0.377 |

## Numerical benchmark

The scalar/reference comparison covered five kinetic and operating-mode cases. Recorded maximum absolute differences were:

- Final mobile/deposited state: **4.85723e-17**.
- Flow or pressure: **2.13163e-14**.
- Conservation balance residual: **1.02696e-15**.
- Pore-volume/time conversion: **7.88702e-12**.

## Scientific data integrity

SHA-256 comparisons confirmed that **376 CSV, JSON and NPZ scientific input/result files** are byte-identical to their archived execution records. Counts by directory: `analysis/`: 86, `numerics/`: 16, `strengthening/`: 242, `validation/`: 32.

The companion `REPRODUCIBILITY_CHECKS.json` records each compared repository-relative filename and hash, the comparison-set digest, the commands, pass types and benchmark values. The original scientific data were not altered by these checks.

## Tested environment and scope

Python 3.12.14; NumPy 2.3.5; SciPy 1.17.0; Matplotlib 3.10.8; mpmath 1.3.0; Pillow 13.0.0; g++ with C++17 support.

The full PDE scenario sweeps and nonlinear calibration searches were not rerun in this audit; their archived results are supplied. The two `--help` invocations test command-line parsing only. Timings describe this execution and are not performance guarantees. These numerical and integrity checks do not add experimental validation.
