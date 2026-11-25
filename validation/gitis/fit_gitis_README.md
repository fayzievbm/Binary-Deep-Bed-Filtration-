# Conditional normalized-shape comparison with Gitis et al 2010

From `validation/gitis/`, run `python ../fit_gitis.py` with the accelerated solver and
the digitized CSV files present. No curve pixels are fitted. The inputs are
the digitized experimental symbols in Figure 6 of Gitis et al., Chemical
Engineering Journal 163 (2010) 78–85, DOI 10.1016/j.cej.2010.07.044.

## Data separation

- Calibration: nine symbols at depth 0.82 m in Figure 6.
- Spatial holdout: eighteen symbols at depths 1.02 and 1.62 m in that same
  experiment. They are never included in the objective function or model-choice
  criterion. These are spatial holdouts, not independent experiments.
- Figure 7: all symbols are evaluated without refitting as an external
  unknown-condition diagnostic. Its run-to-run conditions have not been
  independently matched; the transferred prediction is not called controlled
  external validation.

The `comparison_role` column in every prediction CSV states the role of each
point. Signed residuals are prediction minus observation. Image digitization
envelopes are retained in the CSV and plotted, but they are not experimental
standard deviations and are not used as inverse-variance weights.

## Verified protocol and explicit assumptions

The experimental source describes a 0.44 media porosity. Figure 6 gives approach
velocity 10 m/h, length 1.62 m, and observation depths 0.82, 1.02 and 1.62 m.
The original dimensionless time is `t*=10/1.62 * t_hours`; the factor m0 remains
explicit in the PDE. The scalar specialization uses fixed unit Darcy flow and
the same unchanged balance, fixed-deposit-porosity and hydraulic equations as
the current manuscript.

An experimental inlet mass concentration and particle density sufficient to
determine its volume fraction were not established. Therefore `cin=1e-6` is an
explicit normalized-shape assumption, not a measured concentration. Sensitivity
comparisons use 1e-7 and 1e-5. Constant deposit porosity eps=0.5 is also assumed,
not measured. Absolute retained mass or hydraulic predictions are not validated
against this experiment.

Published time zero contains nonzero particles at all three measurement depths.
The clean-bed model uses an inlet step 0.1 h before the plotted time origin.
This is an explicit clock-offset assumption, not a measured startup time or an
extra fitted parameter. Sensitivity tests use offsets 0.08 and 0.15 h. Main
assumptions are fixed in advance rather than chosen using held-out errors.

## Restricted model fit

The multistage fit estimates only alpha_s, charging multiplier chi, and B1/cin.
Beta_s=0 is an explicit restriction for this comparison. B2/cin=10000 and
Bcap/cin=20000 are fixed far above realized loading. The code checks that aging
and capacity arrest are not reached. Therefore these data assess the charging
to transition change only; they do not validate the aging branch, terminal
capacity, detachment, or mixed-component interactions.

The original linear baseline fits alpha_s and beta_s to the same nine shallow
symbols. All deeper symbols remain unused for both fits. The baseline is not
given different hydraulic or geometric equations.

Unweighted least squares is preceded by a reproducible small differential
evolution search over the three existing multistage kinetic parameters, with
seed 20261125, population multiplier 8 and 45 generations. Both search and
local refinement use shallow symbols only. The search at 80 cells supplies a
starting point for final 160-cell local calibration. Main prediction CSVs,
metrics and plots all use the final 160-cell parameter set. The linear baseline
is also calibrated at 160 cells. No claim of a mathematically proven global
optimum or statistical parameter confidence interval is made.

## Files and interpretation

`fit_gitis_results.json` contains final parameters, calibration and holdout
errors, the coarse-grid initial search, fixed-parameter mesh sensitivity,
assumption sensitivities and external diagnostic errors. Prediction CSVs give
each actual fitted or held-out symbol, the predicted ratio and residual.
The PNG and vector PDF display calibration and spatial holdouts alongside the
same linear baseline. `fit_run.log` records execution output.

Numerical agreement of normalized shapes is conditional on the stated missing
metadata and source-model restrictions. Improvement relative to the linear
baseline does not establish full physical calibration. Model discrepancy,
grid sensitivity, startup uncertainty and transfer failure must be reported
alongside any successful within-experiment spatial prediction.

## Observation mapping

The plotted response R is interpreted as model c/c0 under an assumed proportional measurement response. Related studies report turbidity and oocyst measurements, but the exact 2010 run/observable correspondence is unverified. Values from those related studies have not replaced the declared inlet or timing assumptions. The comparison quantifies normalized-response shape under that mapping, rather than establishing a material-concentration calibration.
