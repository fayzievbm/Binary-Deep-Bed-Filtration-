# Assessment using published observations

The first comparisons yield mixed results. All original 2017 balance and
hydraulic relations remain unchanged; only linear versus multistage deposition
is compared. Both datasets assess an aggregate one-component limit.

| Test and normalized-response or concentration-ratio RMSE | Linear | Multistage |
|---|---:|---:|
| Gitis shallow-depth calibration, 9 symbols | 0.069414 | 0.006710 |
| Gitis spatial holdout, 18 symbols | 0.025765 | 0.013125 |
| Alem constant-flow calibration, 26 symbols | 0.014602 | 0.091086 |
| Alem constant-head holdout, 29 symbols | 0.129193 | 0.142375 |
| Alem constant-flow permeability, 3 unambiguous symbols | 0.121156 | 0.130485 |

Gitis uses an assumed proportional observation mapping from published response R
to model c/c0; its measurement identity/run mapping was not independently verified.
It also uses assumed inlet volume fraction, fixed deposit porosity and startup
clock offset, with explicit sensitivities. Only charging and transition are
visited. The deeper series are spatial holdouts in the same run. Intermediate
depth discrepancies remain, and transfer to Figure 7 fails under incompletely
known run conditions. These results do not identify absolute capacity or
establish hydraulic or binary validation. The multistage fit has three free
kinetic quantities; the linear baseline has two. No complexity-adjusted model
selection or statistical confidence interval is claimed.

Alem uses measured geometry, media porosity, feed concentration and particle
density. Fixed kinetic thresholds, charging multiplier and deposit porosity
are assumptions. Both closures fit attachment and detachment only; the
multistage detachment estimate hits its imposed lower bound. Its capacity is
not reached and its attachment and capacity compensate in sensitivity fits.
The tested multistage specification fits concentrations worse than the linear
closure. This does not prove that all possible multistage specifications fit
worse. Both closures have held-out permeability-ratio RMSE about 0.16–0.17.

The constant-flow linear prediction reproduces the total retained fraction
closely (0.6780 versus reported 0.6795), but gives final effective permeability
0.7207 versus the observed series near 0.94. Reconstructing the nine measured
flow retention sections with assumed deposit porosity 0.5 gives 0.7251. A
separate distribution-independent resistance bound gives maximum permeability
0.8126 at that porosity, or 0.8966 in the compact-deposit limit. This nominal
retention–hydraulic tension is conditional on matching measurement domain,
endpoint and hydraulic normalization. Reported initial conductivity variability
is about 7%; individual ratio errors are unavailable. It is not a formal
statistical rejection. See `alem/results/hydraulic_consistency.json` for the
derivation and `alem/taps_and_profile_check.md` for measurement-domain checks.

Digitization bounds, experimental measurement uncertainty and finite-grid
differences are separate. The 156 digitized symbols comprise 87 Gitis
normalized responses and 69 Alem concentration/hydraulic observations; two aggregate
retained-mass endpoints and nine verified flow section values are transcribed
separately. Source experimental uncertainties have not been invented.

Sources:
- Gitis et al. Chemical Engineering Journal 2010, 163, 78–85.
  https://doi.org/10.1016/j.cej.2010.07.044
- Alem et al. Transport in Porous Media 2015, 106, 303–321.
  https://doi.org/10.1007/s11242-014-0402-8

The manuscript treats these as conditional published-data comparisons. Complete
binary and four-stage validation remains unresolved. No additional laboratory
experiment was conducted and no source-model curve was treated as data.
