# Findings from the conditional Gitis comparison

The charging-to-transition closure reproduces the shallow normalized-concentration
decline much better than the original linear closure. Its spatial transfer is
less accurate and retains a clear discrepancy at the intermediate depth. These
results support a limited comparison of normalized concentration shapes under
explicit assumptions; they do not establish full experimental validation.

The final 160-cell fit uses only the nine experimental symbols at 0.82 m.
Differential evolution followed by local least squares estimates
alpha_s*=20.8462, chi=0.164060 and B1/cin=80.5335. Beta_s*=0 is imposed as a
restriction. At the assumed cin=1e-6, B1=8.05335e-5. Fixed B2=0.01 and
Bcap=0.02 are not reached. The largest retained saturation is 0.00153066.
The inferred absolute loading threshold is conditional on the unmeasured inlet
volume fraction and must not be presented as an identified material property.

| Model | Calibration RMSE at 0.82 m | Spatial holdout RMSE at 1.02 and 1.62 m |
|:---|---:|---:|
| Charging-to-transition closure | 0.006710 | 0.013125 |
| Original linear closure | 0.069414 | 0.025765 |

Errors are in units of c/c0. The multistage fit has calibration R2=0.9907 and
combined spatial-holdout R2=0.7808. Holdout RMSE is 0.017415 at 1.02 m and
0.006422 at 1.62 m. At 1.02 m, the largest absolute concentration discrepancy
is 0.039679, substantially larger than the image-digitization envelope. Thus
high calibration R2 does not imply agreement at every observation depth.

The linear baseline estimates alpha_s*=5.01014 and beta_s* approximately zero.
It gives nearly time-independent normalized concentrations after initial
transport, so it cannot reproduce the observed substantial decline using the
same unchanged framework. The comparison gives that baseline its own two fitted
kinetic coefficients and the same nine training symbols.

With the final parameters fixed, changing from 160 to 320 cells changes sampled
normalized predictions by at most 0.003902. The 320-cell calibration and holdout
RMSE values are 0.007267 and 0.013610. Finite-grid uncertainty is therefore
reported separately from image-extraction uncertainty. The final main CSV,
figure and metrics use 160 cells and the same parameter set.

Changing assumed cin from 1e-7 to 1e-5 leaves the normalized predictions and
fitted charging parameters almost unchanged, while all absolute retained
quantities scale with that assumption. This is evidence that the present
normalized data do not identify absolute loading or inlet volume fraction.
Clock offsets of 0.08 and 0.15 h give refitted shallow RMSE values 0.006984 and
0.006708, and heldout RMSE values 0.013183 and 0.013121. The main offset remains
the preselected 0.1 h; it is not selected using spatial holdouts.

Applying the Figure 6 fit to Figure 7 without refitting gives RMSE values
0.2211, 0.1245 and 0.1174 at 0.82, 1.02 and 1.62 m. This transfer does not
agree with those observations. Because complete run-condition equivalence was
not verified, it is reported only as an external unknown-condition diagnostic.
It supports neither a controlled validation claim nor a conclusion that every
parameter should be transferable between the two runs.

No hydraulic or retained-mass observations from Gitis were fitted. Deposit
porosity remains fixed at the assumed value 0.5. The source provides media
porosity 0.44, but these are different physical quantities. No additional
physics, varying deposit porosity, dispersion or new retained state was added.

The maximum species inventory residual for the final multistage simulation is
6.84e-19 in original-pore-volume units, and the minimum carrier saturation is
0.996938. Numerical conservation is strong, but it does not remove experimental
metadata uncertainty or model discrepancy.

## Observation mapping

The plotted response R is interpreted as model c/c0 under an assumed proportional measurement response. Related studies report turbidity and oocyst measurements, but the exact 2010 run/observable correspondence is unverified. Values from those related studies have not replaced the declared inlet or timing assumptions. The comparison quantifies normalized-response shape under that mapping, rather than establishing a material-concentration calibration.
