# Continuous effluent exceedance

This folder postprocesses the existing matched-capacity numerical controls. It introduces no physical equation, fitted parameter, or new PDE simulation. All 18 required input CSVs are bundled in `inputs/`: three kinetic laws (capped linear, charging followed by linear, multistage), two operating modes, and 80/160/320 cells. The 160-cell source filenames have no grid suffix. The inlet volume fractions are (0.05, 0.01), the endpoint is 14 injected initial pore volumes (PV), and the first PV is omitted from this metric.

Run from any directory:

```bash
python /path/to/quality/analyze_quality.py
```

Dependencies are NumPy and Matplotlib. The script writes `quality_results.json`, `quality_metrics.csv`, `quality_curves.csv`, and `figure_continuous_quality.png`/`.pdf`. Input SHA256 hashes are included in the JSON. No inputs outside this folder are needed.

## Metric and numerical treatment

With $R(N)=\max_i[c_{i,\mathrm{out}}(N)/c_{i,\mathrm{in}}]$, the reported measure is

$$
X_\theta=\int_1^{14}[R(N)-\theta]_+\,dN.
$$

It measures the magnitude and duration of normalized effluent exceedance above an illustrative threshold. Lower values mean less exceedance on this prescribed common horizon. The thresholds are not drinking-water standards, and the horizon does not optimize a stopping policy. This metric uses the maximum of the two component ratios, so it is not an absolute pollutant-mass discharge measure.

Each component ratio is interpolated linearly between the stored samples (spacing 0.02 PV). The integration inserts threshold crossings and intersections of the two component lines, then integrates exactly on each resulting segment. Thus it is exact for the stated saved-data interpolant, not for the unsaved continuous solver trajectory or the continuum solution. Two analytical crossing examples are checked in the script. Across the selected thresholds 0.8, 0.9 and 0.95, the largest difference between this calculation and nodal trapezoidal integration is $3.55\times10^{-4}$ PV. Decimation to 0.04 and 0.08 PV changes the result by at most 0.00101 and 0.00565 PV, respectively. These are sampling-sensitivity estimates, not certified sampling-error bounds.

The metric satisfies

$$
|X_\theta(R)-X_\theta(S)|\leq\int_1^{14}|R-S|\,dN.
$$

The script verifies a conservative version of this inequality between successive saved-grid interpolants, using the sum of componentwise integrated absolute ratio differences. These are bounds between numerical interpolants; they do not bound the unknown continuum error. The largest 160-to-320-cell change in $X_\theta$ across the full 0.5–1 threshold curve is below 0.00360 PV. The CSV also includes `sampled_interpolant_good_PV`; this differs conceptually from the threshold counter accumulated at internal solver steps.

## Findings

The 320-cell results are:

| Mode | Threshold | Capped linear $X_\theta$ | Multistage $X_\theta$ | Multistage minus capped linear |
|---|---:|---:|---:|---:|
| Pressure | 0.80 | 2.241286 | 2.294345 | 0.053059 |
| Pressure | 0.90 | 0.997821 | 1.093268 | 0.095447 |
| Pressure | 0.95 | 0.498050 | 0.539624 | 0.041573 |
| Flow | 0.80 | 2.240738 | 2.282324 | 0.041586 |
| Flow | 0.90 | 0.944741 | 1.017231 | 0.072490 |
| Flow | 0.95 | 0.383893 | 0.466656 | 0.082763 |

The multistage law has greater exceedance than the capped-linear control at all six declared mode/threshold combinations, and the ordering is the same at 80 and 160 cells. This result complements the threshold-defined quality-volume comparison: it does not support a universal performance advantage from multistage kinetics. The small charging-versus-capped-linear differences should be interpreted with the reported grid and sampling sensitivities rather than promoted as a new optimization result.

The figure shows differences relative to the capped-linear control, with 320-cell lines and bands spanning 80/160/320 cells. The complete absolute values and threshold curves are retained in the CSV/JSON outputs. These simulations remain illustrative and uncalibrated.
