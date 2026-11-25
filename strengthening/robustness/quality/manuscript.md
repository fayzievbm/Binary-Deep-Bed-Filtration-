### Continuous effluent-quality sensitivity

Threshold-defined acceptable effluent volume can change substantially when a numerically spread concentration front lies close to the prescribed threshold. We therefore supplemented the quality-volume comparison with a continuous exceedance measure computed from the existing matched-capacity controls. Let $N$ denote cumulative injected initial pore volumes and define

$$
R(N)=\max_{i=1,2}\frac{c_{i,\mathrm{out}}(N)}{c_{i,\mathrm{in}}},
\qquad
X_\theta=\int_1^{14}[R(N)-\theta]_+\,\mathrm dN.
$$

The thresholds are illustrative, and the first injected pore volume is excluded consistently with the displacement convention used for the quality-volume measure. Lower $X_\theta$ indicates less integrated normalized exceedance on the common 14-PV horizon. It does not define an optimized operating endpoint or a regulatory compliance criterion. Unlike an indicator-defined acceptable-volume measure, this functional is continuous in the outlet trajectory:

$$
|X_\theta(R)-X_\theta(S)|
\leq\int_1^{14}|R(N)-S(N)|\,\mathrm dN.
$$

Each component outlet ratio was interpolated linearly between stored samples separated by 0.02 PV. Threshold crossings and intersections of the component interpolants were inserted before integration. The resulting values are exact for these interpolants. Nodal trapezoidal integration differed by at most $3.55\times10^{-4}$ PV at the declared thresholds; decimation to 0.04 and 0.08 PV changed the values by at most 0.00101 and 0.00565 PV. These checks quantify sensitivity to the stored sampling and do not establish an error bound for the continuum solution.

At 320 cells and $\theta=0.90$, the capped-linear and multistage values were respectively 0.997821 and 1.093268 PV under fixed pressure, and 0.944741 and 1.017231 PV under fixed flow. The multistage exceedance was also greater at thresholds 0.80 and 0.95, with the same ordering on the 80- and 160-cell grids. The largest 160-to-320-cell change across the full 0.5–1 threshold curve was below 0.00360 PV. Thus the continuous comparison supports the conclusion that the multistage closure changes transient effluent quality without providing a universal improvement over a control with the same absorbing capacity. All these results remain illustrative and uncalibrated.

**Figure caption.** Continuous effluent exceedance relative to the capped-linear control, $\Delta X_\theta=X_\theta-X_\theta^{\mathrm{capped\ linear}}$, under fixed pressure and fixed flow. Positive values indicate more exceedance on the common PV 1–14 interval. Lines show the 320-cell results; shaded bands span the 80-, 160- and 320-cell calculations. All controls use the same component capacities and feed. Thresholds are illustrative.
