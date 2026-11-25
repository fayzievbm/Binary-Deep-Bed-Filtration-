# Hydraulic compatibility analysis

This analysis derives and evaluates a necessary retained-loading/permeability
condition using the original Makhmudov–Saidullaev–Khuzhayorov permeability law.
It uses the existing Alem extraction and changes no model equation, observation
or fitted parameter. The support constant A is derived from the original
permeability function; it is not an added permeability coefficient.

## Reproduce

With Python, NumPy, SciPy and Matplotlib installed, run:

```bash
python analyze_compatibility.py
```

The script is portable: its inputs and outputs are resolved relative to its own
directory. No network access, source-image download, model run or fitting is
required.

## Mathematical result

For admissible occupied fraction D in [0,1), relative local resistance is
r(D)=(1-sqrt(D))^-3. For all D in this interval,

    r(D) >= 1 + A*D,
    A = 21.269750043077014.

For D>0, substitute s=sqrt(D) and minimize

    H(s) = (s*s - 3*s + 3) / (s*(1-s)**3).

The derivative numerator is P(s)=2s^3-8s^2+12s-3; its denominator is positive.
P is strictly increasing on (0,1), with endpoint values -3 and 3, since
P'(s)=6s^2-16s+12 is at least 2 on [0,1]. H diverges at both endpoints, so the
unique root s*=0.308586020218263... is the global minimizer. The support extends
to D=0 by equality. The contact loading is D*=0.095225331874146....

For a homogeneous clean reference permeability and a serial one-dimensional
bed of constant cross section, the effective intrinsic-permeability ratio is
1/mean(r). Integrating the support inequality gives

    k_eff/k0 <= 1/(1 + A*mean(D)).

For binary deposition, mean(D)=sum(mean(B_i)/(1-eps_i)). For one effective
population, mean(D)=mean(B)/(1-eps). Inverting the bound at a nominal ratio K
gives the necessary loading limit

    B_max = (1-eps)*(1/K - 1)/A.

The linear-support coefficient is optimal for the permeability function alone.
For mean(D)<=D*, profiles taking only 0 and D* can attain equality in this
relaxed functional bound. This construction is not a claim that the kinetic
system can generate such a profile. Capacity, measured-profile or reachability
restrictions can tighten the necessary limit.

## Nominal numerical implication

Existing flow retention gives mean(B)=0.005422932692307692. At nominal mapped
ratio 0.94, the bound permits at most 0.0015004827773 when assumed eps=0.5,
which is 27.6692% of the nominal retained loading. The compact-deposit limiting
case eps=0 permits at most 0.0030009655546, or 55.3384%.

These are maximum loadings not excluded by a necessary bound, not estimates of
the true retained loading or a validated correction. The ratio 0.94 is a rounded
level of the published low-concentration constant-flow series, not a newly
identified permeability measurement exactly at the 83-PV mass endpoint.

The two-panel figure shows the loading bounds and deterministic normalization
boundaries for eps=0,0.3,0.5,0.7. Lambda_B multiplies nominal loading and lambda_K
multiplies the nominal mapped ratio. They are hypothetical sensitivity
coordinates, with no fitted values, probabilities or confidence interpretation.
The plotted nominal marker has no invented error bar. Passing a boundary does
not establish agreement with the full profile or the kinetic equations.

## Assumptions and interpretation limits

- Hydraulic and retained-mass records must cover the same bed and endpoint.
- The clean permeability reference must be homogeneous as assumed in this
  derivation; a heterogeneous reference requires the corresponding weighting.
- Conductivity ratios need the appropriate density/viscosity conversion before
  comparison with effective intrinsic permeability. The original dilute feed
  viscosity correction in the selected case is below 0.03%.
- Head normalization and any external apparatus resistance must be matched.
- Initial conductivity is reported with about 7% variability. This is not a
  reported uncertainty interval for K/K0 and is not used to manufacture ratio
  error bars or to select a normalization multiplier.
- Internal deposit porosity remains an assumption, not a measured quantity.
- The result is a nominal deterministic compatibility diagnostic, not a
  statistical rejection, a new constitutive law or complete model validation.

## Files

- `analyze_compatibility.py`: executable derivation checks, numerical tables and
  plotting.
- `inputs/`: copied existing Table 1 retention records, verified flow section
  records and distilled experimental/nominal parameters.
- `input_provenance.json`: original input paths and SHA256 hashes, primary
  source links and comparison boundaries.
- `compatibility_results.json`: constants, precision checks, inversions and
  interpretation limits.
- `compatibility_inversion.csv`: nominal limits for four fixed porosities.
- `compatibility_curves.csv`: original-law bounds versus mean retained loading.
- `normalization_boundaries.csv`: deterministic multiplier boundary curves.
- `normalization_matrix.csv`: sampled scenarios labeled only as not excluded
  or excluded by the necessary bound, not as probabilities.
- `hydraulic_compatibility.png` and `.pdf`: publication comparison figure.
- `manuscript.md`: ready-to-integrate proposition, proof, inversion, table and
  figure caption with unnumbered displays.

## Numerical checks and provenance

The positive root and A are checked independently against an 80-digit Decimal
Newton calculation and a separate scalar minimization. The code checks 60,000
sampled support values including D=0 and D*, inverse/forward consistency, the
section-weighted versus aggregate mean loading, and unchanged input hashes.
These checks support implementation; the analytic derivative argument proves
the pointwise inequality globally.

Original model: Makhmudov et al., Fluid Dynamics 2017, 52, 299–308,
https://doi.org/10.1134/S0015462817020132.

Experimental records: Alem et al., Transport in Porous Media 2015, 106,
303–321, Table 1 and Figure 6,
https://doi.org/10.1007/s11242-014-0402-8.

The experiment source images and full articles are not copied into this
analysis. Existing observations and the original permeability relation are
preserved.
