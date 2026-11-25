# Validation protocol for the unchanged 2017 framework

Prepared 25 November 2026. This note specifies dimensional mapping, calibration and held-out tests. It does not assert that fitting has succeeded or that a complete binary experiment has been validated. Only the net deposition closure may change; all seven original balances, fixed deposit porosity, liquid exchange, concentration definitions, Darcy law, permeability and viscosity remain unchanged.

## 1 Decisions for the current experimental work

1. Alem et al. (2015) provides the primary aggregate concentration-and-hydraulic test. Begin with a one-component effective representation of its polydisperse kaolinite suspension. Set c2=S22=S32=S42=0. This is a valid limiting case of the binary equations, but cannot establish species-resolved binary validity.
2. Calibrate only on the 0.25 g/L constant-flow concentration series. Keep the 0.25 g/L constant-head concentration series, hydraulic histories in both modes, and final retained-depth profiles out of the objective. A successful concentration fit followed by hydraulic mismatch is a meaningful failure of the fixed framework, not permission to change k or liquid exchange.
3. Gitis et al. (2010) is a conditional normalized-shape and spatial-generalization test unless actual feed concentration, particle density and needed hydraulic metadata are recovered from its underlying experimental sources. Its laboratory curves do not reach terminal breakthrough, so they cannot validate capacity arrest or all four branches.
4. Do not claim that the paper's illustrative deposit porosity 0.5 is an experimental measurement. It may be fixed as a declared modeling assumption, with sensitivity assessed separately. Do not tune it on the hydraulic data and then describe those same data as held out.

## 2 Experimental quantities and their model counterparts

For one particle material of density rho_p (kg/m3), with suspended mass concentration C_m (kg/m3 of suspension), the mobile particle volume fraction is c=C_m/rho_p. This is exact when C_m is defined per suspension volume and is a dilute approximation if the laboratory concentration uses carrier-water volume. Numerical conversions are 1 mg/L=0.001 kg/m3 and 1 g/L=1 kg/m3.

The mobile saturation is S21=a=phi*c, not c. With B=S41 and fixed deposit porosity eps,

    phi=1-B/(1-eps), S31=eps*B/(1-eps), S1=phi-a.

A retained mass concentration sigma_m in kg/m3 of BULK bed maps to B=sigma_m/(m0*rho_p). A retained solid-volume fraction sigma_v per BULK bed maps to B=sigma_v/m0. Alem defines its specific deposit as the latter. Neither quantity equals B without the m0 conversion.

For cross-section A and bed length L,

    M_retained=rho_p*m0*A*L*mean_xstar(B),
    M_mobile=rho_p*m0*A*L*mean_xstar(a).

The exact species balance is

    M_in-M_out = M_retained+M_mobile

for initially clean pores. A deposit estimate obtained by integrating inlet-minus-outlet flux includes the mobile inventory unless that correction is made. Thus Alem's BTC-derived mean-deposit curve is not independent validation of the same BTC. The end-of-run washed/excavated retained-mass profile is the stronger independent target, subject to the reported experimental mass-recovery error.

If two materials are ever measured separately, total mass concentration is rho_p1*c1+rho_p2*c2. A turbidity signal or a particle-number distribution cannot be substituted for this sum without its measurement response and density conversion.

## 3 Time, pore volumes and hydraulic forcing

Keep the original dimensionless variables:

    x*=x/L, t*=U0*t/L, U*=U/U0,
    p*=k0*p/(U0*mu0*L), alpha*=L*alpha/U0, beta*=L*beta/U0.

The factor m0 remains in every saturation time derivative. With the INITIAL bed pore volume Vp0=m0*A*L, cumulative injected pore volumes are

    NPV(t)=integral_0^t Q(s) ds / Vp0
          =(1/m0)*integral_0^tstar U*(s) ds.

Hence NPV=t*/m0 only when U*=1. Under constant head, integrate NPV alongside the model; do not map the experimental NPV axis to time using the initial flow rate throughout. Compare C/C0 at the model's own cumulative NPV. If measured Q(t) is instead imposed, the test is a flow-conditioned transport comparison and cannot simultaneously count as validation of predicted Darcy flow.

Measured flow Q maps to Darcy velocity U=Q/A. Pore velocity is U/(m0*phi); its normalized value is (U/U0)/phi. Do not enter a reported pore velocity directly as U. If a published pore-velocity curve was derived from concentration-based porosity rather than directly measured, it is not an independent hydraulic target.

Hydraulic conductivity K (m/s) differs from intrinsic permeability k (m2): k=K*mu/(rho_w*g). For a vertical column, compare the Darcy driving potential corresponding to measured total hydraulic-head drop, DeltaP_eff=rho_w*g*DeltaH, rather than blindly using a raw pressure-head difference without elevation. This is an experimental-to-potential conversion, not an added filtration mechanism.

The model's predicted normalized column hydraulic conductance is

    G/G0 = [integral_0^1 (mu/mu0)/(k/k0) dx*]^-1.

It equals U/Uclean at fixed head and DeltaP_clean/DeltaP at fixed U. Only for uniform viscosity can one identify it directly with the harmonic effective intrinsic-permeability ratio.

## 4 Alem 2015: verified metadata and unresolved quantities

Primary record: https://link.springer.com/article/10.1007/s11242-014-0402-8 . Author-deposited primary article: https://www.academia.edu/52371161/Alem_et_al_2014_Transp_Porous_Med . The accepted/online year is 2014; the volume106 issue citation is 2015.

Verified experimental anchors: L=0.40m, inner diameter0.045m, m0=0.37 +/-0.01, particle density2600kg/m3, constant flow145mL/min, sand conductivity about7e-4m/s, concentration series0.25/0.50/1.00/1.50g/L, and duration83 initial pore volumes. The lowest feed maps to c0=9.61538e-5. Diameter/flow give A=0.00159043m2 and U0=0.00151953m/s, consistent with the paper's stated0.152cm/s. The initial pore-volume travel time is about97.4s. Fixed-flow83PV corresponds to t*=30.71, approximately2.25h.

The column is vertical. Both modes begin with approximately the same measured flow. A constant-head reservoir setting of105cm is stated. However, combining K0=7e-4m/s, DeltaH=1.05m and L=.40m predicts U0=.0018375m/s, about21% above the flow inferred from145mL/min. Do not silently impose all these rounded/global values as exactly consistent bed-scale measurements. Prefer the initial pressure/potential drop measured between the actual sand-column piezometers. If unavailable, explicitly normalize the fixed-head calculation to its observed clean-bed U0 and disclose that the absolute head/permeability baseline is not independently established by this normalization. Do not claim an absolute hydraulic validation from a normalized baseline.

Not experimentally established by the inspected article: deposit internal porosity eps; two resolved species' inputs/outputs; independent per-species thresholds or capacities; a unique conversion of its broad size distribution into two effective populations. Constant eps=.5 is a model assumption. Temperature-dependent water viscosity/density and the exact bed-scale baseline head should be checked before absolute pressure predictions.

The article reports calibrated turbidimeter-to-mass concentration measurements. Its broad particle distribution and dominant straining/sedimentation mechanisms mean any one-component q is an effective closure. An aggregate fit does not resolve those mechanisms. High-concentration late internal-cake behavior is outside the clean scope; do not fit that interval and then interpret the parameters as depth-only capture. Selection of an early/low-concentration interval must be declared before examining its held-out performance.

## 5 Minimal fitting strategy and identifiability

Use physical geometry, m0, C0 and rho_p from the experiment, not the illustrative manuscript values. Fit the constant-flow concentration trace using only two effective parameters initially: alpha_s and beta_s. Set alpha_d=beta_d=0 as an explicit restricted baseline; fixed U alone cannot separate alpha_s from alpha_d*U. Keep charging multiplier chi and ordered thresholds B1,B2 fixed under a declared protocol, and test their sensitivity afterward rather than silently optimizing them.

Do not fit three or four kinetic parameters simply because an optimizer permits it. Permit cap as a third calibrated parameter only if profile scans show it is constrained by the calibration trace. If the trajectory never reaches its capacity and is insensitive to it, report an assumed cap or bound, not a measured estimate. Conversely, cap enters the aging multiplier alpha*cap/B before arrest: apparent cap identification can therefore be a correlation with alpha and fixed B2. It must survive alternative reasonable fixed threshold choices before receiving a physical interpretation.

Charging not resolved in the data cannot be validated. Start at B=0 at the experimental injection time; never initialize a fitted nonzero B at the first observed sample merely to force agreement. Preset charging thresholds may put that transition within an unobserved early interval, but this must be labeled an untested assumption and profiled for its effect on predictions. Similarly, if neither aging nor arrest is exercised, successful fitting supports only the visited kinetic regime, not all four branches.

Fit several initial guesses and show parameter profiles or a Jacobian/sensitivity diagnosis. A very broad beta profile means deposition and release are not separately identified; a broad cap profile means no justified cap estimate. If a two-parameter fit fails, record its residual structure before expanding the parameter subset. Adding one original dynamic coefficient may be tested using two independently observed operating conditions, but then the fitted condition is calibration rather than held-out validation.

Quantify concentration error and hydraulic error separately. Report concentration RMSE on C/C0, maximum deviation and systematic temporal residuals; avoid relying solely on R2. Weight by experimental/digitization uncertainty only when that uncertainty has actually been quantified. Repeat digitization and report its envelope. Do not invent Gaussian measurement errors from the line thickness.

## 6 Gitis 2010: what its data can and cannot identify

Official DOI: https://doi.org/10.1016/j.cej.2010.07.044 . Author-deposited accepted primary manuscript inspected at https://www.academia.edu/7237022/Deep_bed_filtration_model_with_multistage_deposition_kinetics ; published primary article text also inspected at https://electronicsandbooks.com/edt/manual/Magazine/C/Chemical%20Engineering%20Journal/2010%20Volume%20163/1-2/78-85.pdf .

Figs6-7 concern two laboratory runs with sampling at0.82,1.02,1.62m; they are not six independent experiments. Their media porosity is0.44 and the figure comparison uses U=10m/h. The nominal packed depth in the prose is1.6m, whereas the plotted largest depth is1.62m; use the declared sampling geometry consistently. The paper says these laboratory runs did not reach terminal breakthrough. Its illustrative Fig1 values, including10mg/L feed and1500kg/m3 particle density, are not verified measurements of those EPA runs and must not be transplanted.

The normalized output curves do not, in the inspected account, establish an absolute experimental C0, material density and deposit porosity for reproducing the saturation-storage effects. If the reported ordinate is a turbidity proxy, equating it to mass C/C0 requires the applicable optical calibration. Consult the underlying EPA report or2008 experimental paper before assigning physical mass concentrations.

For fixed U and beta_d=0, normalized equations can be written using y=c/c0, b=B/c0 and zeta=c0/(1-eps). Their storage depends on zeta and kinetic thresholds B1/c0,B2/c0,cap/c0. Thus c0 and eps cannot be identified separately from normalized concentration histories; low-loading fits may be practically insensitive even to zeta. An arbitrary c0 value can define a conditional shape calculation, but it cannot become an inferred experimental fact or validate the original hydraulic/liquid-storage closure.

Use an explicitly conditional minimal fit (effective capture, charging multiplier and normalized charging threshold) only if those features are resolved. Fix or profile detachment if weakly constrained; do not fit unvisited aging/capacity branches. Calibrating on shallow ports and predicting deeper ports tests spatial generalization within a run, not independent experimental replication. Predicting a second run is held out only if its different feed/coagulation/chemistry conditions are known and legitimately accounted for; the source itself used different fitted coefficients between Figs6 and7.

Dimensional rate comparison to Gitis must include normalization. If its bulk retained mass sigma obeys sigma_t=Ka*U*C-Kd*sigma, then with B=sigma/(m0*rho_p), c=C/rho_p the corresponding present coefficients are alpha=Ka*U and beta=m0*Kd, not beta=Kd. Its fitted coefficients are model estimates, not measurements to import as ground truth. Its advection-dispersion model also differs from the authorized dispersion-free equations.

## 7 Held-out tests and declared failure criteria

1. Hold out the entire low-concentration constant-head condition and hydraulic records, rather than random points from one trace. Carry the fitted physical coefficients unchanged; change only the experimentally imposed boundary data. Additional concentration runs may serve as further predictions after defining their admissibility.
2. Use only actual experimental traces/symbols. Published simulation lines, calibrated reference-model parameters, and BTC-derived deposit integrals are not additional independent observations.
3. No dispersion is permitted in the authorized framework. If early arrival or front width is systematically incompatible, report it. A restricted post-front loading fit is permissible only when labeled as such, not as validation of the whole breakthrough curve. Do not hide diffusion by choosing an artificially coarse transport grid.
4. Match final deposit profiles to measured section averages using each section's physical width. Compare corresponding local permeability sections separately. Do not compare local measured permeability with a global column harmonic mean.
5. A direct closure consistency test uses co-located retained fraction and permeability. For local sigma_v per bulk bed, B=sigma_v/m0 and (approximately constant viscosity)

       eps_inferred=1-sigma_v/[m0*(1-(k/k0)^(1/3))^2].

   Correct hydraulic conductivity for the small viscosity ratio if needed before using k/k0. If no single admissible constant0<=eps<1 describes the measurements within their uncertainties, the fixed deposit-volume/permeability closure is incompatible with those measurements independent of q. Since eps=.5 is to remain fixed, a mismatch at that value must be reported rather than repaired on held-out data.
6. At whole-column capacity under continued feed, q=0 everywhere predicts c_out/c_in approaching1. A persistent subunit plateau over a demonstrably capacity-filled bed would contradict the terminal branch. A merely finite-duration experiment below unity does not prove that contradiction: it may not have filled the bed.
7. The absorbing capacity branch cannot release a fully saturated deposit. Do not use backwash or saturated flush release to claim validation of it.
8. Preserve all original identities and mass balances, positivity, and refinement checks at experimental dilute concentrations and long run times. Fitting success does not excuse carrier-liquid negativity, threshold-event artifacts or numerical broadening.

## 8 Reporting language

Use: 'Calibration to digitized aggregate constant-flow observations, with held-out constant-head and hydraulic comparisons in the one-component limit of the unchanged binary model.' If the predictions are inaccurate, say which observable and which regime fail, with errors. For Gitis use 'conditional normalized-concentration comparison' until the missing scale/metadata are established. Never call either exercise validation of all seven experimental saturations, independently resolved binary kinetics, or all four deposition stages unless corresponding observations exist.

## 9 Kinetics-independent mass–hydraulic consistency bound

The unchanged permeability closure implies a useful global necessary condition without assuming a uniform deposit. Put D=B/(1-eps) in the one-component limit and r(D)=(1-sqrt(D))^-3. Define

    h_min=min_{0<D<1} [r(D)-1]/D.

With s=sqrt(D), stationarity gives 2s^3-8s^2+12s-3=0. Its derivative is strictly positive on0<s<1, so the unique root gives the global minimum (the objective diverges at both endpoints):

    s_min=0.3085860202182409,
    D_min=0.09522533187413257,
    h_min=21.269750043077014.

Consequently r(D)>=1+h_min*D pointwise, irrespective of spatial deposit shape. Since mu/mu0>=1, the series resistance yields

    K_app/K0 <= 1/[1+h_min*mean(D)].

Here K_app denotes the hydraulic-conductivity/conductance-derived global ratio measured over the same bed interval. For uniform viscosity it is the effective intrinsic-permeability ratio. The inequality is a lower supporting-line argument; it does not assume global convexity of r, which is in fact not convex near D=0.

For initial-pore-volume count NPV and independently measured retained fraction f_ret=M_retained/M_in,

    mean(B)=c0*NPV*f_ret,
    mean(D)=c0*NPV*f_ret/(1-eps).

The dry-mass retained fraction in Alem's0.25g/L fixed-flow table is f_ret=.6795. If83PV uses initial pore volume and all intervals coincide, c0=.25/2600 gives mean(B)=.0054229326923. The necessary upper ratios are

    .89658403225 for eps=0 (most permissive admissible porosity),
    .81255307015 for fixed eps=.5.

A provisional digitized Figure6 ratio near.936 exceeds both. This would be a kinetics-independent mismatch of the original hydraulic framework, conditional on the measurement definitions and uncertainties. At ratio.936, even eps=0 allows f_ret at most.40280698 with the other inputs fixed; the table gives.6795. No retuning of q can reconcile both exact targets because changing q cannot alter the measured final inventory used in the bound.

Before claiming exclusion, verify that the permeability called 'overall' was measured across exactly the same sand interval used for the nine-section retained mass; entrance-excluding piezometer locations could matter substantially. Also verify the initialPV definition, normalization of initial hydraulic conductivity, and the actual ratio uncertainty. The reported K0 uncertainty alone is not automatically an error bar for K/K0 and does not establish its covariance or repeatability. The unknown-porosity bound differs from the provisional observed ratio by only about.0394, so do not label it statistically robust without a justified measurement interval.

An honest uncertainty test is: compute the smallest physically plausible retained mean B_min from documented input/density/NPV/mass-recovery uncertainties; form the most permissive model upper bound 1/[1+h_min*B_min/(1-eps)] (or eps=0 if seeking exclusion of every nonnegative deposit porosity); compare it to the experimental ratio's documented LOWER bound. Only if the experimental lower bound exceeds this model upper bound is incompatibility established over those uncertainty ranges. If uncertainty information is unavailable, report the nominal discrepancy and the tolerance needed to reconcile it, not a confidence level.
