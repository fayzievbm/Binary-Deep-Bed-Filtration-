# Kinetics-only extension of the 2017 model

Run `python solver.py` using NumPy and Matplotlib. This package does not use the
earlier structural-maturation model. All figures, metrics and results in this
directory refer ONLY to a piecewise modification of net deposition kinetics.
Parameters are illustrative, dimensionless and uncalibrated.

## Unchanged equations

For species i=1,2 the evolved state is `a_i=S_2i` and `B_i=S_4i`. The original
time scale `t*=U0 t/L` is retained, including the factor `m0=0.35`:

`m0 a_i,t* + (u* c_i),x* = -q_i`,

`m0 B_i,t* = q_i`.

Deposit porosities remain constant: `eps=(0.5,0.3)`. Thus
`S_3i=eps_i B_i/(1-eps_i)`, `D=sum[B_i/(1-eps_i)]`, `phi=1-D`,
`c_i=a_i/phi`, and `S1=phi-a1-a2`. These are algebraic consequences of the
original seven saturation balances, not new closures. The original permeability,
viscosity, Darcy equation and concentration definitions are unchanged:

`k/k0=(1-sqrt(D))^3`, `mu/mu0=1+2.5(c1+c2)`,

`u*=Delta p*/mean(mu_rel/k_rel)` for fixed pressure, or
`Delta p*=u* mean(mu_rel/k_rel)` for prescribed flow. The computational column
is 0<=x*<=1, with inlet pressure 10 and outlet pressure 0 for pressure-controlled
examples. Rate-controlled examples prescribe u*=10. Shared scalar hydraulic
velocity is used at every transport face. Local |p*_x*|=u* mu_rel/k_rel enters
the original dynamic coefficient relations if their coefficients are nonzero:

`alpha_i*=alpha_si*+alpha_di* u*`,

`beta_i*=beta_si*+beta_di* |p*_x*|`.

Reported coefficients are alpha_s=(1,1), beta_s=(1,1), alpha_d=beta_d=(0,0).
Original dimensional reference scales may be assigned as L=1 m, U0=10^-5 m/s,
k0=10^-12 m^2 and mu0=10^-3 Pa s, giving a time scale of 10^5 s and a pressure
scale of 10^4 Pa. The introduced thresholds and charging factors are illustrative.

## The only model change

The original linear net rate `q_i=alpha_i c_i-beta_i B_i` is replaced by the
IJAM 2020 Eq.8 form with separate species thresholds:

1. Charging, 0<=B_i<=B1_i: `q_i=chi_i alpha_i c_i`.
2. Transition, B1_i<B_i<=B2_i: `q_i=alpha_i c_i-beta_i B_i`.
3. Aging, B2_i<B_i<Bcap_i: `q_i=alpha_i(Bcap_i/B_i)c_i-beta_i B_i`.
4. Terminal capacity, B_i=Bcap_i: `q_i=0`.

Terminal arrest applies to the full net rate, including detachment, as in the
stated piecewise law. It is not an attachment-only capacity limiter. There are
no new retained states, maturation-transfer rates, varying deposit porosities,
shared-capacity functions or flush-memory experiments.

Defaults are B1=(0.001,0.0002), B2=(0.004,0.001), Bcap=(0.02,0.004), chi=(0.5,0.5).
The independent caps obey `sum[Bcap_i/(1-eps_i)]=0.0457142857<1`, yielding
`phi>=0.9542857143` in the multistage model.

## Discontinuous closure completion

At B1, the charging-side rate can be positive while the transition-side rate
is negative. A literal evaluation of alternating branches can chatter. The
event solver uses the attractive-interface completion: hold B=B1 while
`q_left>=0>=q_right`, then leave the interface when the rates permit crossing.
This explicitly specified interpretation completes the discontinuous kinetic
law; no general classical uniqueness is asserted. The equality value of a
branch at an instant of transversal upward crossing does not consume a finite
time. At an exactly reached B2, a nonpositive lower-branch rate selects the
transition branch (including downward crossing and the repelling interface);
a positive lower-branch rate continues into aging. Terminal capacity is
absorbing, exactly as prescribed. The main constant
feed runs produced zero attractive-interface hold events; a separate local
benchmark verifies holding and terminal arrest under zero suspended concentration.

## Conservative numerical algorithm

The solver uses first-order upwind finite-volume transport followed by a local
reaction update (Lie splitting). The transported mobile concentration and
coefficients are frozen during each source substep. Each kinetic branch is then
integrated analytically, with the remaining source-substep time split exactly
at B1, B2 and Bcap events. In particular, the charging branch is linear in time,
the transition branch is exponential, and the aging branch has

`B(t+dt)^2 = alpha Bcap c/beta + [B(t)^2-alpha Bcap c/beta] exp(-2 beta dt/m0)`.

The deposit increment is applied with the opposite sign to mobile inventory:
`a_new=a_old-(B_new-B_old)`. Thus each species' reaction inventory is conserved.
Transport inlet/outlet ledgers use the very same fluxes as the state update.
The overall method is first order in space and time, despite exact frozen-source
branch integration. Freezing concentration is a numerical approximation, not a
new constitutive law or smoothing of the discontinuous rates.

The time-step safety factor is 0.8. Transport requires
`dt <= m0*min(phi)*dx/u`; a source donor bound is
`dt <= min[m0*phi/(alpha_i*max(1,chi_i,Bcap_i/B2_i))]`.
No independent clipping of mobile or deposit state is used. Threshold and cap
values are set by calculated event crossings. The square-root guard in the
analytical B^2 formula only guards floating-point arithmetic; the evolved
states are checked and any materially negative state triggers an error.

## Experiments and verification

Initial pores contain clean carrier fluid. Constant inlet c=(0.05,0.01) is
supplied through t*=0.8. Seven main 320-cell runs compare original linear and
multistage kinetics under fixed pressure and flow, vary chi1 to 0.25 and 1,
and vary the second-species alpha_s2 to 3. In the parameter comparison, all
other model equations and parameters are fixed.

For frozen inlet concentration, the exact stage-event times are
species 1: 0.014, 0.03611261557, 0.12268445789;
species 2: 0.014, 0.04380523292, 0.12876910068.
Source integration is checked against these analytical branch solutions, and
the original linear source against B_i=alpha_i c_i/beta_i*(1-exp(-beta_i t*/m0)).
These local checks verify source integration, not the complete PDE solution.
An additional downward B2 crossing benchmark checks that an aging trajectory
switches back to the transition branch at the computed event time. A coupled
steady benchmark compares both hydraulic modes against the exact spatially
uniform capped state: Dcap=0.045714285714, c=(0.05,0.01), relative viscosity1.15,
u*=10*kcap/1.15 and Delta p*=10*1.15/kcap.

Space refinement uses 80,160,320 cells against projected cell averages of a
640-cell solution at t*=0.2. The inlet is smoothly ramped over 0<=t*<=0.1.
Time refinement uses 80 cells and step caps 0.0002,0.0001,0.00005 against a
fixed-grid reference step cap 0.0000125. Errors are means of absolute errors
across the four independent state fields. Both studies show decreasing errors
consistent with the first-order method and a finite numerical reference.

`summary.json` contains complete parameters, admissibility extrema, conservation
residuals, fixed-porosity identity residuals, event counts and refinement errors.
`case_metrics.csv` summarizes cases. Per-case CSVs contain breakthrough curves,
hydraulic observables and inventories. NPZ files contain final fields and
snapshots. Four figures are available in PNG and vector PDF.

Every inventory in the data files is normalized by original pore volume. The
balance residual is `integral(a+B) - integral_0^t(u c_in/m0)dt +
integral_0^t(u c_out/m0)dt`. Multiplying by m0 gives bulk-volume storage units.
This is numerical verification, not experimental validation.

For the alpha_s2=3 contrast case, a separate 80/160/320/640-cell refinement checks
the outlet species-2 peak, sampled every 0.0005 through t*=0.25. Individual outlet
concentrations need not obey a bound by their own feed concentration: when that
species is already at terminal capacity, capture of the other species can
immobilize carrier liquid and increase its mobile concentration. Total mobile
concentration and all seven original saturations are monitored independently.
