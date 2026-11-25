# Operating mode and capacity matched controls

These calculations retain the seven saturation balances, fixed deposit
porosities, liquid exchange, Darcy equation, permeability and viscosity of the
2017 model. Only the net deposition closure is varied between the existing
linear/multistage laws and declared ablation controls. No experimental
parameter is inferred here; all numbers are illustrative.

## Reproduction

Python with NumPy and Matplotlib is sufficient. From the parent directory run:

```bash
python run_mode_controls.py
python run_control_refinement.py
```

The first script executes 20 cases and reproduces the two primary figures.
The second executes 12 additional control-refinement cases. Full fields,
transport histories and summary results are stored here. To redraw both figures
without repeating any calculation, from the parent directory run:

```bash
python -c "import run_mode_controls as s; s.redraw_saved_figures()"
```

The included `numerics/reference_solver.py` is an unchanged copy of the preceding
two-component solver. Four static-rate benchmark comparisons give identical
final states before the new cases are run. The original seven-case study and
its previously recorded results remain separate and unchanged.

## Numerical method and API

`numerics/solver.py` advances the original dimensionless time t*, not a reduced
equation written in pore volumes. It integrates injected PV from the common
transport flux. The source stage recomputes hydraulic coefficients after
transport, as in the previous Lie splitting. Consequently, the continuous
operating-mode equivalence is not imposed on the numerical trajectory.

The API `simulate(p,n,tend,inlet,record_dt,...)` retains the earlier arguments.
Additional keywords are:

- `initial`: a four-field array with shape `(2,2,n)`, or a preceding result dict.
- `pv_end`, `record_pv`: stop/sample at cumulative injected pore volumes.
- `inlet_axis='pv'`: interpret the inlet callback against PV instead of t*.
- `store_profiles=True`: return the state at every recorded observation.
- `quality_threshold`: one threshold or two component-specific thresholds.
- `quality_levels`: a list of common thresholds for an admissible-volume curve.

The source donor restriction uses the globally valid attachment bound
alpha_bound = alpha_s + alpha_d*forcing. Under positive fixed pressure,
`k/k0<=1` and `mu/mu0>=1` imply `u*<=Delta p*=forcing`; under imposed flow,
`u*=forcing` exactly. This bound applies to post-transport source coefficients.
The actual source donor ratio is checked at every step. Physical saturation and
mobile-volume bounds are also checked. The latest source includes all-step
derived-z and mobile-fraction diagnostics. For the original 20 archived cases,
only sampled-z extrema are reported because those runs preceded that added
diagnostic. This distinction affects diagnostics, not the state calculation.

## Protocol and observations

All cases start from clean-water-filled pores. Both constant-pressure and
constant-flow experiments start at `Delta p*=10`, `u*=10`, with `mu*=1` initially.
The same constant inlet `(0.05,0.01)`, media porosity `m0=0.35`, fixed deposit
porosities `(0.5,0.3)`, and kinetic thresholds/capacities are used. In particular,
the initial concentration is zero inside the bed; the inlet step occurs at t*=0.

The operating-mode test uses alpha_s=beta_s=0, alpha_d=(0.1,0.1),
beta_d=(0.05,0.05), at a common endpoint of six injected pore volumes. Independent
time integrations use 80, 160 and 320 cells. Maximum state differences decrease
from 9.73e-5 to 3.99e-5; mean trajectory differences decrease from 1.34e-6 to 3.42e-7.
The finite-grid discrepancy is consistent with the independently split time
integrations, rather than being forced to machine zero. The clock-energy
identity defect decreases from 1.32e-4 to 3.25e-5.

The static perturbations add alpha_s=beta_s=eta for eta=0, 0.1, 0.3, 1 at 160 cells,
leaving the same dynamic coefficients. Their final states all coincide after
capacity filling, so the recorded comparison uses the maximum/mean field
differences over the entire PV trajectory. The maximum normalized outlet
difference grows from about 0.000715 at eta=0 to 0.335 at eta=1. This demonstrates
loss of dose collapse for this chosen parameter path, not a universal monotonic
law over all possible static coefficients.

## Kinetic controls

The control runs set alpha_s=beta_s=(1,1), alpha_d=beta_d=0, and stop at 14 PV.
Each operating mode compares:

1. Original linear: q=alpha*c-beta*B, with no terminal capacity.
2. Capacity-matched linear: the same rate for B<C, then q=0 at C; the capacity
   boundary is absorbing, as in the multistage model.
3. Charging without aging: irreversible chi*alpha*c below B1, then the same
   linear reversible rate until the same absorbing capacity. B2 is inactive.
4. Multistage: the complete previously specified charging, transition, aging
   and absorbing-capacity law.

All three capped controls fill both capacities throughout the column by 14 PV.
They consequently have identical terminal resistance 2.3665387823 and deposits
(0.02,0.004). The original uncapped law does not share that terminal inventory.
This removes terminal-capacity confounding when comparing the three capped
models, and isolates transient consequences of charging and aging.

## Quality and energy conventions

The quality criterion requires BOTH component outlet/inlet ratios to be below
the same illustrative threshold theta. Thresholds from 0.5 to1 are supplied;
none is a drinking-water or other regulatory standard. The plotted volume is
the integral of that indicator against injected PV after the first initial PV.
This omits the nominal initial clean-water displacement, but small remaining
windows near its front can still be dominated by numerical spreading.
`control_refinement.json` documents this resolution dependence.

At theta=0.9, the 320-cell multistage/capacity-matched-linear volumes are
1.83673/2.99193 PV under fixed pressure and 0.88450/0.05730 PV under fixed flow.
The small fixed-flow linear volume decreases substantially with refinement;
it should not be described as a resolved physical service window. At theta=0.95,
multistage is worse in BOTH modes: 2.180/3.01670 PV under fixed pressure and
3.17754/5.08477 PV under fixed flow. Thus the result is a threshold- and
operating-mode-dependent tradeoff, not a universal multistage advantage.

`energy_star` is the quadrature integral of u* Delta p* dt*. Dividing it by m0
gives the normalization integral Delta p* dPV. The tabulated
`energy_star_per_good_effluent_PV` uses the former numerator; divide by m0 to
obtain pump work per physical good-volume scale. If admissible volume is zero,
the quotient is undefined and serialized as null/NA. Energy covers the same 14-PV
horizon for all controls, including any off-specification effluent; it is not an
optimized stopping-policy result. Under fixed pressure every case has E*=49.
At 320 cells the flow-controlled multistage and matched-linear values are
107.06961 and 104.80316, respectively. Greater removal at a selected threshold
therefore need not reduce energy.

`figure_operating_mode_invariance` uses the 320-cell pure dynamic trajectories,
the three-grid comparison and the 160-cell static perturbations.
`figure_capacity_matched_controls` uses 160 cells consistently for all four
closures. Its entire threshold curve is supplied rather than choosing a
threshold that favors one closure. The refinement cases are supplementary.

## Verification limits

Maximum core-run species-balance residual is 2.34e-15. All reported physical
saturations remain nonnegative. Across sampled core states, the minimum derived
z is 0.78610 and minimum mobile fraction is 0.88753; the capped-refinement cases
also record nonnegative derived z at every step. Capacity crossing counters are
implementation diagnostics, not a physical observable. A later correction to
the charging-linear counter affects no archived state, flux or metric here;
the control summaries use directly measured final capacity fractions.
