# Accelerated scalar experimental-comparison solver

This is the one-component specialization of the current, unchanged 2017
saturation framework, with either its original linear net deposition law or
the existing multistage net law. It introduces no new dispersion, structure,
or deposition state. Experimental fitting is intentionally separate.

Files:

- `validation_solver.py`: SciPy-compatible Python API.
- `validation_kernel.cpp`: complete accelerated integration source.
- `benchmark_validation_solver.py`: comparison against the current two-species
  Python solver with the second species inert.
- `validation_solver_benchmark.json`: observed benchmark results.

The wrapper builds `validation_kernel.so` with `g++ -O3 -std=c++17 -shared -fPIC`
on first use, or when the C++ source changes. Only NumPy and a C++17 compiler
are required; no Numba or extra physical model is used.

## Calling the API

```python
import numpy as np
from validation_solver import Params, simulate

p = Params(m0=0.37, eps=0.5, mode="rate", forcing=1.0,
           alpha_s=1.0, beta_s=1.0, charging=0.5,
           B1=0.001, B2=0.004, cap=0.02)
r = simulate(p, sample=np.array([0., 1., 10., 40., 83.]),
             sample_axis="pv", cin=0.25/2600,
             depths=(0.25, 0.5, 0.75, 1.0), n=80)
predicted_normalized_concentration = r["c_normalized"]
```

Here 0.25 g/L = 0.25 kg/m3, so particle volume concentration is
`cin=0.25/2600=9.615384615e-5` for particle density 2600 kg/m3. There is **no**
additional factor of 1000. Porosity 0.37 is separate from fixed deposit porosity
`eps`. Setting `eps=0.5` is a structural assumption unless independently measured;
it must not be reported as an experimentally measured deposit property.

The time variable remains `t*=U0*t/L`. The porosity factor is explicit in the
state equations: `m0*a_t*+(u*c)_x*=-q`, `m0*B_t*=q`. Cumulative injected pore
volumes are integrated as `PV=integral(u*/m0 dt*)`. In fixed unit flow,
`t*=m0*PV`; under fixed pressure this relation is not valid and PV is integrated.

For pressure control choose `mode="pressure"`. The scalar `forcing` is the
dimensionless pressure drop. To match the measured initial clear-fluid velocity,
choose the velocity reference as that initial velocity and use forcing 1.0.
This is an explicitly stated normalization, not inference that all separately
published initial head, conductivity and discharge data are mutually exact.

The returned dictionary contains:

- `c_depth` and `c_normalized`: concentrations at requested depth fractions;
  shape is number of observation times by number of depths.
- `time`, `pv`, `flow`, `pressure`: original dimensionless time, cumulative PV,
  Darcy rate and pressure drop.
- `Bmean`: spatial average deposited-particle saturation.
- `B_profile`, `a_profile`, `c_profile`, `phi_profile`, `S1_profile`,
  `S3_profile`, and `x`: full fields at requested observations.
- `k_profile`, `mu_profile`, `pressure_profile`: relative local permeability,
  relative viscosity, and pressure at cell centers with the outlet pressure
  datum set to zero.
- `keff`: intrinsic effective permeability relative to clean-bed permeability,
  computed as `1/mean(1/krel)`.
- `hydraulic_conductance`: `1/mean(mu_rel/krel)`. This includes the suspension
  viscosity and is distinct from intrinsic relative permeability.
- `input`, `output`, `inventory`, `max_balance_error`: particle ledger in
  original-pore-volume units. Influx and efflux are integrals of `u*c/m0` over
  dimensionless time. Multiplication by m0 gives bulk-volume storage units.
- `minimum_saturation`, `minimum_carrier`, `maximum_concentration`, `steps`,
  `hold_events`, `cap_events`: admissibility and event diagnostics.

Full retained profiles can be converted to particle mass per bulk volume as
`m0 * particle_density * B_profile`. Conversion to a mass per unit bed length
also requires the known column cross-sectional area. Never compare these
different experimental units directly without conversion.

## Numerical identity and speed

The integration algorithm is the current first-order conservative upwind
transport plus exact frozen-concentration source branches, with explicit
threshold events, holding completion at the attractive B1 discontinuity and
terminal arrest at capacity. The second species is simply absent. Source
increments cancel exactly between mobile and retained inventories. A failed
donor/admissibility check halves and retries the step; it does not clip states.
Constant deposit porosity and the original permeability/viscosity relations are
unchanged. Original velocity/gradient-dependent kinetic coefficients remain
available through alpha_d and beta_d; default calibration should keep them zero
unless the experiment can separately identify them.

Wrapper tests at 40 cells show agreement with the existing Python solver at
roundoff for linear and multistage kinetics under both fixed-flow and fixed-head
operation, including nonzero original dynamic coefficients. See the JSON
benchmark for exact errors. One illustrative 80-cell run through 83 PV requires
approximately 0.016 seconds after compilation in this environment. This is a
speed measurement, not an experimental fit or a general hardware guarantee.

## Calibration discipline

The existing kinetic parameters are alpha_s, beta_s, charging, B1, B2 and cap,
plus optional original dynamic coefficients alpha_d and beta_d. Fix measured
geometry, m0, feed concentration and operating mode. Fix or sensitivity-test eps
as a constant; fitting eps and capacity simultaneously without independent data
can hide non-identifiability. Begin with only 2–4 independently supported kinetic
parameters. Do not fit every threshold, detachment coefficient and hydraulic
quantity to the same sparse breakthrough curve.

Parameters constrained by one-component data validate only the scalar limit.
Independent-depth profiles, operating-mode changes, hydraulic data or separate
runs should be withheld from fitting when used as validation. Such validation
does not establish two-component interaction accuracy.
