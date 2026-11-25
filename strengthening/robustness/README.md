# Robustness of the multistage kinetic extension

This portable supplement tests whether the four loading-history outcomes survive a continuous transition at the second deposition threshold. The primary sharp multistage law remains the baseline. Only a kinetic sensitivity variant is introduced; every retained 2017 balance, deposit porosity, liquid immobilization relation, Darcy law, permeability law, viscosity law and rate-coefficient relation remains unchanged.

## Reproduction

Requirements: Python 3, NumPy and Matplotlib. From this directory:

```bash
python source_checks.py
python run_smoothing.py
python quality/analyze_quality.py
```

The smoothing command executes 56 history cases (seven widths, two meshes, four histories) and two additional time/source-accuracy refinements. It takes several minutes on the study machine. To reuse complete individual results and redraw/reaggregate without repeating them:

```bash
python run_smoothing.py --reuse
```

For a smaller standalone reproduction:

```bash
python run_smoothing.py --widths 0 0.4 0.6 --grids 80 160
```

The scripts use relative paths and bundled inputs. The optional continuous-quality analysis uses saved control trajectories, with no further PDE simulations; see [quality/README.md](quality/README.md).

## Kinetic sensitivity and numerical method

For $T_i=B_{2i}$ and $delta_i=fT_i$, let $L_i=T_i-\delta_i$ and $H_i=T_i+\delta_i$. The band lies strictly between the first threshold and the capacity. Within and above the transition stage, the alternative kinetic expression is

$$
q_i^\delta=\alpha_i c_i\left[1+w_i^\delta(B_i)\left(\frac{C_i}{B_i}-1\right)\right]-\beta_i B_i,
\qquad B_{1i}<B_i<C_i,
$$

where $s_i=(B_i-L_i)/(2\delta_i)$ and

$$
w_i^\delta=
\begin{cases}
0,&s_i\leq0,\\
3s_i^2-2s_i^3,&0<s_i<1,\\
1,&s_i\geq1.
\end{cases}
$$

The original charging branch, generalized holding at $B_{1i}$ and absorbing capacity are preserved. This smooths only the second-threshold jump, not every discontinuity of the closure. At $f=0$, the solver calls the archived sharp source map directly. No near-threshold absolute tolerance selects a branch in the smoothed expression; its integration tolerances control ODE accuracy only.

The conservative finite-volume transport and the posttransport coefficient/concentration freezing are the same as in the sharp solver. Exterior charging, transition and aging paths are advanced analytically until a step encounters the blend. Adaptive SSP-RK3 with step doubling then advances the smooth scalar ODE for the remainder. The default local accuracy settings are relative $10^{-9}$ and absolute $10^{-14}$. A bound $\mathcal L$ on the smooth source derivative enforces $h\mathcal L\leq1/2$ for each substep. The attachment-factor bound is enlarged to

$$
M_i^\delta=\max\left(1,\chi_i,\frac{C_i}{T_i-\delta_i}\right).
$$

The original outer donor bound uses this factor and a global attachment-coefficient bound. Nonnegative SSP stages and this donor bound control deposit gain. The source exchanges precisely $\Delta B_i$ between mobile and deposited inventories; it never clips a mobile inventory independently. Reaching the first threshold or capacity completes the same physical boundary event as the sharp baseline. The original admissibility diagnostics are retained.

`source_checks.py` compares the source map with an independent fixed-step classical RK4 reference using 4096 and 8192 substeps. It also checks the exact sharp dispatch, convergence of a threshold-crossing map as the smoothing width decreases, nonnegative source donors, and paired conservation. These checks isolate the frozen source integration. Full-PDE CFL refinement separately measures effects of transport, coefficient freezing and Lie splitting.

## Protocol

All histories start with clean water and zero mobile/deposited particles. The pressure drop is 10; $m_0=0.35$, deposit porosities are $(0.5,0.3)$, static attachment/detachment coefficients are one, and all dynamic rate coefficients are zero. Thresholds are $B_1=(0.001,0.0002)$, $B_2=(0.004,0.001)$, capacities $C=(0.02,0.004)$ and charging factors $(0.5,0.5)$.

Until $t^*=2$, each species independently receives either its final feed $(0.002,0.0005)$ or high preconditioning feed $(0.006,0.0015)$. The four LL, HL, LH and HH histories then share the same final feed until $t^*=10$. The relative halfwidths are $f=0,0.01,0.05,0.1,0.2,0.4,0.6$ on 80 and 160 cells. All these bands avoid both the first threshold and capacity. For HH at $f=0.4$ and $0.6$, the 160-cell study is repeated with CFL 0.4 rather than 0.8 and source tolerances tightened to $10^{-11}$ and $10^{-16}$.

The wider $f=0.6$ case is included because the local root calculation predicts loss of the low branch; it is a sensitivity scenario, not a fitted transition width.

## Results

The sharp-limit runs reproduce all eight saved baseline final states and recorded flux/deposit trajectories exactly. All 58 calculations conserve each particle inventory to a maximum residual of $2.51\times10^{-13}$ and retain nonnegative tracked saturations; the minimum mobile fraction is 0.97459.

For the tested $f\leq0.4$, all four histories still select the original four interior combinations. Stable deposits are unchanged because they remain outside the smoothing band. Transients and unstable basin separators do change. At $f=0.4$, the local separators are approximately $(0.002974535,0.000790492)$.

The local low and unstable roots collide at $f_{c1}=0.5218587662$ and $f_{c2}=0.5299812290$. At $f=0.6$, all four full-column histories converge to approximately $(0.006322899161,0.001368017062)$ with flux $6.7532268905$. The 160-cell maximum final source residual is $2.90\times10^{-13}$ and remaining flux spread is $1.14\times10^{-11}$. Thus the memory claim has a measured sensitivity limit: it persists for the tested moderate smoothings and disappears for the broader variant that removes the low roots.

Grid refinement and the representative time/source-accuracy refinements preserve both outcomes. Final-state differences are tiny because the trajectories approach uniform roots; the aggregate file also reports transient differences, which are not zero. The largest representative 80-to-160-cell transient mean-deposit difference is $5.78\times10^{-6}$, and the CFL/source-tolerance comparison gives $6.56\times10^{-6}$.

Independent frozen-source tests give a maximum default-map error $2.26\times10^{-10}$ against tiny-step RK4 for a deliberately long source duration 0.01, decreasing to $7.21\times10^{-12}$ under tighter tolerances. The source tests and full-PDE refinement address different errors and are reported separately.

The saved-control continuous exceedance analysis confirms greater multistage exceedance than the capacity-matched linear control at illustrative thresholds 0.8, 0.9 and 0.95, in both operating modes and on all three grids. No universal multistage quality benefit is claimed.

## Files and provenance

- `numerics/reference_solver.py`: immutable copy of the sharp solver used for the preceding history study.
- `numerics/smooth_solver.py`: numerical fork with smoothing parameters, source dispatcher and the required enlarged donor bound.
- `numerics/smooth_source.py`: the smooth kinetic expression and conservative adaptive source map.
- `inputs/sharp_baseline/`: eight saved sharp-history CSVs and final-field files used for a direct reproduction comparison.
- `source_checks_results.json`: independent frozen-source checks.
- `results/smoothing_results.json`, `smoothing_summary.csv`, `local_root_map.csv`: aggregate metrics, individual summaries and local roots.
- `results/smooth_*.csv`, `*_fields.npz`, `*.json`: every simulated time series, final profile and diagnostic record.
- `results/figure_smoothing_robustness.png` and `.pdf`: local roots and full-column loading-history comparison.
- `quality/`: independent continuous exceedance analysis of previously saved capacity-matched controls.
- `manuscript.md`: ready-to-insert methods, results and interpretation; no assigned equation numbering.

The local source classification uses prescribed final concentration and static rates. It is not a theorem about stability or uniqueness of the coupled PDE. The absorbing capacity remains an additional boundary state. Different preconditioning feeds supply different particle masses, so no matched-inventory hydraulic advantage is inferred.
