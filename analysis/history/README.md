# History-selected interior equilibria

The four calculations start from the same clean-water-filled bed and zero mobile or retained particle saturations. All impose the same dimensionless pressure drop 10, use static attachment and detachment coefficients equal to one, and retain the original seven-saturation balances and fixed-porosity hydraulic relations. Only the preconditioning inlet history differs.

For 0 <= t* < 2, a selected species receives the higher inlet volume fraction 0.006 (type 1) or 0.0015 (type 2). An unselected species receives its final value throughout. At t*=2, all four runs use the same inlet (0.002,0.0005). They are evolved to t*=10 on 80 and 160 cells. The time is the original t*=U0*t/L, not injected pore volumes.

| History | Preconditioning inlet | Final retained type 1 | Final retained type 2 | Final Darcy flux |
|---|---|---:|---:|---:|
| LL | (0.0020, 0.0005) | 0.0020000000 | 0.0005000000 | 8.0281958012 |
| HL | (0.0060, 0.0005) | 0.0063245553 | 0.0005000000 | 6.8744835440 |
| LH | (0.0020, 0.0015) | 0.0020000000 | 0.0014142136 | 7.7994725069 |
| HH | (0.0060, 0.0015) | 0.0063245553 | 0.0014142136 | 6.7466279342 |

At the common final inlet, each species has an interior transition equilibrium B_L = c and an interior aging equilibrium B_H = sqrt(Bcap*c). The four computed final states reproduce the four combinations of these two roots. They are not the only equilibria permitted by the complete source law: terminal capacity is absorbing and permits additional boundary states, which the selected protocols do not visit.

At 160 cells, the largest absolute final local deposition source over either species and the bed is 1.005e-14, the largest deviation of a retained profile from its exact uniform interior root is 1.066e-14, and the largest final Darcy-flux error is 1.975e-12. No capacity event or B1 holding event occurs in any run. The largest accumulated species-inventory residual is 2.482e-13. The minimum reconstructed carrier saturation over all runs is 0.967284.

The maximum 80-versus-160 final-state difference is 2.802e-16, while the largest difference during the recorded transients is 0.003623 in Darcy flux and 6.054e-06 in retained saturation. Agreement of the terminal branch selection is therefore much stronger than transient pointwise agreement; it is supported by the independently calculated uniform roots. No 320-cell run was needed to distinguish the branches in this example.

The preconditioning protocols introduce different particle inputs and leave different retained inventories. The four hydraulic plateaus follow the original permeability and viscosity relations and are not an improvement at matched retained mass. The local source phase portraits establish local fixed-concentration signs around the interior roots; these finite-column trajectories are not a general stability proof for the coupled PDE. No experimental validation or first-ever claim for history dependence is made.

Figure history_multistability shows the allowed inlet waveforms, both retained means and the shared-pressure Darcy response. Shading marks preconditioning. LL denotes neither species preconditioned, HL type 1 only, LH type 2 only, and HH both. Coincident retained traces are intentionally overlaid.

Figure history_phase_portraits shows the two local sources at the final prescribed concentration, with the axes scaled as 10^3*B and 10^3*q*. The dotted line marks B2, L and H mark the two interior roots. The final panel compares exact plateau fluxes (bars) and numerical final fluxes (crosses). Capacity boundary states are outside the displayed source intervals.

`run_history.py` regenerates the eight calculations, per-run CSV histories, final-field NPZ files, summary JSON and two PNG/vector-PDF figures.
