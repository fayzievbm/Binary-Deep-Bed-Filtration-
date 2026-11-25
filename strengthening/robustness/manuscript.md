## Robustness to smoothing the second kinetic transition

The upward jump at the second deposition threshold creates a repelling interface in the primary multistage closure. To determine whether the observed loading-history effect requires an exactly discontinuous switch, we considered a kinetic sensitivity variant with a continuous transition between the transition and aging expressions. Every balance and nonkinetic constitutive relation inherited from the 2017 filtration framework was retained. The sharp closure remained the primary model; smoothing was not calibrated to experimental data.

Let $T_i=B_{2i}$, $\delta_i=fT_i$, and $s_i=(B_i-T_i+\delta_i)/(2\delta_i)$. For $B_{1i}<B_i<C_i$, the sensitivity variant is

$$
q_i^\delta=\alpha_i c_i\left[1+w_i^\delta(B_i)\left(\frac{C_i}{B_i}-1\right)\right]-\beta_iB_i,
$$

where

$$
w_i^\delta(B_i)=
\begin{cases}
0,&s_i\leq0,\\
3s_i^2-2s_i^3,&0<s_i<1,\\
1,&s_i\geq1.
\end{cases}
$$

Thus the original transition and aging laws apply outside $[T_i-\delta_i,T_i+\delta_i]$, and their blend is continuously differentiable at its edges. The charging law, generalized holding at $B_{1i}$ and absorbing capacity remain unchanged. Every tested band lies strictly between the first threshold and capacity. At $f=0$, the original sharp source map is used directly.

The four clean-start LL, HL, LH and HH inlet histories were repeated at $f=0,0.01,0.05,0.1,0.2,0.4,0.6$ on 80 and 160 cells. Their pressure drop, preconditioning duration, common final feed and final time were the same as in the primary history study: $\Delta p^*=10$, preconditioning until $t^*=2$, final feed $(0.002,0.0005)$ and termination at $t^*=10$. Static attachment and detachment coefficients equal one; dynamic coefficients equal zero. Source concentrations and rate coefficients were frozen after transport, as in the sharp calculation. Exterior source paths were integrated analytically until encountering the blend; adaptive SSP-RK3 with step doubling then advanced the smooth ODE for the remainder of that source step. The relative and absolute local tolerances were $10^{-9}$ and $10^{-14}$, with an additional source-derivative step bound. The mobile and deposited inventories exchange precisely the same increment. The donor bound uses the attachment factor

$$
M_i^\delta=\max\left(1,\chi_i,\frac{C_i}{T_i-\delta_i}\right),
$$

which accounts for the larger attachment factor at the lower edge of a smoothed band. No independent clipping of mobile inventory or tolerance-based selection of a second-threshold branch was used. Two 160-cell HH calculations at $f=0.4$ and $0.6$ were repeated with the CFL factor reduced from 0.8 to 0.4 and source tolerances tightened to $10^{-11}$ and $10^{-16}$.

The eight sharp-limit repetitions reproduced the previously saved final states, recorded fluxes and mean deposits exactly. An independent fixed-step classical RK4 calculation of the frozen smooth source gave maximum errors of $2.26\times10^{-10}$ at the default tolerances and $7.21\times10^{-12}$ at the tighter tolerances in source stress tests with duration $0.01$. The independent RK4 calculations with 4096 and 8192 substeps differed by at most $2.14\times10^{-16}$. A separate threshold-crossing source test approached the sharp result as the smoothing halfwidth decreased. These source checks are distinct from transport/source-freezing and splitting errors, which were assessed by the full-PDE refinements.

For every tested width through $f=0.4$, the two stable local deposit levels remain outside the blend, so their values are unchanged. The four clean-start histories continued to select all four combinations. At 160 cells, the spreads across histories remained $0.00432455532$ and $0.000914213562$ for the two mean deposits, and $1.28156786704$ for the Darcy flux. This persistence concerns the selected final plateaus: smoothing changes the transient kinetics and the basin separator. For example, at $f=0.4$ the local unstable roots are approximately $0.002974535$ and $0.000790492$, below the corresponding sharp thresholds $0.004$ and $0.001$.

The same calculation also identifies a loss of this multiplicity under sufficiently broad smoothing. With the final concentration prescribed, the low stable root and the intermediate unstable root meet in saddle-node collisions at

$$
f_{c1}=0.5218587662,\qquad f_{c2}=0.5299812290.
$$

These values follow from simultaneously imposing $q_i^\delta=0$ and $\partial q_i^\delta/\partial B_i=0$; they describe the local static-rate source, not a general bifurcation or stability theorem for the coupled PDE. At $f=0.6$, only one stable interior local root remains for each species. All four full-column histories then converged by $t^*=10$ to the same deposits, approximately $(0.006322899161,0.001368017062)$, and flux $U^*=6.7532268905$ on both meshes. On 160 cells, the remaining history-to-history flux spread was $1.14\times10^{-11}$ and the maximum absolute final net source was $2.90\times10^{-13}$.

Across all 58 calculations, the maximum species inventory-balance residual was $2.51\times10^{-13}$, all tracked saturations were nonnegative, and the minimum mobile fraction was 0.97459. The maximum final deposit distance to a stable local root was $2.66\times10^{-13}$. For the representative $f=0.4$ and $0.6$ histories, the maximum transient mean-deposit difference between 80 and 160 cells was $5.78\times10^{-6}$; their final-state differences were below $4.67\times10^{-15}$. Halving the CFL factor and tightening the source tolerance changed the representative transient mean deposits by at most $6.56\times10^{-6}$ and their final states by less than $1.93\times10^{-15}$. Final-state agreement alone does not imply that the complete transients are identical.

The history effect therefore persists for the tested moderate smoothings and is not solely a consequence of evaluating an exactly discontinuous second-threshold switch. Its survival nevertheless depends on the kinetic shape: broader admissible smoothing removes the low roots and the tested histories lose their distinct final plateaus. The defensible interpretation is history selection within a specified multistage closure and its tested smooth variants. These calculations do not establish experimental hysteresis, general coupled-PDE stability, or a hydraulic benefit at matched retained inventory. The absorbing capacity remains an additional boundary state.

**Suggested figure caption.** Kinetic smoothing robustness at unchanged filtration physics. (a,b) Stable and unstable local interior roots at the prescribed final concentration; open symbols show final mean deposits from the four 160-cell clean-start histories. Vertical dotted lines and black points identify the local saddle nodes. (c) The four histories retain distinct flux plateaus at relative smoothing halfwidth 0.4. (d) They converge to a common plateau at halfwidth 0.6. Shading marks the preconditioning interval. The source-law calculation classifies local roots; full-column trajectories provide numerical history-selection evidence under the stated protocols.

The complementary continuous effluent-exceedance analysis is supplied separately in `quality/manuscript.md`; its metric avoids relying solely on the duration of small threshold-crossing windows.
