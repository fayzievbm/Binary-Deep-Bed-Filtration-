## Coupled uniform equilibria with the original rate coefficients

The local roots in the preceding subsection prescribe the kinetic coefficients. Uniform equilibria can also be characterized while retaining their original Darcy dependence. Write $b_i=B_{i1}$, $T_i=B_{i2}$ and $\mu_\gamma=1+(5/2)\sum_i\gamma_i$ for a fixed positive uniform feed $c_i=\gamma_i$. Retain nonnegative original rate coefficients with $(\beta_{si},\beta_{di})\ne(0,0)$, the capacity condition and $\sum_iw_i\gamma_i\le1$. Let L denote $b_i<B_i<T_i$ and H denote $T_i<B_i<C_i$; all inequalities here are strict.

**Proposition (branchwise uniform-equilibrium uniqueness).** For a prescribed pressure drop $P>0$ or Darcy flux $u_f>0$, each L/H assignment has at most one uniform equilibrium inside its two assigned branches. Define

$$r_i^P(\kappa)=\frac{\gamma_i(\alpha_{si}+\alpha_{di}P\kappa/\mu_\gamma)}{\beta_{si}+\beta_{di}P},\qquad r_i^u(\kappa)=\frac{\gamma_i(\alpha_{si}+\alpha_{di}u_f)\kappa}{\beta_{si}\kappa+\beta_{di}u_f\mu_\gamma}.$$

For the chosen control and assignment $\sigma$, put

$$\mathcal B_i^{\sigma_i}(\kappa)=\begin{cases}r_i(\kappa),&\sigma_i=L,\\\sqrt{C_ir_i(\kappa)},&\sigma_i=H,\end{cases}\qquad G_\sigma(\kappa)=\sum_iw_i\mathcal B_i^{\sigma_i}(\kappa)-(1-\kappa^{1/3})^2.$$

A zero of $G_\sigma$ in $(0,1)$ gives the selected equilibrium if and only if its deposits satisfy both branch inequalities. **Proof.** At a uniform steady state $q_i=0$. Substitution of $u=P\kappa/\mu_\gamma$ or $|p_\xi|=u_f\mu_\gamma/\kappa$ gives the displayed ratios. The L/H zero-source equations give $\mathcal B_i$, and the retained permeability law gives $G_\sigma=0$. Each ratio, and hence each deposited candidate, is continuous and nondecreasing in $\kappa$. The term $-(1-\kappa^{1/3})^2$ is strictly increasing, so $G_\sigma$ has at most one zero. An accepted root reconstructs all seven saturations as in Section 3.5, makes every source vanish, and determines a linear pressure profile. This proves existence and uniqueness for that branch assignment, without a PDE stability or overall uniqueness claim.

For purely dynamic rates with $\beta_{di}>0$ and $\alpha_{di}\gamma_i>0$, both controls reduce to

$$d_i=\frac{\alpha_{di}\gamma_i}{\beta_{di}\mu_\gamma},\qquad \sum_{i:\sigma_i=L}w_id_i\kappa+\sum_{i:\sigma_i=H}w_i\sqrt{C_id_i\kappa}=(1-\kappa^{1/3})^2.$$

Here $G_\sigma(0)=-1<G_\sigma(1)$, so each assignment has exactly one algebraic candidate. Actual branch admissibility remains essential. Define $\ell_\sigma$ as the maximum of $b_i/d_i$ for L components and $T_i^2/(C_id_i)$ for H components, and $h_\sigma$ as the minimum of 1, $T_i/d_i$ for L components and $C_i/d_i$ for H components. Then existence inside the selected open branches is equivalent to

$$\ell_\sigma<h_\sigma,\qquad G_\sigma(\ell_\sigma)<0<G_\sigma(h_\sigma).$$

If any $d_i$ vanishes, the divided test is not used and that component has no positive L/H interior equilibrium. With general static/dynamic rates, endpoint signs must instead be checked for the appropriate ratio; an unconstrained candidate need not always exist. Replacing L by H raises the deposited candidate at an accepted L root and therefore lowers the permeability of a corresponding accepted H assignment. Thus LL has the highest and HH the lowest permeability when all four states coexist; LH and HL have no universal ordering. Threshold, absorbing-capacity and nonuniform states are not counted. The supplementary derivation gives the complete assumptions, boundary parameterizations and numerical root enclosures.

## Suggested results insertion: Darcy feedback shifts coexistence regions

The new equilibrium calculation uses the existing porosities, thresholds and capacities, final feed $(0.002,0.0005)$, and dynamic coefficients $\alpha_{di}=0.1$, $\beta_{di}=0.05$. The strictly interior candidates below show that all four L/H assignments can coexist under the shared permeability feedback. Their fixed-pressure fluxes and fixed-flow pressure drops differ because their retained inventories differ.

| Assignment | Relative permeability $\kappa$ | Darcy flux at $P=10$ | Pressure drop at $u_f=10$ |
|:--|--:|--:|--:|
| LL | 0.766671784 | 7.619098472 | 13.124912398 |
| LH | 0.749184529 | 7.445312086 | 13.431270421 |
| HL | 0.670565723 | 6.664007184 | 15.005986225 |
| HH | 0.659429028 | 6.553331954 | 15.259413181 |

![Suggested figure. Admissible homogeneous interior branch combinations under coupled Darcy feedback: (a) purely dynamic rates, (b) mixed original rates at pressure drop $P=10$, and (c) mixed rates at Darcy flux $u_f=10$. Both map axes vary the original dynamic attachment coefficients at fixed feed, viscosity and $\beta_{di}=0.05$. Mixed panels use $\alpha_{si}=\beta_{si}=0.1$. Solid, dashed, dotted and dash-dotted curves delimit the LL, LH, HL and HH branch families. Stars mark $\alpha_{di}=0.1$. Colors count only strictly interior homogeneous assignments; zero does not mean that the closure has no equilibrium. Threshold, capacity and nonuniform states are excluded.](results/equilibrium_regime_map.png)

The map evaluates each of the four scalar candidates on a $321\times321$ grid over $0\le d_i/T_i\le2$. Shared hydraulic feedback shifts the admissibility boundaries but does not remove the four-state region. Mixed static/dynamic rates also make the branch families depend on operating mode: the two mixed panels have different admissible sets at 11,110 of 103,041 sampled parameter points. These counts characterize this parameter grid, not probabilities. For example, at $(d_1/T_1,d_2/T_2)=(1.4625,0.33125)$ all four interior assignments exist under pressure control, while only HL is accepted under flow control. These are equilibrium-existence differences, not demonstrated transient selection.

For the twelve example states spanning the three panels, the maximum original-source residual is $8.68\times10^{-19}$ and the smallest strict branch margin is $2.38\times10^{-4}$. Interval-checked root brackets have width approximately $2\times10^{-45}$ and bound the error in the computed, unrounded double-precision permeabilities by $1.95\times10^{-16}$; the table adds decimal rounding. The grid itself uses ordinary floating-point bisection; the interval checks apply only to the example states. This addition establishes coupled branchwise uniqueness and coexistence. It does not infer PDE stability, clean-bed reachability, or experimental identification of these equilibria.
