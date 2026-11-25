# Uniform equilibria with the original Darcy-dependent rates

This note extends the existing prescribed-rate/local branch calculation to **uniform steady states of the full retained 2017 balance and hydraulic system**. It changes no physical equation. A branchwise monotone scalar equation is a reduction of this particular constitutive model, not a new general root-finding method. It proves neither PDE stability nor transient reachability.

## Setting and admissibility

Use the manuscript's dimensionless variables and suppress stars. Write \(b_i=B_{i1}\), \(T_i=B_{i2}\), \(C_i=B_{i,\max}\), \(w_i=(1-\varepsilon_i)^{-1}\), and \(a_i=S_{2i}\). Assume

\[
0\le\varepsilon_i<1,\quad 0<b_i<T_i<C_i,\quad
D_{\max}=\sum_iw_iC_i<1,\quad
\gamma_i>0,\quad \sum_iw_i\gamma_i\le1.
\]

Here \(c_i=\gamma_i\) is the fixed uniform final mobile feed. Put \(\mu_\gamma=1+(5/2)\sum_i\gamma_i\). All original coefficients \(\alpha_{si},\alpha_{di},\beta_{si},\beta_{di}\) are nonnegative, and \((\beta_{si},\beta_{di})\ne(0,0)\) for each component. Prescribe either a positive uniform pressure drop \(P\) over the unit interval or a positive Darcy flux \(u_f\). Thus the original \(\beta_i\) is positive for \(\kappa>0\). If a component has zero attachment, it cannot have an equilibrium in either positive open branch considered here.

The selected open branches are L (the reversible transition branch, \(b_i<B_i<T_i\)) and H (the aging branch, \(T_i<B_i<C_i\)). Threshold equalities, generalized charging-floor states, absorbing-capacity states and nonuniform equilibria are outside this proposition. Strict inequalities avoid any dependence on a tie convention at the repelling interface.

## Proposition 1: general static/dynamic rates

For each branch assignment \(\sigma\in\{L,H\}^2\), there is at most one uniform equilibrium whose two deposits lie strictly inside their assigned branches.

For \(0<\kappa<1\), define

\[
r_i^P(\kappa)=
\frac{\gamma_i(\alpha_{si}+\alpha_{di}P\kappa/\mu_\gamma)}
{\beta_{si}+\beta_{di}P},\qquad
r_i^u(\kappa)=
\frac{\gamma_i(\alpha_{si}+\alpha_{di}u_f)\kappa}
{\beta_{si}\kappa+\beta_{di}u_f\mu_\gamma}.
\]

Use the expression corresponding to the prescribed control. Put

\[
\mathcal B_i^{\sigma_i}(\kappa)=
\begin{cases}
r_i(\kappa),&\sigma_i=L,\\
\sqrt{C_ir_i(\kappa)},&\sigma_i=H,
\end{cases}
\qquad
G_\sigma(\kappa)=\sum_iw_i\mathcal B_i^{\sigma_i}(\kappa)
-(1-\kappa^{1/3})^2.
\]

The equilibrium exists if and only if a zero of \(G_\sigma\) in \((0,1)\) satisfies all selected strict branch inequalities. In that case the zero and the uniform deposited state are unique for that assignment. Automatically \(\kappa>(1-\sqrt{D_{\max}})^3>0\).

**Proof.** Uniform steady balances require \(q_i=0\); the mobile equations are then satisfied by the prescribed constant concentrations. Uniform Darcy flow gives \(u=P\kappa/\mu_\gamma\) under pressure control, or \(|p_\xi|=u_f\mu_\gamma/\kappa\) under flow control. Substitution into the original coefficient definitions yields the displayed ratios. The zero-source equations are \(B_i=r_i\) on L and \(B_i^2=C_ir_i\) on H. The retained permeability law \(\kappa=(1-\sqrt{D})^3\), with \(D=\sum_iw_iB_i\), is equivalent to \(G_\sigma=0\).

Both ratios are nondecreasing:

\[
\frac{d r_i^P}{d\kappa}=
\frac{\gamma_i\alpha_{di}P}
{\mu_\gamma(\beta_{si}+\beta_{di}P)}\ge0,
\qquad
\frac{d r_i^u}{d\kappa}=
\frac{\gamma_i(\alpha_{si}+\alpha_{di}u_f)\beta_{di}u_f\mu_\gamma}
{(\beta_{si}\kappa+\beta_{di}u_f\mu_\gamma)^2}\ge0.
\]

Therefore each \(\mathcal B_i^{\sigma_i}\) is continuous and nondecreasing. The term \(-(1-\kappa^{1/3})^2\) is strictly increasing on \((0,1)\). Hence \(G_\sigma\) is strictly increasing and has at most one zero. An accepted root determines \(\phi=1-\sum_iw_iB_i>0\), \(a_i=\phi\gamma_i\), \(S_1=\phi(1-\sum_i\gamma_i)\) and \(S_{3i}=(w_i-1)B_i\). The assumed feed and capacity bounds make all reconstructed saturations nonnegative and ensure their sum is one. The source terms vanish and a linear pressure profile completes an exact steady solution of every retained balance. Conversely any selected uniform equilibrium must satisfy these same relations. This establishes the equivalence and branchwise uniqueness. It does not claim uniqueness among different branch assignments. \(\square\)

For general rates, a root is not guaranteed just by selecting L/H labels. In particular, \(G_\sigma(0+)\) need not equal \(-1\). With nonzero static attachment, the pressure-controlled low-permeability limit can have a positive deposited candidate. For fixed flow and \(\beta_{di}=0\), the ratio is constant in \(\kappa\), with continuous limit \(\gamma_i(\alpha_{si}+\alpha_{di}u_f)/\beta_{si}\). General existence must be checked using endpoint signs and branch admissibility, rather than imported from the purely dynamic special case.

A sharper formulation is available when desired. The set \(I_\sigma\subset(0,1)\) of permeabilities satisfying all strict branch inequalities is an open interval, possibly empty, because each ratio is nondecreasing. If \(I_\sigma=(\ell_\sigma,h_\sigma)\), an admissible root exists exactly when its continuous one-sided endpoint values satisfy \(G_\sigma(\ell_\sigma+)<0<G_\sigma(h_\sigma-)\). Equality at an endpoint is excluded from these open branches.

## Proposition 2: purely dynamic rates and an explicit existence test

Suppose \(\alpha_{si}=\beta_{si}=0\), \(\beta_{di}>0\), and \(\alpha_{di}\gamma_i>0\) for each component. Define

\[
d_i=\frac{\alpha_{di}\gamma_i}{\beta_{di}\mu_\gamma}>0.
\]

Both controls give \(r_i=d_i\kappa\). For every branch assignment the unconstrained scalar equation is

\[
\sum_{i:\sigma_i=L}w_id_i\kappa+
\sum_{i:\sigma_i=H}w_i\sqrt{C_id_i\kappa}
=(1-\kappa^{1/3})^2.
\]

It has exactly one candidate zero in \((0,1)\), since \(G_\sigma(0)=-1\), \(G_\sigma(1)>0\), and \(G_\sigma\) is strictly increasing. Define

\[
\ell_\sigma=\max_i
\begin{cases}b_i/d_i,&\sigma_i=L,\\ T_i^2/(C_id_i),&\sigma_i=H,\end{cases}
\qquad
h_\sigma=\min\left(1,\min_i
\begin{cases}T_i/d_i,&\sigma_i=L,\\ C_i/d_i,&\sigma_i=H.\end{cases}\right).
\]

An accepted open-branch equilibrium exists if and only if

\[
\ell_\sigma<h_\sigma,
\qquad G_\sigma(\ell_\sigma)<0<G_\sigma(h_\sigma).
\]

These statements follow by translating the branch bounds into bounds on \(\kappa\), then using strict monotonicity. First check \(\ell_\sigma<h_\sigma\); no endpoint evaluation outside \((0,1]\) is needed if this fails. If a \(d_i\) vanishes, the divided bounds are not used: that component has no positive open-branch root. If all \(d_i=0\), the only unconstrained endpoint root is \(\kappa=1\), outside this positive-deposit interior family.

The purely dynamic candidates are the same under both controls, in agreement with the existing pore-volume correspondence. This result adds the equilibrium branch test and its coupled permeability feedback; it does not add a second operating-mode-equivalence claim.

## Corollary: partial hydraulic ordering

At fixed coefficients and one control mode, suppose two accepted assignments differ only by replacing L with H for at least one component, without replacing any H with L. Then the new assignment has strictly smaller \(\kappa\). At the accepted original root, an L component satisfies \(0<r_i=B_i<T_i<C_i\), hence \(\sqrt{C_ir_i}>r_i\). The replacement strictly raises \(G\) at that root. The new strictly increasing \(G\), if it has an accepted root, must cross zero at a lower permeability. The same argument covers general original static/dynamic rates.

Consequently, when all four states are accepted,

\[
\kappa_{HH}<\min(\kappa_{LH},\kappa_{HL})
\le\max(\kappa_{LH},\kappa_{HL})<\kappa_{LL}.
\]

There is no universal ordering between LH and HL; symmetry can make them equal. Under prescribed pressure the lower permeability means lower flux; under prescribed flow it means a larger pressure drop. This is a comparison of different deposited inventories, not a resistance advantage at matched retained mass.

## Regime-map construction

The map fixes \(\gamma=(0.002,0.0005)\), \(\mu_\gamma=1.00625\), \(\beta_{di}=0.05\), porosities, thresholds and capacities. Its axes are \(d_1/T_1\) and \(d_2/T_2\), varying **only the original attachment-rate coefficients** \(\alpha_{di}=d_i\beta_{di}\mu_\gamma/\gamma_i\). Static coefficients are zero in the first panel and \(\alpha_{si}=\beta_{si}=0.1\) in the two mixed panels. The latter compare \(P=10\) with \(u_f=10\). These controls match the clean-water Darcy state; they need not have the same final flow, pressure or deposited state.

At each grid point each assignment is counted once if its unique scalar root is strictly inside both branches. The count ranges from zero to four **within this homogeneous interior family**. Zero does not mean that the full closure has no equilibrium: absorbing caps always supply additional boundary possibilities, and this proposition excludes threshold and spatially varying states.

The overlaid boundary curves are parameterized independently of the grid classification. Fix one \(B_i\) at a selected branch endpoint, vary the other deposit along its selected interval, and set

\[
\kappa=(1-\sqrt{\sum_iw_iB_i})^3,\qquad
\rho_i=\begin{cases}B_i,&L,\\B_i^2/C_i,&H.\end{cases}
\]

The corresponding map coordinates are obtained from

\[
d_i^P=\frac{\rho_i(\beta_{si}+\beta_{di}P)-\gamma_i\alpha_{si}}
{\beta_{di}P\kappa},
\qquad
d_i^u=\frac{\rho_i}{\kappa}
+\frac{\rho_i\beta_{si}-\gamma_i\alpha_{si}}
{\beta_{di}\mu_\gamma u_f}.
\]

Only nonnegative dynamic coefficients are plotted. These are analytic boundary parameterizations evaluated numerically; threshold endpoints themselves are not asserted to be open-branch equilibria.

## Numerical checks and their limits

The script computes ordinary double-precision bisection roots for the maps and records the actual final bracket widths. For the twelve displayed example states, it also uses 100-digit calculations to place brackets of width approximately \(2\times10^{-45}\), then independently evaluates \(G\) at both finite-decimal endpoints with **80-digit outward-rounded mpmath interval arithmetic**. Strict negative/positive interval signs and analytical monotonicity enclose the unique root. These interval checks apply only to the examples; the full grid is explicitly not interval-certified. The saved endpoint strings, sign checks and error bounds allow the distinction to be audited.

The reported source residuals are evaluated from the **original** coefficient formulas and the chosen L/H source at each double-precision equilibrium, not just from the transformed scalar equation. Reconstructed saturation sums, positive branch margins, constitutive permeability and Darcy relations are checked separately. None of these algebraic checks demonstrates transient solver convergence, PDE stability, or clean-start reachability for the new dynamic-rate equilibria.
