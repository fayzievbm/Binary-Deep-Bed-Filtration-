# Analytical results within the unchanged saturation framework

These results concern the existing model with the multistage deposition closure.
They add no physical state, balance equation, porosity evolution, dispersion,
aggregation or external cake. The admissibility reserve introduced below is an
algebraic combination of the original variables, used for analysis only.

## Notation and standing assumptions

Write $\xi=x^*$, $\tau=t^*$, $u=U^*$, and $p=p^*$. The rate coefficients and sources
below are dimensionless; stars on them are suppressed for readability. Let

$$a_i=S_{2i},\quad B_i=S_{4i},\quad C_i=B_{i,\max},\quad
w_i=\frac1{1-\varepsilon_i},\quad \eta_i=w_i-1,$$

$$D=\sum_iw_iB_i,\qquad\phi=1-D,\qquad c_i=\frac{a_i}{\phi},\qquad
\kappa(D)=(1-\sqrt D)^3,\qquad\mu(c)=1+\frac52\sum_i c_i.$$

The retained equations are

$$m_0a_{i,\tau}+\partial_\xi(uc_i)=-q_i,\qquad
m_0B_{i,\tau}=q_i,\qquad u_\xi=0,\qquad
p_\xi=-u\frac{\mu}{\kappa},$$

with $\alpha_i=\alpha_{si}+\alpha_{di}u$ and
$\beta_i=\beta_{si}+\beta_{di}|p_\xi|$. The source branches and threshold
completion are those already stated in the manuscript: attractive charging
interfaces may hold, the prescribed convention selects a side at a repelling
aging interface, and the terminal capacity is absorbing. Parameters satisfy

$$0<B_{i1}<B_{i2}<C_i,\quad\chi_i>0,\quad0\le\varepsilon_i<1,\quad
D_{\max}=\sum_iw_iC_i<1.$$

All rate coefficients are nonnegative. The column is $0\le\xi\le1$, flow
is positive, and the same inlet concentration programme is used when comparing
operating modes. Statements about continuous solutions apply on their interval
of existence with sufficient regularity for the displayed balances, traces and
characteristic comparisons. They are not global existence or uniqueness
theorems for a discontinuous PDE system.

## Proposition 1 Operating modes are related by a positive clock change

Suppose $\alpha_{si}=\beta_{si}=0$ for both components. Define injected
initial pore volumes by

$$N(\tau)=\frac1{m_0}\int_0^\tau u(s)\,ds.$$

Because $u>0$, this clock is strictly increasing on its range. Put
$W=\mu/\kappa$. Darcy's law gives $|p_\xi|=uW$, and therefore

$$q_i=uF_i,\qquad
F_i=\begin{cases}
\chi_i\alpha_{di}c_i,&0\le B_i<B_{i1},\\
\alpha_{di}c_i-\beta_{di}WB_i,&B_{i1}<B_i\le B_{i2},\\
\alpha_{di}(C_i/B_i)c_i-\beta_{di}WB_i,&B_{i2}<B_i<C_i,\\
0,&B_i=C_i.
\end{cases}$$

The corresponding one-sided interface completion is inherited by $F_i$.
The state equations in the injected-volume clock are

$$a_{i,N}+\partial_\xi c_i=-F_i,\qquad B_{i,N}=F_i.$$

Consequently, any admissible state path solving these equations generates
corresponding constant-flow and constant-pressure solutions by changing its
clock. This assertion is a solution correspondence; it does not assert that
every admissible solution with given initial data is unique.

**Proof.** On each open kinetic branch, factor $u$ out of attachment and
detachment. Multiplication by a positive scalar preserves the signs and convex
weights of one-sided source values at the charging threshold. It therefore
preserves the zero-rate holding selection. It also preserves the selected side
at the aging threshold and the absorbing zero source at capacity. Next use
$m_0\partial_\tau=u\partial_N$, $u_\xi=0$, and divide the particle
balances by $u$. All remaining rate functions depend on the state and not on
the operating speed. Algebraic reconstruction gives the same seven saturation
fields along the state path. Reversing this calculation constructs solutions in
the original time variable. □

For such a common state path define its dimensionless hydraulic resistance

$$\mathcal R(N)=\int_0^1\frac{\mu(c(\xi,N))}{\kappa(D(\xi,N))}\,d\xi.$$

At constant prescribed flow $u_f>0$, and at constant pressure drop $P>0$,
the corresponding responses are

$$\Delta p_f(N)=u_f\mathcal R(N),\qquad
u_h(N)=\frac P{\mathcal R(N)},$$

$$\tau_f(N)=\frac{m_0N}{u_f},\qquad
\tau_h(N)=\frac{m_0}{P}\int_0^N\mathcal R(s)\,ds.$$

In particular,

$$u_h(N)\Delta p_f(N)=Pu_f.$$

Thus breakthrough responses, deposit profiles, saturation fields and resistance
coincide as functions of injected pore volume for corresponding solutions;
physical loading times and the controlled hydraulic observables differ. A
time-dependent inlet programme must be matched as a function of $N$, not
merely assigned the same frequency in $\tau$. Constant feeds satisfy this
requirement directly. The result holds on the common range of $N$ if one
protocol only reaches a finite injected volume.

For the general original coefficients, on every open branch one instead has

$$\frac{q_i}{u}=F_i+\frac{H_i}{u},$$

where $H_i$ is the same branch formula with $\alpha_{di}$ replaced by
$\alpha_{si}$, and $\beta_{di}W$ replaced by $\beta_{si}$. At a
generalized interface, the **same convex weight** must be used for these two
one-sided decompositions; separate completions of $F_i$ and $H_i$ need
not add to the completed full source. This identity isolates how static rates
can break the operating-mode collapse. Nonzero static coefficients do not
necessarily cause a difference on every special trajectory, but they introduce
the explicit factor $1/u$. Different initial states or inlet programmes can
also invalidate a proposed comparison.

## Proposition 2 An invariant region for the seven saturations

For the arrested multistage closure, consider initial data satisfying

$$a_i\ge0,\qquad0\le B_i\le C_i,\qquad
z:=1-\sum_iw_i(a_i+B_i)\ge0,$$

and inlet concentrations satisfying
$\sum_iw_ic_{i,\mathrm{in}}\le1$. Under the declared threshold/arrest
completion these inequalities are preserved. Consequently all seven reconstructed
saturations are nonnegative and $\phi\ge\phi_{\min}:=1-D_{\max}>0$.

**Proof.** At $a_i=0$, $c_i=0$, and every active branch satisfies
$q_i\le0$; charging, holding and capacity arrest give zero. Thus the mobile
particle equation is quasipositive. At $B_i=0$, charging gives
$q_i=\chi_i\alpha_ic_i\ge0$. Event crossing and absorbing arrest keep
$B_i$ in its prescribed interval. Summing the weighted particle inventories
gives

$$m_0\partial_\tau\sum_iw_i(a_i+B_i)+\partial_\xi(uY)=0,
\qquad Y:=\sum_iw_ic_i.$$

Because $z=\phi(1-Y)$ and $u_\xi=0$, this is equivalent to

$$m_0z_\tau+\partial_\xi\left(u\frac z\phi\right)=0.$$

This is a conservative transport equation for the derived reserve $z$, with
nonnegative initial and inflow data. Positivity follows along regular
characteristics; the same argument extends to admissible weak formulations for
which the corresponding positivity comparison is available. Equivalently,

$$m_0\phi Y_\tau+uY_\xi=(Y-1)\sum_iw_iq_i.$$

Along a common characteristic, $Z=1-Y$ obeys

$$\frac{dZ}{d\tau}=\frac{\sum_iw_iq_i}{m_0\phi}Z,$$

so its sign is preserved, independently of the sign of net deposition. Finally,

$$S_1=\phi-\sum_i a_i=z+\sum_i\eta_i a_i\ge0,$$

and $S_{3i}=\eta_iB_i\ge0$. The original saturation sum remains one. □

This region is a sufficient invariant subset of the physical saturation
simplex, not a characterization of every physically admissible initial state.
It contains the clean-bed initial condition and the dilute inlets used here.
If $\varepsilon_{\min}=\min_i\varepsilon_i>0$, then

$$\sum_i c_i\le(1-\varepsilon_{\min})Y\le1-\varepsilon_{\min},\qquad
S_1\ge\varepsilon_{\min}\phi_{\min}.$$

The proposition guarantees volume admissibility. It does **not** prove that a
dilute constitutive approximation remains accurate near the edge of this region.

### Corollary 2a A sufficient discrete invariant-domain condition

For uniform cells of width $h$, let the transport step use one positive Darcy
flux $u^n$. It preserves the region in Proposition 2 if

$$\lambda_j:=\frac{\Delta\tau\,u^n}{m_0h\phi_j^n}\le1.$$

Indeed,

$$a_{i,j}^{\mathrm{tr}}=(1-\lambda_j)a_{i,j}^n+
\frac{\Delta\tau u^n}{m_0h}\frac{a_{i,j-1}^n}{\phi_{j-1}^n},$$

$$z_j^{\mathrm{tr}}=(1-\lambda_j)z_j^n+
\frac{\Delta\tau u^n}{m_0h}\frac{z_{j-1}^n}{\phi_{j-1}^n}.$$

At the inlet the latter upwind value is
$1-\sum_iw_ic_{i,\mathrm{in}}\ge0$. Deposits and $\phi$ do not change
in this substep.

For the frozen-concentration event source map put

$$M_i=\max\left(1,\chi_i,\frac{C_i}{B_{i2}}\right).$$

Every branch satisfies $q_i\le\alpha_i^{\mathrm{src}}M_i c_i^{\mathrm{tr}}$.
Therefore its net deposit increment satisfies

$$\delta B_{i,j}\le\frac{\Delta\tau}{m_0}
\alpha_i^{\mathrm{src}}M_i\frac{a_{i,j}^{\mathrm{tr}}}{\phi_j^n}.$$

If $\alpha_i^{\mathrm{src}}\le\overline\alpha_i$ and

$$\frac{\Delta\tau\,\overline\alpha_iM_i}{m_0\phi_j^n}\le1,$$

then $a_{i,j}^{n+1}=a_{i,j}^{\mathrm{tr}}-\delta B_{i,j}\ge0$.
Paired exchange preserves $a_i+B_i$ and hence $z$ exactly. Event handling
preserves deposit bounds. The two substeps consequently preserve the same
convex region, without independent clipping of mobile/deposited inventories.

A globally safe bound for the newly supplied solver fork is

$$\overline\alpha_i=\alpha_{si}+\alpha_{di}\overline u,\qquad
\overline u=\begin{cases}u_f,&\text{prescribed constant flow},\\
P_{\max},&\text{positive prescribed head bounded by }P_{\max}.
\end{cases}$$

For pressure control this follows from $0<\kappa\le1$, $\mu\ge1$, and
$u=P/\int\mu/\kappa\le P_{\max}$. It remains safe when the source step
recomputes Darcy flow after transport. A condition based only on the earlier
attachment rate need not bound the recomputed dynamic rate. A safety factor
$0<\nu<1$ can be applied to both restrictions. This is an invariant-domain
result for the event-split discretization, not a convergence theorem or a claim
that the complete coupled method is exact or second order.

## Proposition 3 Two stable nonterminal local branches and four binary steady states

For one component prescribe constant positive local $c,\alpha,\beta$. Write
$r=\alpha c/\beta$, $b=B_1$, $T=B_2$, and $C=B_{\max}$.
If

$$\max(b,T^2/C)<r<T,$$

there are two locally asymptotically stable **nonterminal interior** source
equilibria,

$$B_L=r,\qquad B_H=\sqrt{Cr}.$$

Their local source relaxation derivatives are $-\beta/m_0$ and
$-2\beta/m_0$, respectively. The interface $B=T$ separates their open
basins because

$$q(T^-)=\beta(r-T)<0,\qquad
q(T^+)=\beta(Cr/T-T)>0.$$

**Proof.** The transition equilibrium lies in $(b,T)$ exactly when
$b<r<T$; its derivative is negative. The aging equilibrium lies in $(T,C)$
exactly when $T^2/C<r<C$; differentiating its source at the root gives
$-2\beta/m_0$. On either side of either equilibrium the scalar source points
toward that root. The two interface signs give separation of the open basins.
The value exactly on the repelling interface requires the stated tie convention;
no unrestricted Filippov uniqueness assertion is used. □

The broader interval $T^2/C<r<T$ also has a stable low branch when
$r\le b$, but that branch is the generalized charging-floor state $B=b$,
not an interior smooth equilibrium. The absorbing state $B=C$ is an
additional boundary equilibrium for every feed; it is deliberately not counted
among these two nonterminal branches.

The interior branch admissibility limits are

$$r_{\uparrow}=T,\qquad r_{\downarrow}=T^2/C,\qquad
\frac{r_{\uparrow}}{r_{\downarrow}}=\frac CT>1.$$

They predict a below-capacity hysteresis interval under slow forcing cycles.
They are not exact switching values for finite-rate ramps. For prescribed
$r(\tau)$, the transition tracking error $e=B-r$ obeys
$e'=-\beta e/m_0-r'$. In aging, $e_2=B^2-Cr$ obeys
$e_2'=-2\beta e_2/m_0-Cr'$. On an interval that remains in the relevant open
branch, put $\Delta=\tau-\tau_0$. The finite-time estimates are

$$|e(\tau)|\le |e(\tau_0)|e^{-\beta\Delta/m_0}
+\frac{m_0}{\beta}(1-e^{-\beta\Delta/m_0})\|r'\|_\infty,$$

$$|e_2(\tau)|\le |e_2(\tau_0)|e^{-2\beta\Delta/m_0}
+\frac{m_0C}{2\beta}(1-e^{-2\beta\Delta/m_0})\|r'\|_\infty,$$

where the norm is taken over $[\tau_0,\tau]$. These estimates support the
quasistatic branch-loss limits away from a
switch. Full coupled PDE stability and an exact finite-rate hysteresis loop do
not follow from these local estimates.

Now take constant original static rates, $\alpha_{di}=\beta_{di}=0$, and a
constant uniform inlet $c_i=\gamma_i$. If both
$r_i=\alpha_{si}\gamma_i/\beta_{si}$ satisfy the strict interior window,
each of the four combinations $B_i\in\{r_i,\sqrt{C_ir_i}\}$ gives an exact
homogeneous steady state of the full unchanged model, provided the inlet
admissibility condition holds. Set $\phi=1-\sum_iw_iB_i$,
$a_i=\phi\gamma_i$, $S_1=\phi(1-\sum_i\gamma_i)$ and
$S_{3i}=\eta_iB_i$. Then all sources vanish, all concentration gradients
vanish and Darcy's law determines the uniform pressure gradient. At pressure
drop $P$, $u=P\kappa/\mu(\gamma)$; at prescribed flow $u_f$,
$\Delta p=u_f\mu(\gamma)/\kappa$. The four combinations generally give
different hydraulic responses under the same final feed and forcing.

This proves **existence of four interior-branch combinations**, not that these
are the only steady states or that all are globally stable in the PDE. Additional
absorbing-capacity states and spatial branch mixtures may also exist. Histories
starting from the same clean bed must be demonstrated numerically if
reachability under particular inlet programmes is claimed.

For a below-capacity history experiment keep the preconditioning ratio below
$C$. In fact, under prescribed constant forcing, capacity is reached from
below in finite time only when $r>C$; $r=C$ gives asymptotic approach in
the aging branch. Once capacity has been reached, the imposed absorbing branch
does not release it during a later downward cycle. That terminal memory is
different from the reversible below-capacity branch coexistence analyzed here.

## Proposition 4 Matched uniform terminal deposits give matched terminal hydraulics

Suppose two kinetic closures retain the same balances and hydraulic relations,
both arrest at the same uniform capacities $C_i$, and reach a uniform final
mobile feed $\gamma_i$. Then their terminal saturation and hydraulic fields
are identical, regardless of different loading trajectories. With

$$D_C=\sum_iw_iC_i,\qquad\kappa_C=(1-\sqrt{D_C})^3,\qquad
\mu_C=1+\frac52\sum_i\gamma_i,$$

their terminal responses are

$$u_h=\frac{P\kappa_C}{\mu_C},\qquad
\Delta p_f=\frac{u_f\mu_C}{\kappa_C}.$$

**Proof.** The specified deposit and mobile concentration fields determine
$\phi$, all seven saturations, $\kappa$ and $\mu$ algebraically. Darcy's
law then gives the displayed responses. No kinetic coefficient remains. □

This statement is conditional on the terminal states actually being reached.
It does not assert that linear kinetics reach a prescribed cap for every
parameter set. It motivates comparison with a capped-linear baseline when the
purpose is to isolate trajectory changes from the effect of imposing a finite
capacity. The larger terminal flow relative to an uncapped linear baseline is
not, by itself, a distinct hydraulic constitutive effect.

## Checks and integration notes

`theory_checks.py` uses independent algebra and the stated constitutive laws;
it does not import or execute a filtration solver. It records the four
homogeneous states, their source residuals, the local branch-loss levels,
terminal hydraulic values, exact rational finite-volume reserve identities and
an independently prescribed positive-resistance clock identity. These checks
corroborate the formulas and bookkeeping. They are not a substitute for the
proofs, numerical refinement, or full coupled PDE stability analysis.

The supporting numerical solver is the new fork supplied by the benchmark
agent, with the globally safe attachment bound above; the previously delivered
solver remains unmodified. No first-ever claim is implied by these derivations.
Their distinction from prior multistage-filtration publications must be made
against the inspected literature.
