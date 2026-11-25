## Hydraulic compatibility from the unchanged permeability relation

Consider the original homogeneous reference permeability, a constant-cross-section one-dimensional bed, and admissible deposited saturations. With the normalized bed coordinate $x^*\in[0,1]$, write $\overline f=\int_0^1 f(x^*)\,dx^*$ and

$$D=\sum_i\frac{B_i}{1-\varepsilon_i},\qquad 0\le D<1,\qquad r(D)=\frac{k_0}{k}=(1-\sqrt D)^{-3},\qquad \mathcal K=\frac{k_{\mathrm{eff}}}{k_0}=\frac{1}{\overline{r(D)}}.$$

**Proposition.** The original permeability function satisfies the global pointwise support inequality

$$r(D)\ge1+AD,\qquad A=\min_{0<D<1}\frac{(1-\sqrt D)^{-3}-1}{D}=21.269750043077\ldots.$$

Consequently,

$$\mathcal K\le\left[1+A\overline D\right]^{-1}=\left[1+A\sum_i\frac{\overline B_i}{1-\varepsilon_i}\right]^{-1}.$$

For a single effective deposited population of fixed porosity $\varepsilon$, this becomes $\mathcal K\le[1+A\overline B/(1-\varepsilon)]^{-1}$. The result does not contain deposition-rate coefficients or assume a uniform retained profile.

**Proof.** Set $s=\sqrt D$. For $0<s<1$,

$$H(s)=\frac{(1-s)^{-3}-1}{s^2}=\frac{s^2-3s+3}{s(1-s)^3},\qquad H'(s)=\frac{2s^3-8s^2+12s-3}{s^2(1-s)^4}.$$

The numerator $P(s)=2s^3-8s^2+12s-3$ is strictly increasing on $(0,1)$: its derivative is positive there, with $P'(s)\ge2$ on the closed interval. Since $P(0)=-3$ and $P(1)=3$, it has exactly one root $s_*=0.308586020218263\ldots$. Moreover, $H$ diverges at both endpoints, so this root is the global minimizer. Thus $A=H(s_*)$, with contact loading $D_*=s_*^2=0.095225331874146\ldots$. The inequality follows for every $D>0$ and extends to $D=0$ by the exact equality $r(0)=1$. Integrating the pointwise inequality and taking reciprocals proves the effective-permeability bound. The coefficient $A$ is optimal for a linear support through $(D,r)=(0,1)$ because any larger coefficient fails at $D_*$. This optimality concerns the permeability function; it does not assert that a kinetic trajectory can realize every supporting profile.

### Nominal loading and hydraulic normalization

The existing constant-flow Alem observations give the nominal mean deposited-solid saturation $\overline B=83(0.25/2600)(0.6795)=0.0054229326923$. Here the feed concentration is expressed in $\mathrm{kg\,m^{-3}}$ and the particle density in $\mathrm{kg\,m^{-3}}$. For a mapped target ratio $0<\mathcal K_{\mathrm{nom}}\le1$, a necessary compatibility condition is

$$\overline B\le B_{\max}(\varepsilon;\mathcal K_{\mathrm{nom}})=\frac{1-\varepsilon}{A}\left(\mathcal K_{\mathrm{nom}}^{-1}-1\right).$$

Using the rounded series value $\mathcal K_{\mathrm{nom}}=0.94$ gives the following deterministic limits. The rounded ratio and the 83-pore-volume recovered mass are a nominal mapped pair; their exact run/endpoint correspondence is not independently established.

| Assumed deposit porosity | Upper ratio bound at nominal loading | Maximum loading not excluded at ratio 0.94 | Fraction of nominal loading |
|:---|---:|---:|---:|
| 0 | 0.896584 | 0.003000966 | 55.34% |
| 0.3 | 0.858533 | 0.002100676 | 38.74% |
| 0.5 | 0.812553 | 0.001500483 | 27.67% |
| 0.7 | 0.722292 | 0.000900290 | 16.60% |

These are necessary limits, not estimated corrections to retained mass. At $\varepsilon=0.5$, a loading above 27.67% of the nominal value cannot attain a mapped ratio 0.94 under this bound; even the compact-deposit limit permits only 55.34%. Passing the bound is not sufficient for agreement with the observed spatial profile, the selected kinetic law or the complete experimental record.

A deterministic normalization sensitivity can be stated without assigning observation errors. Let $\lambda_B$ multiply the nominal loading and let $\lambda_K$ multiply the nominal mapped ratio. The necessary condition is

$$\lambda_K\le\frac{1}{\mathcal K_{\mathrm{nom}}\left[1+A\lambda_B\overline B_{\mathrm{nom}}/(1-\varepsilon)\right]}.$$

At nominal loading ($\lambda_B=1$), the largest ratio multipliers not excluded by the bound are 0.953813 for $\varepsilon=0$ and 0.864418 for $\varepsilon=0.5$. They indicate how far a normalization interpretation would have to move the nominal pair to satisfy this necessary inequality. No multiplier is fitted, proposed as a correction, or assigned a probability. The comparison requires the same bed interval, loading endpoint, clean reference and fluid-property mapping. The reported initial conductivity variability of approximately 7% does not supply an uncertainty interval for the normalized ratio, whose numerator and denominator may be correlated. It therefore cannot be used as a confidence region or as evidence of statistical rejection. The original dilute viscosity correction is below 0.03% at the selected feed, but the distinction between conductivity and intrinsic permeability remains explicit.

**Suggested figure caption.** Necessary hydraulic compatibility under the unchanged 2017 permeability relation. (a) Upper effective intrinsic-permeability bounds for four assumed fixed deposit porosities and the nominal mapped Alem loading–ratio pair, with no invented measurement error bars. (b) Deterministic normalization boundaries for loading and ratio multipliers; points below a curve are not excluded by that necessary bound. The multipliers are sensitivity coordinates, not fitted corrections or confidence intervals. Exact pairing of the rounded ratio with the recovered-mass endpoint is not independently verified.

**Source integration note.** Cite the existing Makhmudov et al. (2017) reference for the permeability relation and the existing Alem et al. (2015) reference for the measurements. The proof and inversion are new mathematical consequences of that unchanged relation. The supplied displays are deliberately unnumbered for integration into the manuscript's equation sequence.
