"""Kinetics-independent nominal consistency checks for Alem2015 Table1.

These are model implications, not experimental confidence bounds. Measurement
interval/reference and uncertainty caveats are explicitly retained.
"""
from pathlib import Path
import json,csv
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parent
CIN=.25/2600;PV=83.
edges=np.array([0,.02,.07,.12,.17,.22,.27,.32,.37,.4])
fractions=np.array([9.86,16.94,11.31,8.34,6.08,5.85,3.73,3.69,2.15])/100
relative_width=np.diff(edges)/.4
B=CIN*PV*fractions/relative_width
assert abs(fractions.sum()-.6795)<1e-12
with (ROOT/'alem/data_retention_profile_table1_flow.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['study','source','mode','x_left_m','x_right_m','fraction_of_injected_mass_in_bin','B_mean_in_bin_per_original_pore_volume','method'])
    for a,b,v,z in zip(edges[:-1],edges[1:],fractions,B):w.writerow(['Alem2015','Table1','constant_flow',a,b,v,z,'transcribed printed Table1; converted with reported density/feed/PV/sectionwidth'])
s=brentq(lambda s:2*s**3-8*s*s+12*s-3,.1,.5)
h=((1-s)**-3-1)/(s*s)
Bbar=float(np.sum(B*relative_width))
profiles={}
for eps in [0.,.3,.5]:
    D=B/(1-eps);k=(1-np.sqrt(D))**3
    profiles[str(eps)]={'eps_assumed':eps,'piecewise_constant_profile_keff':float(1/np.sum(relative_width/k)),'mean_occupied_fraction':Bbar/(1-eps),'kinetics_independent_maximum_keff':1/(1+h*Bbar/(1-eps))}
data={'source':'Alem etal2015Table1 andFigure6; original2017permeability relation only','feed_volume_fraction':CIN,'pore_volumes':PV,'retained_mass_fraction':float(fractions.sum()),'mean_B_original_pore_volume':Bbar,'bound_constant_h_min':h,'bound_minimizer_sqrt_D':s,'bound_derivation':'For 0<=D<1, r(D)=(1-sqrt(D))^-3 >=1+h_min*D. Integrate resistance over the same full bed to obtain Keff/K0<=1/(1+h_min*mean(D)). h_min is the minimum of [r(D)-1]/D. Stationarity gives2s^3-8s^2+12s-3=0.','assumed_porosity_checks':profiles,'limitations':['Nominal consistency check; not a statistical rejection test.','Mass/pressure measurements must cover the same bed and same loading endpoint; published term overall is used but precise external pressure reference requires author clarification.','Reported initial hydraulicconductivity variability is about7%; individual ratio measurement uncertainties are not supplied.','Hydraulic conductivity and intrinsicpermeability ratios coincide only at matched fluidviscosity/density; diluteviscosity correction here is small but not silently omitted from simulations.','Finite-bin resistance estimate assumes B uniform within each measured bin; global bound does not require this assumption.','Reported drymass and turbidimetric mass balances have experimental errors; observations are not forced to close exactly.']}
(ROOT/'alem/results/hydraulic_consistency.json').write_text(json.dumps(data,indent=2))
print(json.dumps(data,indent=2))
