"""Independent checks of the frozen kinetic map (no PDE sweeps).

Run: python source_checks.py. NumPy only beyond the solver's Matplotlib import.
The independent reference is fixed-step classical RK4 on the smooth scalar ODE,
not the SSP-RK3 method used by the production source map.
"""
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'numerics'))
import reference_solver as reference
import smooth_solver as solver
import smooth_source


def independent_rhs(B,c,p):
    lo=np.asarray(p.B1)[:,None]
    T=np.asarray(p.B2)[:,None]
    cap=np.asarray(p.cap)[:,None]
    delta=T*p.smoothing_fraction
    s=np.clip((B-T+delta)/(2*delta),0,1)
    w=3*s*s-2*s*s*s
    return (c*(1+w*(cap/B-1))-B)/p.m0


def fixed_rk4(B,c,dt,p,steps):
    b=B.copy();h=dt/steps
    for _ in range(steps):
        k1=independent_rhs(b,c,p)
        k2=independent_rhs(b+h*k1/2,c,p)
        k3=independent_rhs(b+h*k2/2,c,p)
        k4=independent_rhs(b+h*k3,c,p)
        b+=h*(k1+2*k2+2*k3+k4)/6
    return b


def main():
    rng=np.random.default_rng(20261125)
    p=solver.Params()
    B=np.asarray(p.cap)[:,None]*rng.uniform(0,.999,(2,31))
    B[:,0]=p.B1;B[:,1]=p.B2;B[:,2]=p.cap
    c=rng.uniform(0,.03,B.shape)
    alpha=np.ones_like(B);beta=np.ones_like(B)
    sharp=[]
    for dt in (.0001,.01,.2):
        a=solver.source_advance(B,c,alpha,beta,dt,p)
        b=reference.source_advance(B,c,alpha,beta,dt,p)
        error=float(abs(a[0]-b[0]).max())
        assert error==0 and a[1:]==b[1:]
        sharp.append({'dt':dt,'maximum_difference':error,'identical_event_counts':True})
    smooth=[]
    for f in (.01,.1,.4,.6):
        p=solver.Params(smoothing_fraction=f)
        T=np.asarray(p.B2)[:,None]
        delta=f*T
        initial=T+delta*np.linspace(-.9,.9,9)[None,:]
        cc=np.broadcast_to(np.array([.002,.0005])[:,None],initial.shape)
        got=solver.source_advance(initial,cc,np.ones_like(cc),np.ones_like(cc),.01,p)[0]
        ref=fixed_rk4(initial,cc,.01,p,4096)
        ref2=fixed_rk4(initial,cc,.01,p,8192)
        error=float(abs(got-ref2).max())
        ref_error=float(abs(ref-ref2).max())
        tight_p=replace(p,source_rtol=1e-11,source_atol=1e-16)
        tight=solver.source_advance(initial,cc,np.ones_like(cc),np.ones_like(cc),.01,tight_p)[0]
        tight_error=float(abs(tight-ref2).max())
        assert error<5e-10 and tight_error<1e-11 and ref_error<2e-12
        smooth.append({'relative_halfwidth':f,'dt':.01,'SSP_vs_independent_RK4_error':error,
                       'tight_source_vs_independent_RK4_error':tight_error,
                       'RK4_4096_vs8192_error':ref_error})
    crossing=[]
    initial=.9*np.asarray(solver.Params().B2)[:,None]
    cc=np.array([.006,.0015])[:,None]
    ones=np.ones_like(cc)
    sharp_cross=solver.source_advance(initial,cc,ones,ones,.1,solver.Params())[0]
    for f in (.1,.05,.02,.01):
        value=solver.source_advance(initial,cc,ones,ones,.1,solver.Params(smoothing_fraction=f))[0]
        crossing.append({'relative_halfwidth':f,'dt':.1,'maximum_difference_from_sharp':float(abs(value-sharp_cross).max())})
    assert all(crossing[j+1]['maximum_difference_from_sharp']<crossing[j]['maximum_difference_from_sharp'] for j in range(len(crossing)-1))
    positivity=[]
    for f in (.01,.4,.6):
        p=solver.Params(smoothing_fraction=f)
        B=np.asarray(p.cap)[:,None]*rng.uniform(0,.99,(2,40))
        D=np.sum(B/(1-np.asarray(p.eps)[:,None]),axis=0)
        phi=1-D
        c=rng.uniform(0,.02,B.shape)
        a=phi*c
        boost=np.maximum(1,np.maximum(np.asarray(p.charging),np.asarray(p.cap)/(np.asarray(p.B2)*(1-f))))
        dt=.8*np.min(p.m0*phi[None,:]/boost[:,None])
        newB,_,_=solver.source_advance(B,c,np.ones_like(c),np.ones_like(c),dt,p)
        newa=a-(newB-B)
        assert newB.min()>=0 and newa.min()>=-1e-16
        assert np.max(abs((newa+newB)-(a+B)))<1e-17
        positivity.append({'relative_halfwidth':f,'dt':float(dt),'minimum_mobile_after_source':float(newa.min()),
                           'minimum_deposit_after_source':float(newB.min()),
                           'paired_exchange_error':float(abs((newa+newB)-(a+B)).max())})
    out={'sharp_limit':sharp,'narrow_width_crossing_limit':crossing,'smooth_map_independent_reference':smooth,
         'source_donor_admissibility':positivity,
         'reference_sha256':hashlib.sha256((ROOT/'numerics'/'reference_solver.py').read_bytes()).hexdigest(),
         'interpretation':'Checks isolate the frozen source map; full PDE CFL refinement separately measures transport/source freezing and Lie splitting error.'}
    (ROOT/'source_checks_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
