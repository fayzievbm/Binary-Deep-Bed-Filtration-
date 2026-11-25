"""Independent algebraic checks for the new propositions; no solver imports.

Run with standard Python3. The recorded values corroborate the analytical
formulas, not convergence or stability of a discretized filtration problem.
"""
from fractions import Fraction as F
from itertools import product
from math import sqrt,isclose
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent

def four_uniform_states():
    eps=(.5,.3);w=tuple(1/(1-e) for e in eps)
    b=(.001,.0002);T=(.004,.001);C=(.02,.004)
    gamma=(.002,.0005);alpha=(1.,1.);beta=(1.,1.);m0=.35
    r=tuple(a*c/d for a,c,d in zip(alpha,gamma,beta))
    lower=tuple(t*t/c for t,c in zip(T,C))
    assert all(max(bb,ll)<rr<tt for bb,ll,rr,tt in zip(b,lower,r,T))
    mu=1+2.5*sum(gamma)
    rows=[]
    for branches in product(('L','H'),repeat=2):
        B=tuple(rr if state=='L' else sqrt(cc*rr)
                for state,rr,cc in zip(branches,r,C))
        q=[];der=[]
        for i,state in enumerate(branches):
            if state=='L':
                assert b[i]<B[i]<T[i]
                q.append(alpha[i]*gamma[i]-beta[i]*B[i])
                der.append(-beta[i]/m0)
            else:
                assert T[i]<B[i]<C[i]
                q.append(alpha[i]*C[i]*gamma[i]/B[i]-beta[i]*B[i])
                der.append((-alpha[i]*C[i]*gamma[i]/B[i]**2-beta[i])/m0)
        D=sum(ww*bb for ww,bb in zip(w,B));phi=1-D;k=(1-sqrt(D))**3
        a=tuple(phi*c for c in gamma)
        reserve=1-sum(ww*(aa+bb) for ww,aa,bb in zip(w,a,B))
        carrier=phi-sum(a)
        u_head=k/mu;dp_flow=mu/k
        assert max(abs(v) for v in q)<1e-16
        assert carrier>0 and reserve>0
        assert isclose(u_head*dp_flow,1.,rel_tol=0,abs_tol=2e-15)
        rows.append(dict(branches=''.join(branches),B=B,D=D,phi=phi,
          mobile_saturations=a,carrier=carrier,reserve=reserve,
          k_relative=k,mu_relative=mu,u_at_unit_head=u_head,
          pressure_at_unit_flow=dp_flow,source_residual=q,
          local_source_relaxation_eigenvalues=der))
    return dict(eps=eps,B1=b,B2=T,capacity=C,inlet=gamma,r=r,
        lower_switch_limits=lower,upper_switch_limits=T,
        maximum_source_residual=max(abs(v) for x in rows for v in x['source_residual']),
        rows=rows,note='Exact homogeneous steady-state formulas. Local source eigenvalues are not coupled PDE stability eigenvalues.')

def rational_reserve_identities():
    # Exact rational arithmetic tests the conservative identities independently
    # of any numerical solver or event-source implementation.
    w=(F(2),F(10,7));eta=tuple(v-1 for v in w)
    m0=F(7,20);u=F(3,2);h=F(1,3)
    C=(F(1,50),F(1,250))
    B=[(F(1,1000),F(1,2000)),(F(3,1000),F(1,1000)),(F(5,1000),F(3,2000))]
    a=[(F(1,100),F(1,200)),(F(1,50),F(1,500)),(F(1,200),F(1,100))]
    phi=[1-sum(v*b for v,b in zip(w,bj)) for bj in B]
    z=[1-sum(v*(aa+bb) for v,aa,bb in zip(w,aj,bj)) for aj,bj in zip(a,B)]
    cin=(F(1,50),F(1,100));zin=1-sum(v*c for v,c in zip(w,cin))
    dt=m0*h*min(phi)/(2*u);lam=dt*u/(m0*h)
    atr=[];ztr=[];transport_residuals=[]
    for j in range(3):
        upstream=cin if j==0 else tuple(a[j-1][i]/phi[j-1] for i in range(2))
        new=tuple(a[j][i]-lam*(a[j][i]/phi[j]-upstream[i]) for i in range(2))
        zdef=1-sum(v*(aa+bb) for v,aa,bb in zip(w,new,B[j]))
        zflux=(1-lam/phi[j])*z[j]+lam*(zin if j==0 else z[j-1]/phi[j-1])
        transport_residuals.append(zdef-zflux)
        assert all(v>=0 for v in new) and zdef>=0
        atr.append(new);ztr.append(zdef)
    exchange_residuals=[];carrier_residuals=[];minimum=F(1)
    for j in range(3):
        # Mixed capture/release increments within component availability bounds.
        capture=lambda i:min(atr[j][i]/4,(C[i]-B[j][i])/2)
        delta=(capture(0),-B[j][1]/3) if j!=1 else (-B[j][0]/5,capture(1))
        anew=tuple(aa-dd for aa,dd in zip(atr[j],delta))
        bnew=tuple(bb+dd for bb,dd in zip(B[j],delta))
        pnew=1-sum(v*bb for v,bb in zip(w,bnew))
        znew=1-sum(v*(aa+bb) for v,aa,bb in zip(w,anew,bnew))
        carrier=pnew-sum(anew)
        exchange_residuals.append(znew-ztr[j])
        carrier_residuals.append(carrier-znew-sum(e*aa for e,aa in zip(eta,anew)))
        assert all(v>=0 for v in anew) and all(0<=bnew[i]<=C[i] for i in range(2))
        assert znew>=0 and carrier>=0
        minimum=min(minimum,*anew,*bnew,znew,carrier)
    assert all(v==0 for v in transport_residuals+exchange_residuals+carrier_residuals)
    return dict(arithmetic='exact rational',cells=3,maximum_transport_CFL=float(max(lam/v for v in phi)),
        transport_reserve_identity_residuals=[str(v) for v in transport_residuals],
        source_reserve_identity_residuals=[str(v) for v in exchange_residuals],
        reconstructed_carrier_identity_residuals=[str(v) for v in carrier_residuals],
        minimum_checked_nonnegative_quantity=float(minimum),
        note='Checks algebraic positivity bookkeeping with hand-chosen admissible conservative exchanges, not an event-source solver.')

def terminal_and_clock_checks():
    eps=(.5,.3);C=(.02,.004);gamma=(.05,.01)
    D=sum(c/(1-e) for c,e in zip(C,eps));k=(1-sqrt(D))**3
    mu=1+2.5*sum(gamma);P=10.;uf=10.;m0=.35
    terminal=dict(D=D,phi=1-D,k_relative=k,mu_relative=mu,
       flow_at_pressure10=P*k/mu,pressure_at_flow10=uf*mu/k,
       label='Common uniform terminal field for any two arresting closures with these capacities and feed.')
    # Independently prescribe a positive resistance path. The closed-form
    # pressure clock must obey the defining pore-volume derivative exactly.
    Nvals=(0.,.2,1.,3.)
    rows=[]
    for N in Nvals:
        resistance=1+.7*N+.2*N*N
        tau_head=m0/P*(N+.35*N*N+(.2/3)*N**3)
        dtau_dN=m0/P*resistance
        u_head=P/resistance;dp_flow=uf*resistance
        conversion=u_head*dtau_dN/m0
        reciprocity=u_head*dp_flow-P*uf
        assert isclose(conversion,1.,rel_tol=0,abs_tol=2e-15)
        assert abs(reciprocity)<2e-14
        rows.append(dict(N=N,resistance=resistance,tau_head=tau_head,
            tau_flow=m0*N/uf,pore_volume_derivative_identity=conversion,
            hydraulic_reciprocity_residual=reciprocity))
    return dict(terminal=terminal,synthetic_positive_resistance_clock=rows,
        note='The resistance polynomial is an independent clock identity check, not a simulated filtration history.')

def main():
    out=dict(purpose='Independent checks supporting the analytical propositions; no filtration-solver imports or runs.',
        uniform_steady_states=four_uniform_states(),
        exact_reserve_algebra=rational_reserve_identities(),
        terminal_and_clock=terminal_and_clock_checks(),
        limitations=['No proof of global PDE existence or uniqueness.','Local source stability is distinct from full coupled PDE stability.','No first-ever novelty claim is established by these checks.'])
    target=ROOT/'theory_checks_results.json'
    target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'output':str(target),'steady_state_combinations':4,
        'maximum_source_residual':out['uniform_steady_states']['maximum_source_residual'],
        'exact_reserve_residuals_zero':True,'terminal':out['terminal_and_clock']['terminal']},indent=2))

if __name__=='__main__':main()
