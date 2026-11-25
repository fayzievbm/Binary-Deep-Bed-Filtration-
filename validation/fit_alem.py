"""Conditional kinetic calibration and held-out aggregate tests, Alem et al2015.

Digitized markers are observations; the supplied curves are never sampled.
All 2017 balances/hydraulics are unchanged. Only a one-component specialization
and the existing original-linear/multistage deposition choices are used.
"""
from pathlib import Path
from dataclasses import replace,asdict
import json,csv,time
import numpy as np
from scipy.optimize import differential_evolution,minimize
from validation_solver import Params,simulate

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'alem/results';OUT.mkdir(exist_ok=True)
CIN=.25/2600
def read(name):
    with (ROOT/'alem'/name).open() as f:return list(csv.DictReader(f))
con=read('data_breakthrough.csv');hyd=read('data_hydraulics.csv')
train=[r for r in con if r['mode']=='constant_flow']
x=np.array([float(r['pore_volumes']) for r in train]);y=np.array([float(r['value']) for r in train])

def model_from(v,kind,cap=.02,b1=.001,b2=.1):
    return Params(m0=.37,eps=.5,alpha_s=10**v[0],beta_s=10**v[1],
        charging=.5,B1=b1*CIN,B2=b2*CIN,cap=cap,mode='rate',forcing=1.,kinetics=kind)

def loss(v,kind,cap=.02,b1=.001,b2=.1,n=60):
    try:
        r=simulate(model_from(v,kind,cap,b1,b2),sample=x,cin=CIN,n=n)
        return float(np.mean((r['c_normalized'][:,0]-y)**2))
    except (RuntimeError,ValueError,FloatingPointError):return 1e4

def fit(kind,cap=.02,b1=.001,b2=.1,start=None,global_search=False):
    bounds=[(-3.,.5),(-7.,.3)] if kind=='multistage' else [(-.7,1.),(-7.,.3)]
    args=(kind,cap,b1,b2,60)
    if global_search:
        r=differential_evolution(loss,bounds,args=args,seed=1937,popsize=8,maxiter=35,tol=.002,polish=False,workers=1)
        start=r.x
    elif start is None:start=np.array([-1.5,-2.])
    r=minimize(loss,start,args=args,method='Nelder-Mead',bounds=bounds,options={'maxiter':130,'xatol':2e-4,'fatol':1e-8})
    # A short refinement at a finer mesh does not change the observation split.
    rr=minimize(loss,r.x,args=(kind,cap,b1,b2,120),method='Nelder-Mead',bounds=bounds,options={'maxiter':70,'xatol':3e-4,'fatol':2e-8})
    p=model_from(rr.x,kind,cap,b1,b2)
    print(kind,cap,b1,b2,rr.fun,asdict(p),flush=True)
    return p,rr.x,float(rr.fun)

def score(p,n=240):
    scores={};rows=[]
    for mode in ['constant_flow','constant_head']:
        q=replace(p,mode='rate' if mode=='constant_flow' else 'pressure')
        for target,source,col in [('C_over_C0',con,'c_normalized'),('K_over_K0',hyd,'keff')]:
            data=[d for d in source if d['mode']==mode]
            obsx=np.array([float(d['pore_volumes']) for d in data]);obs=np.array([float(d['value']) for d in data])
            pred=simulate(q,sample=obsx,cin=CIN,n=n)[col]
            if pred.ndim==2:pred=pred[:,0]
            error=pred-obs
            scores[mode+'_'+target]={'n':len(data),'RMSE':float(np.sqrt(np.mean(error**2))),'MAE':float(np.mean(abs(error))),'bias':float(error.mean()),'maximum_absolute_error':float(abs(error).max())}
            for d,v,e in zip(data,pred,error):rows.append(dict(d,model=p.kinetics,predicted=float(v),residual=float(e),role='calibration' if mode=='constant_flow' and target=='C_over_C0' else 'held_out'))
        r=simulate(q,sample=np.array([83.]),cin=CIN,n=n)
        scores[mode+'_end']={k:float(r[k][-1]) for k in ['time','pv','flow','pressure','Bmean','input','output','inventory','keff','max_balance_error','minimum_saturation','minimum_carrier','maximum_concentration','cap_events','hold_events']}
        scores[mode+'_end']['retained_fraction_of_input']=float(r['Bmean'][-1]/r['input'][-1])
    return scores,rows

def main():
    start=time.perf_counter();fits={};scans=[]
    for kind in ['linear','multistage']:
        p,z,cost=fit(kind,global_search=True)
        scores,rows=score(p);fits[kind]={'parameters':asdict(p),'fit_mean_squared_error_N120':cost,'scores':scores}
        with (OUT/(kind+'_observation_predictions.csv')).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        pv=np.r_[0.,np.linspace(.05,83,251)]
        for mode in ['rate','pressure']:
            r=simulate(replace(p,mode=mode),sample=pv,cin=CIN,n=240)
            np.savez_compressed(OUT/f'{kind}_{mode}_fields.npz',**{k:r[k] for k in ['time','pv','flow','pressure','Bmean','input','output','keff','hydraulic_conductance','c_normalized','B_profile','x','max_balance_error']})
        if kind=='multistage':best=p;bestz=z
    # Conditional scans change only explicitly unmeasured fixed kinetic choices.
    for cap,b1,b2 in [(.01,.001,.1),(.04,.001,.1),(.02,.0005,.05),(.02,.002,.2)]:
        guess=bestz.copy();guess[0]+=np.log10(.02/cap)
        p,z,cost=fit('multistage',cap,b1,b2,start=guess)
        ss,rows=score(p,n=120)
        scans.append({'cap_assumed':cap,'B1_over_c0_assumed':b1,'B2_over_c0_assumed':b2,'parameters':asdict(p),'fit_RMSE':float(np.sqrt(cost)),'head_concentration_RMSE':ss['constant_head_C_over_C0']['RMSE'],'rate_keff83':ss['constant_flow_end']['keff'],'head_keff83':ss['constant_head_end']['keff']})
    refinement=[]
    for n in [60,120,240,480]:
        ss,_=score(best,n=n)
        refinement.append({'n':n,'calibration_RMSE':ss['constant_flow_C_over_C0']['RMSE'],'head_concentration_RMSE':ss['constant_head_C_over_C0']['RMSE'],'rate_keff83':ss['constant_flow_end']['keff'],'head_keff83':ss['constant_head_end']['keff']})
    result={'scope':'Exploratory aggregate one-component limiting test; not binary or fullfourstage validation.','source_doi':'10.1007/s11242-014-0402-8','concentration_volume_fraction':CIN,'calibration':'Only0.25g/LconstantflowC/C0 observations; all available unambiguous poststart points.','held_out':'ConstantheadC/C0; all hydraulic observations; dry/excavated retained mass endpoints.','fixed_assumptions':{'deposit_porosity':.5,'charging':.5,'baseline_B1_over_c0':.001,'baseline_B2_over_c0':.1,'baseline_capacity':.02,'dynamic_coefficients':0.,'pressure_forcing':'Normalized to initial measured Darcyflux rather than simultaneously imposing inconsistenthead/K0/Q values.'},'fits':fits,'conditional_sensitivity':scans,'refinement':refinement,'elapsed_seconds':time.perf_counter()-start,'uncertainty_note':'Pixel bounds describe digitization only; observational experimental error bars unavailable. No confidence intervals or uniqueness of kineticparameters claimed.'}
    (OUT/'fit_results.json').write_text(json.dumps(result,indent=2))
    print('Finished',result['elapsed_seconds'],flush=True)

if __name__=='__main__':main()
