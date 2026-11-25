"""Independent wrapper-consistency check against existing 2017-framework solver."""
from pathlib import Path
import importlib.util,json,time
import numpy as np
from validation_solver import Params,simulate

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('original_solver',ROOT/'reference_solver.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
results=[]
for mode in ['rate','pressure']:
 for kinetics in ['linear','multistage']:
    p=Params(mode=mode,kinetics=kinetics)
    q=old.Params(mode=mode,kinetics=kinetics,alpha_s=(p.alpha_s,0.),alpha_d=(p.alpha_d,0.),beta_s=(p.beta_s,1.),beta_d=(p.beta_d,0.),eps=(p.eps,.3))
    r=old.simulate(q,n=40,tend=.2,inlet=lambda t:np.array([.05,0.]),record_dt=.002)
    fast=simulate(p,sample=r['time'],sample_axis='time',cin=.05,n=40)
    results.append(dict(mode=mode,kinetics=kinetics,
        maximum_final_a_error=float(abs(fast['a_profile'][-1]-r['final'][0,0]).max()),
        maximum_final_B_error=float(abs(fast['B_profile'][-1]-r['final'][1,0]).max()),
        maximum_flow_error=float(abs(fast['flow']-r['u']).max()),
        maximum_pressure_error=float(abs(fast['pressure']-r['dp']).max()),
        maximum_balance_error=float(fast['max_balance_error'].max())))
# Original velocity/pressure-gradient coefficients remain supported.
p=Params(mode='pressure',alpha_d=.02,beta_d=.001)
q=old.Params(alpha_s=(1.,0.),alpha_d=(.02,0.),beta_s=(1.,1.),beta_d=(.001,0.))
r=old.simulate(q,n=40,tend=.2,inlet=lambda t:np.array([.05,0.]),record_dt=.002)
f=simulate(p,sample=r['time'],sample_axis='time',cin=.05,n=40)
results.append(dict(mode='pressure_dynamic_coefficients',kinetics='multistage',maximum_final_a_error=float(abs(f['a_profile'][-1]-r['final'][0,0]).max()),maximum_final_B_error=float(abs(f['B_profile'][-1]-r['final'][1,0]).max()),maximum_flow_error=float(abs(f['flow']-r['u']).max()),maximum_pressure_error=float(abs(f['pressure']-r['dp']).max())))
# Observation-axis conversion is exact for constant unit imposed flow.
PV=np.array([0.,1.,10.,83.]);p=Params(m0=.37,eps=.5,mode='rate',forcing=1.)
start=time.perf_counter();f=simulate(p,sample=PV,sample_axis='pv',cin=.25/2600,n=80);elapsed=time.perf_counter()-start
check=dict(PV_time_conversion_error=float(abs(f['time']-.37*PV).max()),
    maximum_balance_error=float(f['max_balance_error'].max()),
    maximum_all_saturation_negativity=float(-min(0.,f['minimum_saturation'].min())),
    wall_seconds=elapsed,steps=int(f['steps'][-1]),note='Timing excludes one-time C++ compilation; no experimental observations are fitted.')
out=dict(wrapper_consistency=results,PV_axis_and_speed=check)
(ROOT/'validation_solver_benchmark.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
