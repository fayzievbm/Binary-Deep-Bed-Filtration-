"""Accelerated scalar specialization of the unchanged 2017 saturation framework.

Only net deposition kinetics can be linear or the existing multistage law.
No experimental fit is performed in this module. A C++ source kernel is built
on demand using g++; NumPy and SciPy-compatible arrays are returned.
"""
from dataclasses import dataclass,asdict
from pathlib import Path
import ctypes,subprocess
import numpy as np

ROOT=Path(__file__).resolve().parent
COLS=('time','pv','flow','pressure','Bmean','inventory','input','output','max_balance_error',
      'minimum_saturation','minimum_carrier','maximum_concentration','keff','hydraulic_conductance','steps','hold_events','cap_events')

@dataclass
class Params:
    m0:float=.35
    eps:float=.5
    alpha_s:float=1.
    alpha_d:float=0.
    beta_s:float=1.
    beta_d:float=0.
    charging:float=.5
    B1:float=.001
    B2:float=.004
    cap:float=.02
    mode:str='pressure'
    forcing:float=10.
    kinetics:str='multistage'

_LIB=None
def _kernel():
    global _LIB
    if _LIB is None:
        src=ROOT/'validation_kernel.cpp';lib=ROOT/'validation_kernel.so'
        if not lib.exists() or src.stat().st_mtime>lib.stat().st_mtime:
            subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(src),'-o',str(lib)],check=True,capture_output=True)
        _LIB=ctypes.CDLL(str(lib));fp=_LIB.filtration_run
        arr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        fp.argtypes=[arr,ctypes.c_int,arr,ctypes.c_int,ctypes.c_int,ctypes.c_double,arr,arr,arr];fp.restype=ctypes.c_int
    return _LIB.filtration_run

def simulate(params=None,sample=None,sample_axis='pv',cin=.05,depths=(1.,),n=80,cfl=.8,dtmax=np.inf):
    """Return observables sampled at original t* or cumulative injected PV.

    Args:
      params: Params with fixed m0 and fixed deposit porosity eps.
      sample: nonnegative, strictly increasing 1D observation times/PV values.
      sample_axis: 'time' (t*=U0 t/L) or 'pv' (integral u*/m0 dt*).
      cin: inlet particle VOLUME fraction; convert mass concentration using
           independently justified particle density before calling.
      depths: measurement positions x/L in [0,1]. Inlet value is prescribed;
              internal values are linearly interpolated between cell centers;
              downstream value is the final upwind cell concentration.

    All inventory quantities refer to original pore volume. keff is
    1/mean(1/krel); hydraulic_conductance is 1/mean(mu_rel/krel).
    The latter includes the original suspension viscosity; they are distinct.
    """
    p=params or Params();sample=np.ascontiguousarray(np.asarray(sample,dtype=float))
    depths=np.asarray(depths,dtype=float)
    if sample.ndim!=1 or not len(sample) or np.any(sample<0) or np.any(np.diff(sample)<=0):raise ValueError('sample must be nonnegative and strictly increasing')
    if sample_axis not in ('pv','time'):raise ValueError('sample_axis must be pv or time')
    if n<2 or not (0<p.m0<1 and 0<=p.eps<1 and 0<cin<1 and p.forcing>0):raise ValueError('invalid geometry/forcing/concentration')
    if np.any((depths<0)|(depths>1)):raise ValueError('depths must lie in [0,1]')
    if min(p.alpha_s,p.alpha_d,p.beta_s,p.beta_d)<0:raise ValueError('kinetic coefficients must be nonnegative')
    if p.kinetics not in ('linear','multistage') or p.mode not in ('pressure','rate'):raise ValueError('invalid kinetics or hydraulic mode')
    if p.kinetics=='multistage' and not (0<p.B1<p.B2<p.cap<1-p.eps and p.charging>0):raise ValueError('require 0<B1<B2<cap<1-eps and positive charging')
    par=np.ascontiguousarray([p.m0,p.eps,p.alpha_s,p.alpha_d,p.beta_s,p.beta_d,p.charging,p.B1,p.B2,p.cap,p.forcing,float(p.mode=='rate'),float(p.kinetics=='multistage'),cfl,min(dtmax,1e100)],dtype=float)
    out=np.empty((len(sample),len(COLS)));a=np.empty((len(sample),n));B=np.empty_like(a)
    status=_kernel()(par,n,sample,len(sample),int(sample_axis=='pv'),float(cin),out,a,B)
    if status:raise RuntimeError(f'integration failed ({status}); parameters may violate admissibility')
    result={k:out[:,j] for j,k in enumerate(COLS)}
    phi=1-B/(1-p.eps);c=a/phi;x=(np.arange(n)+.5)/n
    krel=(1-np.sqrt(B/(1-p.eps)))**3;mu=1+2.5*c
    cell_resistance=mu/krel/n
    pressure_profile=out[:,3,None]-out[:,2,None]*(np.cumsum(cell_resistance,axis=1)-.5*cell_resistance)
    cdepth=np.array([np.interp(depths,np.r_[0,x,1],np.r_[cin,row,row[-1]]) for row in c])
    result.update(c_depth=cdepth,c_normalized=cdepth/cin,a_profile=a,B_profile=B,c_profile=c,
        phi_profile=phi,S1_profile=phi-a,S3_profile=p.eps/(1-p.eps)*B,
        k_profile=krel,mu_profile=mu,pressure_profile=pressure_profile,x=x,depths=depths,
        parameters=asdict(p),sample_axis=sample_axis,cin=cin)
    return result

def fitting_parameter_notes():
    """Return the existing physical parameters; this does not choose a fit."""
    return dict(fixed_from_experiment=['m0','mode','forcing','cin','column_length','depths'],
        fixed_structure_or_sensitivity=['eps'],kinetic_parameters=['alpha_s','beta_s','charging','B1','B2','cap'],
        optional_original_dynamic_coefficients=['alpha_d','beta_d'],
        recommended_initial_fit='Fit only 2–4 identifiable kinetic parameters; fix dynamic coefficients at zero unless protocol independently supports them. Do not fit porosity and cap simultaneously without independent structural data.',
        warning='Fitting scalar specialization does not validate mixed-component interactions. No dispersion or new mechanism is added.')

if __name__=='__main__':
    r=simulate(sample=np.array([0.,1.,5.,10.]),sample_axis='pv',n=80)
    print({k:r[k].tolist() for k in ['time','pv','flow','pressure','Bmean','max_balance_error']})
