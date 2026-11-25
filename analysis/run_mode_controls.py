"""Existing-framework operating-mode and capacity-matched kinetic comparisons.

Every PDE run advances the original t* balances. Injected PV is integrated from
the numerical transport flux; source hydraulics are recomputed after transport.
No reduced-PV equation is used to force agreement between the two modes.
"""
from pathlib import Path
from dataclasses import replace,asdict
import sys,json,csv,time,hashlib,importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'numerics'))
from solver import Params,simulate,fields

OUT=ROOT/'mode_controls';OUT.mkdir(exist_ok=True)
LEVELS=np.unique(np.r_[np.linspace(.5,1.,101),.8,.9,.95])
all_results={}

def clean_json(x):
    if isinstance(x,dict):return {k:clean_json(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean_json(v) for v in x]
    if isinstance(x,np.ndarray):return clean_json(x.tolist())
    if isinstance(x,np.generic):return clean_json(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x

def summary(r):
    D,phi,c,k,u,dp,al,be=fields(r['final'],Params(**r['p']))
    cap=np.asarray(r['p']['cap'])[:,None]
    masks=np.isclose(r['final'][1],cap,rtol=0,atol=1e-11)
    quality=[]
    for theta in [.5,.8,.9,.95,1.]:
        j=int(np.argmin(abs(r['quality_levels']-theta)))
        good=float(r['quality_good_pv_after_initial_pv'][-1,j])
        quality.append(dict(illustrative_threshold=theta,
            good_effluent_PV_after_first_initial_PV=good,
            energy_star_per_good_effluent_PV=(float(r['energy_star'][-1])/good if good>0 else None)))
    return dict(parameters=r['p'],n=r['n'],steps=r['steps'],
        final_time_star=float(r['time'][-1]),final_PV=float(r['pv'][-1]),
        final_retained=r['retained'][-1].tolist(),
        final_flow_star=float(r['u'][-1]),final_pressure_drop_star=float(r['dp'][-1]),
        final_resistance=float(r['dp'][-1]/r['u'][-1]),
        final_fraction_of_cells_at_capacity=masks.mean(axis=1).tolist(),
        max_species_balance_error=r['maxerr'].tolist(),minimum_saturation=float(r['minstate']),
        minimum_carrier=float(r['min_carrier']),maximum_occupancy=float(r['maxD']),
        maximum_source_donor_ratio=float(r['maximum_source_donor_ratio']),
        minimum_derived_z=r.get('minimum_derived_z'),minimum_mobile_fraction=r.get('minimum_mobile_fraction'),
        energy_star=float(r['energy_star'][-1]),
        energy_per_injected_PV=float(r['energy_per_pore_scale'][-1]/r['pv'][-1]),
        max_normalized_outlet=np.max(r['outlet']/np.array([.05,.01]),axis=0).tolist(),
        min_normalized_outlet_after_first_PV=np.min(r['outlet'][r['pv']>=1]/np.array([.05,.01]),axis=0).tolist(),
        quality=quality)

def save(name,r):
    data=np.column_stack([r['pv'],r['hist'],r['energy_star'],r['good_pv'],r['good_pv_after_initial_pv']])
    header='PV,t_star,c1_out,c2_out,u_star,delta_p_star,B1_integral,B2_integral,inventory1,inventory2,input1,input2,output1,output2,max_balance1,max_balance2,energy_star,good_PV_theta05,good_PV_theta05_after_first_PV'
    np.savetxt(OUT/f'{name}.csv',data,delimiter=',',header=header,comments='')
    np.savez_compressed(OUT/f'{name}_fields.npz',PV=r['pv'],time=r['time'],
        state_history=r['state_history'],final=r['final'],x=r['x'],quality_levels=r['quality_levels'],
        quality_good_pv=r['quality_good_pv'],quality_good_pv_after_initial_pv=r['quality_good_pv_after_initial_pv'])
    all_results[name]=summary(r)

def run(name,p,n,pv_end):
    start=time.perf_counter()
    r=simulate(p,n=n,tend=10,pv_end=pv_end,record_pv=.02,
        inlet_axis='pv',store_profiles=True,quality_threshold=.5,quality_levels=LEVELS)
    save(name,r)
    print(name,'seconds',round(time.perf_counter()-start,3),
        'PV',r['pv'][-1],'caps',all_results[name]['final_fraction_of_cells_at_capacity'],flush=True)
    return r

def compare(p,q):
    assert np.allclose(p['pv'],q['pv'],rtol=0,atol=1e-9)
    diff=p['state_history']-q['state_history']
    cerror=(p['outlet']-q['outlet'])/np.array([.05,.01])
    product=(p['u']/p['u'][0])*(q['dp']/q['dp'][0])
    return dict(maximum_state_trajectory_difference=float(abs(diff).max()),
        mean_absolute_state_trajectory_difference=float(abs(diff).mean()),
        maximum_normalized_outlet_trajectory_difference=float(abs(cerror).max()),
        final_state_difference=float(abs(diff[-1]).max()),
        maximum_normalized_hydraulic_product_error=float(abs(product-1).max()),
        pressure_to_flow_duration_ratio=float(p['time'][-1]/q['time'][-1]),
        flow_to_pressure_energy_ratio=float(q['energy_star'][-1]/p['energy_star'][-1]),
        clock_energy_identity_defect=float(abs(p['time'][-1]/q['time'][-1]-q['energy_star'][-1]/p['energy_star'][-1])))

def reference_check():
    ref=ROOT/'numerics/reference_solver.py'
    if not ref.exists():return {'note':'Original reference path unavailable; archived execution record remains supplied.'}
    spec=importlib.util.spec_from_file_location('_unchanged_reference',ref)
    old=importlib.util.module_from_spec(spec);sys.modules[spec.name]=old;spec.loader.exec_module(old)
    records=[]
    for mode in ['pressure','rate']:
        for kinetic in ['linear','multistage']:
            r=old.simulate(old.Params(mode=mode,kinetics=kinetic),n=80,tend=.2,record_dt=.01)
            s=simulate(Params(mode=mode,kinetics=kinetic),n=80,tend=.2,record_dt=.01)
            error=float(abs(r['final']-s['final']).max())
            assert error<2e-14,(mode,kinetic,error)
            records.append(dict(mode=mode,kinetics=kinetic,maximum_final_state_error=error))
    return dict(original_solver_sha256=hashlib.sha256(ref.read_bytes()).hexdigest(),cases=records)

def plots(pure,sweep,controls,refinement):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,
        'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':350})
    fig,axs=plt.subplots(2,2,figsize=(10.4,7.1),layout='constrained')
    p,q=pure[320];pv=p['pv']
    for i,color in enumerate(['#2267A0','#D16A24']):
        axs[0,0].plot(pv,p['outlet'][:,i]/[.05,.01][i],color=color,label=f'Type {i+1}, pressure')
        axs[0,0].plot(pv[::10],q['outlet'][::10,i]/[.05,.01][i],ls='none',marker='o',ms=3,mfc='white',color=color,label=f'Type {i+1}, flow')
    axs[0,0].set(xlabel='Injected pore volumes',ylabel='Outlet / inlet concentration',title='(a) Pure dynamic rates, 320 cells')
    axs[0,0].legend(fontsize=9.5,loc='best')
    grids=np.array([r['n'] for r in refinement]);errors=np.array([r['comparison']['maximum_state_trajectory_difference'] for r in refinement])
    axs[0,1].loglog(1/grids,errors,'o-',color='#2267A0')
    axs[0,1].set(xlabel='Cell width',ylabel='Maximum state difference',title='(b) Independent time integrations')
    eta=np.array([r['eta'] for r in sweep]);derr=np.array([r['comparison']['maximum_state_trajectory_difference'] for r in sweep])
    axs[1,0].plot(eta,derr,'s-',color='#D16A24')
    axs[1,0].set(xlabel='Static-rate multiplier',ylabel='Maximum state difference',title='(c) Effect of static rates, 160 cells')
    axs[1,1].plot(pv,p['u']/p['u'][0],color='#2267A0',label='Flow ratio at fixed pressure')
    axs[1,1].plot(pv,(q['dp']/q['dp'][0])**-1,'--',color='#D16A24',label='Inverse pressure ratio at fixed flow')
    axs[1,1].set(xlabel='Injected pore volumes',ylabel='Normalized hydraulic response',title='(d) Reciprocal hydraulic relation')
    axs[1,1].legend(fontsize=9.5)
    for ax in axs.flat:ax.grid(alpha=.2)
    for ext in ['png','pdf']:fig.savefig(OUT/f'figure_operating_mode_invariance.{ext}')
    plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(10.4,7.1),layout='constrained')
    styles={'linear':('#777777',':','Original linear'),
        'capped_linear':('#2267A0','--','Capacity-matched linear'),
        'charging_linear':('#009E73','-.','Charging without aging'),
        'multistage':('#D16A24','-','Multistage')}
    for j,mode in enumerate(['pressure','rate']):
        for kin,(color,ls,label) in styles.items():
            r=controls[(mode,kin)]
            axs[0,j].plot(r['pv'],np.max(r['outlet']/np.array([.05,.01]),axis=1),color=color,ls=ls,label=label)
            axs[1,j].plot(r['quality_levels'],r['quality_good_pv_after_initial_pv'][-1],color=color,ls=ls,label=label)
        axs[0,j].set(xlabel='Injected pore volumes',ylabel='Maximum normalized outlet',
            title=f'({chr(97+j)}) Fixed '+('pressure' if mode=='pressure' else 'flow'))
        axs[0,j].axhline(.95,color='#CCCCCC',lw=.8);axs[0,j].set_ylim(0,1.1)
        axs[1,j].set(xlabel='Illustrative component threshold',ylabel='Admissible volume, PV',
            title=f'({chr(99+j)}) Quality volume after first PV')
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=2,fontsize=10.5)
    for ax in axs.flat:ax.grid(alpha=.2)
    for ext in ['png','pdf']:fig.savefig(OUT/f'figure_capacity_matched_controls.{ext}')
    plt.close(fig)

def redraw_saved_figures():
    """Render without recomputing any PDE trajectory or diagnostic."""
    report=json.loads((OUT/'results.json').read_text())
    def load(name):
        a=np.loadtxt(OUT/f'{name}.csv',delimiter=',',skiprows=1)
        with np.load(OUT/f'{name}_fields.npz') as f:
            return dict(pv=a[:,0],time=a[:,1],outlet=a[:,2:4],u=a[:,4],dp=a[:,5],
                quality_levels=f['quality_levels'].copy(),
                quality_good_pv_after_initial_pv=f['quality_good_pv_after_initial_pv'].copy())
    pure={320:(load('pure_dynamic_pressure_n320'),load('pure_dynamic_rate_n320'))}
    controls={(mode,kin):load(f'control_{mode}_{kin}') for mode in ['pressure','rate'] for kin in ['linear','capped_linear','charging_linear','multistage']}
    plots(pure,report['static_rate_sweep'],controls,report['pure_dynamic_refinement'])

def main():
    started=time.perf_counter();reference=reference_check();pure={};refinement=[]
    dynamic=Params(alpha_s=(0.,0.),beta_s=(0.,0.),alpha_d=(.1,.1),beta_d=(.05,.05))
    for n in [80,160,320]:
        p=run(f'pure_dynamic_pressure_n{n}',dynamic,n,6.)
        q=run(f'pure_dynamic_rate_n{n}',replace(dynamic,mode='rate'),n,6.)
        pure[n]=(p,q);refinement.append(dict(n=n,comparison=compare(p,q)))
    sweep=[]
    for eta in [0.,.1,.3,1.]:
        if eta==0:p,q=pure[160]
        else:
            pars=replace(dynamic,alpha_s=(eta,eta),beta_s=(eta,eta))
            p=run(f'static_eta{eta:g}_pressure',pars,160,6.)
            q=run(f'static_eta{eta:g}_rate',replace(pars,mode='rate'),160,6.)
        sweep.append(dict(eta=eta,comparison=compare(p,q)))
    controls={}
    for mode in ['pressure','rate']:
        for kin in ['linear','capped_linear','charging_linear','multistage']:
            controls[(mode,kin)]=run(f'control_{mode}_{kin}',Params(mode=mode,kinetics=kin),160,14.)
    plots(pure,sweep,controls,refinement)
    results=dict(scope='Only kinetic closure controls and numerical diagnostics; 2017 physical balances and hydraulic relations unchanged.',
        common_initial_forcing='Delta p*=10 and U*=10 give the same clean-water initial state; mu*=1 initially, constant feed starts at t*=0.',
        integration='Independent original t* integrations; PV tracked from transport flux, source hydraulics recomputed after transport.',
        quality_definition='Both component outlet/inlet ratios must be <= illustrative theta. Curves show effluent PV excluding first initial PV; thresholds are not water-quality standards. Energy is total pump work through the common endpoint; dividing by zero admissible volume is undefined.',
        mode_test_endpoint_PV=6.,control_endpoint_PV=14.,reference_solver_check=reference,
        pure_dynamic_refinement=refinement,static_rate_sweep=sweep,cases=all_results,
        wall_seconds=time.perf_counter()-started)
    (OUT/'results.json').write_text(json.dumps(clean_json(results),indent=2))
    with (OUT/'quality_threshold_curves.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['mode','kinetics','illustrative_threshold','good_effluent_PV_after_first_PV','energy_star_per_good_effluent_PV'])
        for (mode,kin),r in controls.items():
            for theta,good in zip(r['quality_levels'],r['quality_good_pv_after_initial_pv'][-1]):w.writerow([mode,kin,theta,good,r['energy_star'][-1]/good if good>0 else 'NA'])
    print(json.dumps(clean_json({'refinement':refinement,'static_sweep':sweep,'elapsed':results['wall_seconds']}),indent=2),flush=True)

if __name__=='__main__':main()
