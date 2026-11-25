"""Reproduce kinetic B2-smoothing robustness without changing 2017 physics.

Run: python run_smoothing.py
Reuse completed runs: python run_smoothing.py --reuse
NumPy and Matplotlib only. All paths are relative to this portable folder.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict, replace
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import numpy as np
from numpy.polynomial import Polynomial
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / 'numerics'))
import smooth_solver as solver
import smooth_source

FEED = np.array([.002, .0005])
HIGH = np.array([.006, .0015])
PRE_END = 2.
END = 10.
CASES = {'LL': (False, False), 'HL': (True, False),
         'LH': (False, True), 'HH': (True, True)}
COLORS = {'LL': '#2267A0', 'HL': '#D16A24', 'LH': '#22885B', 'HH': '#8A56A3'}
WIDTHS = (0., .01, .05, .1, .2, .4, .6)


def inlet_for(case):
    pre = np.where(CASES[case], HIGH, FEED)
    return lambda t: pre.copy() if t < PRE_END - 1e-12 else FEED.copy()


def params(f):
    return solver.Params(smoothing_fraction=f)


def interior_roots(p, species, c=None):
    """All local interior roots; the absorbing cap is listed separately.

    Multiplication by B gives a quartic in the normalized blend coordinate;
    NumPy polynomial roots enumerate every blend root, including unstable ones.
    """
    c = FEED[species] if c is None else c
    T = p.B2[species]
    C = p.cap[species]
    lo = p.B1[species]
    g = p.alpha_s[species] * c
    beta = p.beta_s[species]
    r = g / beta
    high = np.sqrt(C * r)
    if p.smoothing_fraction == 0:
        return [dict(B=r, derivative=-beta, classification='stable', branch='transition'),
                dict(B=T, derivative=None, classification='repelling_interface', branch='sharp_interface'),
                dict(B=high, derivative=-2*beta, classification='stable', branch='aging')]
    d = p.smoothing_fraction*T
    L, H = T-d, T+d
    found = []
    if lo < r <= L:
        found.append((r, 'transition'))
    polyB = Polynomial([L, 2*d])
    weight = Polynomial([0, 0, 3, -2])
    polynomial = g*polyB - beta*polyB*polyB + g*weight*(C-polyB)
    for z in polynomial.roots():
        if abs(z.imag) < 1e-9 and 0 < z.real < 1:
            found.append((L+2*d*z.real, 'blend'))
    if H <= high < C:
        found.append((high, 'aging'))
    result = []
    for B, branch in sorted(found):
        if result and abs(B-result[-1]['B']) < 1e-12*T:
            continue
        s = np.clip((B-L)/(2*d), 0, 1)
        w = s*s*(3-2*s)
        wp = 6*s*(1-s)/(2*d)
        derivative = g*(wp*(C/B-1)-w*C/B**2)-beta
        classification = 'stable' if derivative < -1e-8 else 'unstable' if derivative > 1e-8 else 'degenerate'
        result.append(dict(B=float(B), derivative=float(derivative), classification=classification, branch=branch))
    return result


def saddle_node(species):
    """Solve Q=Q_x=0 with Newton in dimensionless x=B/B2,f=delta/B2."""
    p = params(.1)
    R = FEED[species] / p.B2[species]
    K = p.cap[species] / p.B2[species]
    def residual(v):
        x, f = v
        s = (x-1+f)/(2*f)
        w = s*s*(3-2*s)
        wx = 6*s*(1-s)/(2*f)
        return np.array([R*(1+w*(K/x-1))-x,
                         R*(wx*(K/x-1)-w*K/x**2)-1])
    v = np.array([.53, .525])
    for _ in range(30):
        r = residual(v)
        if np.linalg.norm(r, np.inf) < 1e-13:
            break
        h = 1e-6
        J = np.column_stack([(residual(v+np.eye(2)[j]*h)-residual(v-np.eye(2)[j]*h))/(2*h) for j in range(2)])
        v -= np.linalg.solve(J, r)
    if np.linalg.norm(residual(v), np.inf) > 1e-10:
        raise RuntimeError('Saddle-node calculation failed')
    return {'species': species+1, 'B': float(v[0]*p.B2[species]),
            'relative_halfwidth': float(v[1]), 'dimensionless_residual': residual(v).tolist(),
            'interpretation': 'Local prescribed-final-concentration saddle node; not a coupled PDE bifurcation theorem.'}


def stem(f, n, case, suffix=''):
    return f'smooth_f{f:g}_n{n}_{case}'+suffix


def save_run(f, n, case, r, p, suffix='', elapsed=None):
    name = stem(f,n,case,suffix)
    feed_values = np.array([inlet_for(case)(t) for t in r['time']])
    table = np.column_stack([r['hist'], r['pv'], feed_values])
    header = ('t_star,c1_out,c2_out,u_star,delta_p_star,B1_mean,B2_mean,inventory1,inventory2,'
              'input1,input2,output1,output2,max_balance1,max_balance2,PV,c1_in,c2_in')
    np.savetxt(OUT/(name+'.csv'), table, delimiter=',', header=header, comments='')
    y = r['final']
    D,phi,c,k,u,dp,alpha,beta = solver.fields(y,p)
    q = smooth_source.net_source(y[1],c,alpha,beta,p)
    roots = [interior_roots(p,i) for i in range(2)]
    stable = [np.array([v['B'] for v in rr if v['classification']=='stable']) for rr in roots]
    distance = [float(np.min(abs(y[1,i,:,None]-stable[i][None,:]),axis=1).max()) for i in range(2)]
    selected = [[int(np.argmin(abs(B-stable[i]))) for B in y[1,i]] for i in range(2)]
    np.savez_compressed(OUT/(name+'_fields.npz'), final=y, x=r['x'], c=c, q=q,
                        occupancy=D, mobile_fraction=phi, time=r['time'],
                        retained=r['retained'], u=r['u'], outlet=r['outlet'], pv=r['pv'])
    metrics = {
        'relative_halfwidth':f,'n':n,'case':case,'suffix':suffix,'parameters':asdict(p),
        'preconditioning_feed':np.where(CASES[case],HIGH,FEED).tolist(),
        'final_mean_retained':np.mean(y[1],axis=1).tolist(),
        'final_retained_range':[[float(row.min()),float(row.max())] for row in y[1]],
        'maximum_distance_to_any_stable_local_root':distance,
        'selected_stable_root_indices_by_species':[sorted(set(row)) for row in selected],
        'maximum_absolute_final_net_source':np.max(abs(q),axis=1).tolist(),
        'maximum_final_concentration_error':np.max(abs(c-FEED[:,None]),axis=1).tolist(),
        'final_flow':float(u),'final_occupancy_range':[float(D.min()),float(D.max())],
        'maximum_species_balance_error':np.asarray(r['maxerr']).tolist(),
        'minimum_saturation':float(r['minstate']),'minimum_carrier_saturation':float(r['min_carrier']),
        'minimum_derived_z':float(r['minimum_derived_z']),
        'minimum_mobile_fraction':float(r['minimum_mobile_fraction']),
        'maximum_source_donor_ratio':float(r['maximum_source_donor_ratio']),
        'source_diagnostics':smooth_source.DIAGNOSTICS.copy(),
        'holding_events':int(r['holds']),'capacity_events':int(r['cap_hits']),
        'steps':int(r['steps']),'wall_seconds':elapsed,
    }
    if f==0:
        original=ROOT/'inputs'/'sharp_baseline'/f'history_{case}_n{n}_fields.npz'
        original_csv=ROOT/'inputs'/'sharp_baseline'/f'history_{case}_n{n}.csv'
        if original.exists():
            ref=np.load(original)
            csv=np.genfromtxt(original_csv,delimiter=',',names=True)
            metrics['sharp_baseline_comparison']={
                'maximum_final_state_difference':float(abs(y-ref['final']).max()),
                'maximum_recorded_flow_difference':float(abs(r['u']-csv['u_star']).max()),
                'maximum_recorded_retained_difference':float(abs(r['retained']-np.column_stack([csv['B1_mean'],csv['B2_mean']])).max())}
    (OUT/(name+'.json')).write_text(json.dumps(metrics,indent=2,allow_nan=False)+'\n')
    return metrics


def read_run(f,n,case,suffix=''):
    return np.load(OUT/(stem(f,n,case,suffix)+'_fields.npz'))


def comparisons(widths,grids):
    comparisons=[]
    if 80 in grids and 160 in grids:
        for f in widths:
            for case in CASES:
                a,b=read_run(f,80,case),read_run(f,160,case)
                projected=b['final'].reshape(2,2,80,2).mean(axis=-1)
                comparisons.append({'relative_halfwidth':f,'case':case,'type':'grid80_vs160',
                    'maximum_final_state_difference':float(abs(a['final']-projected).max()),
                    'maximum_recorded_flow_difference':float(abs(np.interp(b['time'],a['time'],a['u'])-b['u']).max()),
                    'maximum_recorded_retained_difference':float(max(abs(np.interp(b['time'],a['time'],a['retained'][:,i])-b['retained'][:,i]).max() for i in range(2)))})
    for n in grids:
        for f in widths:
            if f==0:continue
            for case in CASES:
                baseline=OUT/(stem(0,n,case)+'_fields.npz')
                target=OUT/(stem(f,n,case)+'_fields.npz')
                if baseline.exists() and target.exists():
                    a,b=np.load(baseline),np.load(target)
                    comparisons.append({'relative_halfwidth':f,'case':case,'n':n,'type':'sharp_vs_smoothed',
                        'maximum_recorded_flow_difference':float(abs(np.interp(b['time'],a['time'],a['u'])-b['u']).max()),
                        'maximum_recorded_retained_difference':float(max(abs(np.interp(b['time'],a['time'],a['retained'][:,i])-b['retained'][:,i]).max() for i in range(2)))})
    for f in (.4,.6):
        file=OUT/(stem(f,160,'HH','_refined')+'_fields.npz')
        if file.exists() and (OUT/(stem(f,160,'HH')+'_fields.npz')).exists():
            a,b=read_run(f,160,'HH'),read_run(f,160,'HH','_refined')
            comparisons.append({'relative_halfwidth':f,'case':'HH','type':'cfl.8_vs.4_and_tighter_source_tolerance',
                'maximum_final_state_difference':float(abs(a['final']-b['final']).max()),
                'maximum_recorded_flow_difference':float(abs(np.interp(b['time'],a['time'],a['u'])-b['u']).max()),
                'maximum_recorded_retained_difference':float(max(abs(np.interp(b['time'],a['time'],a['retained'][:,i])-b['retained'][:,i]).max() for i in range(2)))})
    return comparisons


def make_figures(widths):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':11,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'savefig.dpi':300})
    fig,axes=plt.subplots(2,2,figsize=(9.2,6.5),layout='constrained')
    fgrid=np.linspace(.001,.74,450)
    for i in range(2):
        ax=axes[0,i]
        for cls,color,ls in [('stable','#2267A0','-'),('unstable','#B35C29','--')]:
            xx=[];yy=[]
            for f in fgrid:
                roots=[r for r in interior_roots(params(float(f)),i) if r['classification']==cls]
                for r in roots:xx.append(f);yy.append(1e3*r['B'])
            ax.scatter(xx,yy,s=2,color=color,label=cls.capitalize())
        for f in widths:
            for case in CASES:
                path=OUT/(stem(f,160,case)+'_fields.npz')
                if path.exists():
                    r=np.load(path)
                    ax.scatter(f,1e3*r['final'][1,i].mean(),s=26,facecolors='none',edgecolors=COLORS[case],linewidths=1,zorder=4)
        saddle=saddle_node(i)
        critical=saddle['relative_halfwidth']
        ax.axvline(critical,color='#777777',ls=':',lw=1)
        ax.scatter([critical],[1e3*saddle['B']],color='#333333',s=15,zorder=5)
        ax.set_xlabel(r'Relative halfwidth $\delta_i/B_{2i}$')
        ax.set_ylabel(f'Deposit $10^3 B_{i+1}$')
        ax.set_title(f'({chr(97+i)}) Type {i+1}: local roots and final loading',loc='left')
        ax.grid(alpha=.15)
        if i==0:ax.legend(fontsize=9,loc='center right')
    for j,f in enumerate((.4,.6)):
        ax=axes[1,j]
        for case in CASES:
            path=OUT/(stem(f,160,case)+'_fields.npz')
            if not path.exists():continue
            r=np.load(path)
            ax.plot(r['time'],r['u'],color=COLORS[case],label=case,
                    ls={'LL':'-','HL':'--','LH':'-.','HH':':'}[case])
        ax.axvline(PRE_END,color='#777777',ls=':',lw=.8)
        ax.axvspan(0,PRE_END,color='#D8D8D8',alpha=.25,lw=0)
        ax.set_xlabel(r'Time $t^*$')
        ax.set_ylabel(r'Darcy flux $U^*$')
        ax.set_title(f'({chr(99+j)}) Flux histories, relative halfwidth {f:g}',loc='left')
        ax.grid(alpha=.15)
        ax.legend(ncol=2,fontsize=9,loc='lower right')
    for ext in ('png','pdf'):fig.savefig(OUT/f'figure_smoothing_robustness.{ext}')
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--reuse',action='store_true')
    parser.add_argument('--widths',nargs='+',type=float,default=list(WIDTHS))
    parser.add_argument('--grids',nargs='+',type=int,default=[80,160])
    parser.add_argument('--skip-refinement',action='store_true')
    args=parser.parse_args()
    start=time.perf_counter()
    metrics=[]
    for f in args.widths:
        p=params(f)
        for n in args.grids:
            for case in CASES:
                name=stem(f,n,case)
                if args.reuse and (OUT/(name+'.json')).exists():
                    m=json.loads((OUT/(name+'.json')).read_text())
                else:
                    tick=time.perf_counter()
                    smooth_source.reset_diagnostics()
                    r=solver.simulate(p,n=n,tend=END,inlet=inlet_for(case),record_dt=.025)
                    elapsed=time.perf_counter()-tick
                    # For a sharp run no smooth source is evaluated.
                    if not np.isfinite(smooth_source.DIAGNOSTICS['minimum_internal_deposit']):
                        smooth_source.DIAGNOSTICS['minimum_internal_deposit']=None
                    m=save_run(f,n,case,r,p,elapsed=elapsed)
                metrics.append(m)
                print(json.dumps({'f':f,'n':n,'case':case,'B':m['final_mean_retained'],
                                  'flow':m['final_flow'],'seconds':m['wall_seconds']}),flush=True)
    if not args.skip_refinement and 160 in args.grids:
        for f in (.4,.6):
            if f not in args.widths:continue
            name=stem(f,160,'HH','_refined')
            if args.reuse and (OUT/(name+'.json')).exists():
                m=json.loads((OUT/(name+'.json')).read_text())
            else:
                p=replace(params(f),source_rtol=1e-11,source_atol=1e-16)
                tick=time.perf_counter();smooth_source.reset_diagnostics()
                r=solver.simulate(p,n=160,tend=END,inlet=inlet_for('HH'),record_dt=.025,cfl=.4)
                m=save_run(f,160,'HH',r,p,'_refined',time.perf_counter()-tick)
            metrics.append(m)
            print(json.dumps({'refinement':name,'seconds':m['wall_seconds']}),flush=True)
    roots={str(f):[interior_roots(params(f),i) for i in range(2)] for f in args.widths}
    plateaus=[]
    for f in args.widths:
        for n in args.grids:
            rows=[m for m in metrics if m['relative_halfwidth']==f and m['n']==n and not m['suffix']]
            B=np.array([m['final_mean_retained'] for m in rows]);U=np.array([m['final_flow'] for m in rows])
            plateaus.append({'relative_halfwidth':f,'n':n,
                             'retained_spread_across_histories':np.ptp(B,axis=0).tolist(),
                             'flow_spread_across_histories':float(np.ptp(U))})
    output={'scope':'Kinetic sensitivity variant only; all 2017 balances and hydraulic relations unchanged.',
            'widths':args.widths,'grids':args.grids,'parameters_at_sharp_limit':asdict(params(0)),
            'source_method':'Analytic exterior branches; adaptive SSP-RK3 step doubling in blend, h*L<=0.5; frozen posttransport source coefficients as in sharp reference.',
            'local_roots':roots,'local_saddle_nodes':[saddle_node(i) for i in range(2)],
            'plateau_separation':plateaus,'runs':metrics,'refinement':comparisons(args.widths,args.grids),
            'reference_sha256':hashlib.sha256((ROOT/'numerics'/'reference_solver.py').read_bytes()).hexdigest(),
            'wall_seconds_this_invocation':time.perf_counter()-start}
    (OUT/'smoothing_results.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    import csv
    root_rows=[]
    for f in sorted(set(args.widths+np.linspace(.001,.74,450).tolist())):
        for i in range(2):
            for rr in interior_roots(params(f),i):
                root_rows.append([f,i+1,rr['B'],rr['classification'],rr['branch'],rr['derivative']])
    with (OUT/'local_root_map.csv').open('w') as file:
        writer=csv.writer(file);writer.writerow(['relative_halfwidth','species','B','classification','branch','dq_dB']);writer.writerows(root_rows)
    summary=[]
    for m in metrics:
        summary.append([m['relative_halfwidth'],m['n'],m['case'],m['suffix'],*m['final_mean_retained'],m['final_flow'],
                        max(m['maximum_species_balance_error']),max(m['maximum_absolute_final_net_source'])])
    with (OUT/'smoothing_summary.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['relative_halfwidth','n','case','suffix','B1_mean','B2_mean','flow','max_balance_error','max_final_q']);writer.writerows(summary)
    make_figures(args.widths)


if __name__=='__main__':main()
