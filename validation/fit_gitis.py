"""Conditional normalized-shape fit to digitized Gitis et al. (2010), Figure 6.

Only 0.82-m symbols are calibrated. Deeper Figure 6 symbols are spatial holdouts.
Figure 7 is only an unverified-condition external diagnostic. Unknown inlet
volume fraction, deposit porosity and time origin remain explicit assumptions.
"""
from pathlib import Path
import json,csv,time
import numpy as np
from scipy.optimize import least_squares,differential_evolution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validation_solver import Params,simulate

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'gitis';DEPTHS=np.array([.82,1.02,1.62]);L=1.62;U=10.;M0=.44

def read_points(fig=6):
    with (OUT/f'figure{fig}_experimental_points.csv').open() as f:rows=list(csv.DictReader(f))
    return rows

def makeparams(theta,cin=1e-6,kind='multistage'):
    if kind=='linear':
        return Params(m0=M0,eps=.5,alpha_s=float(theta[0]),beta_s=float(theta[1]),mode='rate',forcing=1.,kinetics='linear')
    alpha,chi,r=np.exp(theta)
    return Params(m0=M0,eps=.5,alpha_s=float(alpha),beta_s=0.,charging=float(chi),B1=float(r*cin),B2=1e4*cin,cap=2e4*cin,mode='rate',forcing=1.,kinetics='multistage')

def predict(rows,theta,cin=1e-6,offset=.1,n=80,kind='multistage'):
    times=np.array([float(r['time_h']) for r in rows]);unique,idx=np.unique(times,return_inverse=True)
    p=makeparams(theta,cin,kind)
    result=simulate(p,sample=U/L*(unique+offset),sample_axis='time',cin=cin,depths=DEPTHS/L,n=n)
    depth_idx=np.array([int(np.argmin(abs(DEPTHS-float(r['depth_m'])))) for r in rows])
    pred=result['c_normalized'][idx,depth_idx]
    return pred,result

def metrics(rows,pred):
    obs=np.array([float(r['concentration_ratio']) for r in rows]);err=pred-obs
    return dict(n=len(obs),RMSE=float(np.sqrt(np.mean(err**2))),MAE=float(abs(err).mean()),bias=float(err.mean()),max_absolute_error=float(abs(err).max()),R2=float(1-np.sum(err**2)/np.sum((obs-obs.mean())**2)))

def fit(rows,cin=1e-6,offset=.1,n=80,kind='multistage',start=None):
    training=[r for r in rows if abs(float(r['depth_m'])-.82)<1e-5];obs=np.array([float(r['concentration_ratio']) for r in training])
    starts=[np.log([8.,.45,15.]),np.log([8.,.45,40.]),np.log([12.,.3,60.]),np.log([6.,.6,10.])] if start is None else [np.asarray(start)]
    bounds=(np.log([.1,.01,.01]),np.log([100.,1.,1000.]))
    if kind=='linear':starts=[np.array([5.,0.]),np.array([5.,.01])];bounds=([.01,0.],[100.,10.])
    best=None
    for x0 in starts:
        def residual(theta):
            try:return predict(training,theta,cin,offset,n,kind)[0]-obs
            except (RuntimeError,ValueError):return np.full(len(obs),10.)
        opt=least_squares(residual,x0,bounds=bounds,diff_step=1e-3,ftol=1e-9,xtol=1e-8,gtol=1e-9,max_nfev=180)
        if best is None or sum(opt.fun**2)<sum(best.fun**2):best=opt
    pred,result=predict(rows,best.x,cin,offset,n,kind)
    bydepth={str(d):metrics([r for r in rows if abs(float(r['depth_m'])-d)<1e-5],pred[np.array([abs(float(r['depth_m'])-d)<1e-5 for r in rows])]) for d in DEPTHS}
    held=np.array([abs(float(r['depth_m'])-.82)>1e-5 for r in rows]);trainingmask=~held
    p=makeparams(best.x,cin,kind)
    notes=dict(kinetics=kind,n=n,cin_assumed=cin,eps_assumed=.5,clock_offset_h_assumed=offset,
        alpha_s=p.alpha_s,beta_s=p.beta_s,charging=p.charging if kind!='linear' else None,
        B1_over_cin=p.B1/cin if kind!='linear' else None,B2_over_cin=p.B2/cin if kind!='linear' else None,cap_over_cin=p.cap/cin if kind!='linear' else None,
        training=metrics([r for r,k in zip(rows,trainingmask) if k],pred[trainingmask]),
        spatial_holdout=metrics([r for r,k in zip(rows,held) if k],pred[held]),by_depth=bydepth,
        objective_SSE=float(sum(best.fun**2)),optimizer_success=bool(best.success),optimizer_nfev=int(best.nfev),
        max_balance_error=float(result['max_balance_error'].max()),minimum_carrier=float(result['minimum_carrier'].min()),
        maximum_B=float(result['B_profile'].max()),aging_threshold=float(p.B2),aging_threshold_reached=bool(result['B_profile'].max()>=p.B2),
        cap_events=int(result['cap_events'][-1]),hold_events=int(result['hold_events'][-1]),
        fitted_fields=['alpha_s','beta_s'] if kind=='linear' else ['alpha_s','charging','B1_over_cin'])
    return best.x,pred,result,notes

def save_predictions(filename,rows,pred,role='figure6'):
    fields=list(rows[0])+['comparison_role','prediction_ratio','residual_prediction_minus_observation']
    with (OUT/filename).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for row,y in zip(rows,pred):
            role_i=('calibration' if abs(float(row['depth_m'])-.82)<1e-5 else 'spatial_holdout_same_experiment') if role=='figure6' else 'external_unknown_condition_diagnostic'
            w.writerow(dict(row,comparison_role=role_i,prediction_ratio=y,residual_prediction_minus_observation=y-float(row['concentration_ratio'])))

def plots(rows,theta,baseline,metrics_all):
    plt.rcParams.update({'font.size':12,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
    times=np.linspace(0,6.8,300);sample=U/L*(times+.1)
    curves={}
    for name,th,kind in [('Multistage',theta,'multistage'),('Linear',baseline,'linear')]:
        curves[name]=simulate(makeparams(th,kind=kind),sample=sample,sample_axis='time',cin=1e-6,depths=DEPTHS/L,n=160)['c_normalized']
    with (OUT/'figure6_plot_curves.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['published_time_h','depth_m','multistage_ratio','linear_ratio'])
        for j,d in enumerate(DEPTHS):
            for t,y,z in zip(times,curves['Multistage'][:,j],curves['Linear'][:,j]):w.writerow([t,d,y,z])
    fig,axs=plt.subplots(1,3,figsize=(9,3.6),layout='constrained')
    for j,(ax,d) in enumerate(zip(axs,DEPTHS)):
        data=[r for r in rows if abs(float(r['depth_m'])-d)<1e-5]
        ax.errorbar([float(r['time_h']) for r in data],[float(r['concentration_ratio']) for r in data],xerr=[float(r['time_digitization_envelope_h']) for r in data],yerr=[float(r['ratio_digitization_envelope']) for r in data],fmt='o',color='black',ms=4,capsize=2,label='Digitized symbols')
        ax.plot(times,curves['Multistage'][:,j],color='#D16A24',label='Multistage')
        ax.plot(times,curves['Linear'][:,j],color='#2267A0',ls='--',label='Linear')
        ax.set_xlabel('Published time (h)');ax.set_ylabel('Normalized response');ax.set_title(f'{d:.2f} m'+('\nCalibration' if j==0 else '\nSpatial holdout'),fontsize=11);ax.set_ylim(0,.22);ax.grid(alpha=.15)
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=10)
    fig.savefig(OUT/'gitis_figure6_conditional_fit.png');fig.savefig(OUT/'gitis_figure6_conditional_fit.pdf');plt.close(fig)

def main():
    rows=read_points();print('fitting main multistage',flush=True)
    coarse,cp,cr,coarseinfo=fit(rows,n=80)
    print('global shallow-only search',flush=True)
    training=rows[:9];observed=np.array([float(v['concentration_ratio']) for v in training])
    def objective(theta):
        try:return float(sum((predict(training,theta,n=80)[0]-observed)**2))
        except (RuntimeError,ValueError):return 100.
    glob=differential_evolution(objective,list(zip(np.log([.1,.01,.01]),np.log([100.,1.,1000.]))),
        popsize=8,maxiter=45,seed=20261125,polish=False,x0=coarse,tol=1e-7,workers=1)
    start=glob.x if glob.fun<coarseinfo['objective_SSE'] else coarse
    coarse,cp,cr,coarseinfo2=fit(rows,n=80,start=start)
    print('final N160 shallow-only fit',flush=True)
    th,pred,r,info=fit(rows,n=160,start=coarse)
    save_predictions('figure6_multistage_predictions.csv',rows,pred)
    save_predictions('figure6_multistage_predictions_n160.csv',rows,pred)
    print('fitting linear baseline N160',flush=True)
    lin,lp,lr,li=fit(rows,n=160,kind='linear');save_predictions('figure6_linear_predictions.csv',rows,lp)
    print('fixed-parameter refinement',flush=True)
    refinements=[]
    for n in [80,160,320]:
        yp,rr=predict(rows,th,n=n);refinements.append(dict(n=n,training=metrics(rows[:9],yp[:9]),spatial_holdout=metrics(rows[9:],yp[9:]),max_change_from_n160=float(abs(yp-pred).max())))
    print('assumption sensitivities',flush=True)
    sensitivities=[]
    for cin,offset in [(1e-7,.1),(1e-5,.1),(1e-6,.08),(1e-6,.15)]:
        tt,pp,rr,ss=fit(rows,cin=cin,offset=offset,n=160,start=th);sensitivities.append(ss)
    rows7=read_points(7);pred7,r7=predict(rows7,th,n=160);save_predictions('figure7_unknown_condition_diagnostic.csv',rows7,pred7,role='diagnostic')
    diagnostic={str(d):metrics([v for v in rows7 if abs(float(v['depth_m'])-d)<1e-5],pred7[np.array([abs(float(v['depth_m'])-d)<1e-5 for v in rows7])]) for d in DEPTHS}
    out=dict(source='Gitis et al. Chemical Engineering Journal163(2010)78–85, DOI10.1016/j.cej.2010.07.044; experimental symbols digitized from Figure6',
        purpose='Conditional normalized-shape comparison, not complete physical calibration or validation of four stages.',
        fixed_protocol=dict(m0=.44,column_length_m=1.62,Darcy_rate_m_per_h=10.,mode='rate',dimensionless_forcing=1.,tstar_per_hour=U/L),
        assumptions=['Unknown physical inlet concentration replaced by assumed volume fraction1e-6; sensitivities1e-7/1e-5.','Deposit porosity eps=.5 is fixed and assumed, not measured.','Published time zero contains particles at all depths; assume observations start0.1h after clean-bed model inlet step; sensitivities.08/.15h.','Fit only9 shallow-depth Figure6 symbols. Hold out18 deeper Figure6 symbols.','beta_s=0 is an explicit fitted-model restriction. Aging and capacity thresholds fixed far above realized loading; no aging/capacity stage is validated.','Figure7 conditions not independently confirmed; no refit and no controlled-validation claim.','Digitization envelopes are image extraction bounds, not measurement standard deviations. Unweighted least squares used.'],
        main_multistage=info,linear_baseline=li,refinement_fixed_parameters=refinements,
        coarse_n80_initial_fit=coarseinfo,coarse_n80_after_global_search=coarseinfo2,
        global_search=dict(seed=20261125,popsize=8,maxiter=45,number_function_evaluations=int(glob.nfev),SSE=float(glob.fun),fit_data='Only9 Figure6 symbols at0.82m; no holdout or Figure7 symbols'),
        assumption_sensitivity=sensitivities,figure7_external_diagnostic_by_depth=diagnostic,
        fitted_theta_log=th.tolist(),fitted_linear_theta=lin.tolist())
    (OUT/'fit_gitis_results.json').write_text(json.dumps(out,indent=2));plots(rows,th,lin,out)
    with (OUT/'fit_metrics_summary.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['model','data_role','n','RMSE','MAE','bias','R2'])
        for name,source in [('multistage',info),('linear',li)]:
            for role in ['training','spatial_holdout']:
                v=source[role];w.writerow([name,role,v['n'],v['RMSE'],v['MAE'],v['bias'],v['R2']])
            for depth,v in source['by_depth'].items():w.writerow([name,'depth_'+depth,v['n'],v['RMSE'],v['MAE'],v['bias'],v['R2']])
    with (OUT/'fit_assumption_sensitivity.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['cin_assumed','clock_offset_h','alpha_s','charging','B1_over_cin','training_RMSE','spatial_holdout_RMSE'])
        for v in sensitivities:w.writerow([v['cin_assumed'],v['clock_offset_h_assumed'],v['alpha_s'],v['charging'],v['B1_over_cin'],v['training']['RMSE'],v['spatial_holdout']['RMSE']])
    print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
