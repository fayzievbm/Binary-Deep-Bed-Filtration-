"""2017 two-component filtration equations with ONLY piecewise deposition kinetics changed.

NumPy/Matplotlib are the only dependencies. Parameters are illustrative,
dimensionless and uncalibrated. Run python solver.py to regenerate this study.
"""
from dataclasses import dataclass, replace, asdict
from pathlib import Path
import numpy as np
import json,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent

@dataclass
class Params:
    m0:float=.35
    eps:tuple=(.5,.3)
    alpha_s:tuple=(1.,1.)
    alpha_d:tuple=(0.,0.)
    beta_s:tuple=(1.,1.)
    beta_d:tuple=(0.,0.)
    charging:tuple=(.5,.5)
    B1:tuple=(.001,.0002)
    B2:tuple=(.004,.001)
    cap:tuple=(.02,.004)
    mode:str='pressure'
    forcing:float=10.
    kinetics:str='multistage'

def fields(y,p):
    B=y[1];D=np.sum(B/(1-np.asarray(p.eps)[:,None]),axis=0);phi=1-D
    c=y[0]/phi; k=(1-np.sqrt(np.maximum(D,0)))**3
    mu=1+2.5*c.sum(axis=0);res=np.mean(mu/k)
    u=p.forcing/res if p.mode=='pressure' else p.forcing
    dp=u*res;grad=u*mu/k
    alpha=np.asarray(p.alpha_s)[:,None]+np.asarray(p.alpha_d)[:,None]*u
    beta=np.asarray(p.beta_s)[:,None]+np.asarray(p.beta_d)[:,None]*grad
    return D,phi,c,k,u,dp,alpha,beta

def safe_log_time(start,target,equilibrium,rate):
    """Positive event time for y'=rate*(equilibrium-y), or +inf."""
    ratio=np.divide(target-equilibrium,start-equilibrium,out=np.ones_like(start),where=abs(start-equilibrium)>1e-300)
    valid=(ratio>0)&(ratio<1)&(rate>0)
    out=np.full_like(start,np.inf)
    out[valid]=-np.log(ratio[valid])/rate[valid]
    return out

def source_advance(B,c,alpha,beta,dt,p):
    """Exact branch integration with c,alpha,beta frozen; conservative donor transfer.

    The attractive B1 discontinuity is completed by holding at B1 when
    q_left>=0>=q_right. The terminal cap arrests the entire net exchange.
    """
    g=alpha*c; rate=beta/p.m0
    eq=np.divide(g,beta,out=np.zeros_like(g),where=beta>0)
    if p.kinetics=='linear':
        out=np.where(beta>0,eq+(B-eq)*np.exp(-rate*dt),B+g*dt/p.m0)
        return out,0,0
    b=B.copy();rem=np.full_like(B,dt);hold_count=0;cap_hits=0
    lo=np.broadcast_to(np.asarray(p.B1)[:,None],B.shape)
    mid=np.broadcast_to(np.asarray(p.B2)[:,None],B.shape)
    cap=np.broadcast_to(np.asarray(p.cap)[:,None],B.shape)
    chi=np.broadcast_to(np.asarray(p.charging)[:,None],B.shape)
    rate=np.broadcast_to(rate,B.shape);eq=np.broadcast_to(eq,B.shape)
    for it in range(8):
        active=rem>1e-15
        if not active.any():break
        terminal=active&(b>=cap)
        rem[terminal]=0
        # State exactly at B1: Filippov/one-sided holding of irreversible charged layer.
        hold=active&(b==lo)&(g-beta*lo<=0)
        hold_count+=int(hold.sum());rem[hold]=0
        active=rem>1e-15
        charging=active&(b<lo)
        # At an exactly reached B2 interface, use the lower (printed equality)
        # branch when its derivative is nonpositive; upward crossings continue
        # immediately into the aging branch. This also resolves repelling B2.
        mid_down=(b==mid)&(g-beta*mid<=0)
        transition=active&(b>=lo)&((b<mid)|mid_down)
        aging=active&((b>mid)|((b==mid)&~mid_down))&(b<cap)
        # Charging is linear in time at frozen concentration.
        v=chi*g/p.m0
        event=np.divide(lo-b,v,out=np.full_like(b,np.inf),where=v>0)
        step=np.minimum(rem,event)
        crossed=charging&(event<=rem)
        b[charging]+=v[charging]*step[charging];b[crossed]=lo[crossed]
        rem[charging]-=step[charging]
        # Transition branch: exponential approach to g/beta.
        q=g-beta*b
        target=np.where(q>=0,mid,lo)
        event=safe_log_time(b,target,eq,rate)
        zero_rate=(rate==0)&(g>0)
        event[zero_rate]=(target[zero_rate]-b[zero_rate])*p.m0/g[zero_rate]
        step=np.minimum(rem,event);crossed=transition&(event<=rem)
        out=np.where(beta>0,eq+(b-eq)*np.exp(-rate*step),b+g*step/p.m0)
        b[transition]=out[transition];b[crossed]=target[crossed]
        rem[transition]-=step[transition]
        # Aging branch: B^2 obeys a linear ODE at frozen concentration.
        eq2=eq*cap;target=np.where(g*cap-beta*b*b>=0,cap,mid)
        event=safe_log_time(b*b,target*target,eq2,2*rate)
        event[zero_rate]=(target[zero_rate]**2-b[zero_rate]**2)*p.m0/(2*g[zero_rate]*cap[zero_rate])
        step=np.minimum(rem,event);crossed=aging&(event<=rem)
        out2=np.where(beta>0,eq2+(b*b-eq2)*np.exp(-2*rate*step),b*b+2*g*cap*step/p.m0)
        out=np.sqrt(np.maximum(out2,0))
        b[aging]=out[aging];b[crossed]=target[crossed]
        cap_hits+=int((crossed&(target==cap)).sum())
        rem[aging]-=step[aging]
    if np.max(rem)>1e-12:raise RuntimeError('Unresolved source event cycle')
    return b,hold_count,cap_hits

def feed(t):return np.array([.05,.01])
def smooth_feed(t):return np.array([.05,.01])*.5*(1-np.cos(np.pi*min(t/.1,1)))

def simulate(p=Params(),n=160,tend=.8,inlet=feed,cfl=.8,dtmax=np.inf,record_dt=.002):
    dx=1/n;y=np.zeros((2,2,n));incoming=np.zeros(2);outgoing=np.zeros(2)
    maxerr=np.zeros(2);minstate=0.;min_carrier=1.;maxD=0.;maxc=0.;holds=0;caphits=0;steps=0;t=0.
    rec=sorted(set([0.,tend]+[float(v) for v in np.arange(0,tend+record_dt/2,record_dt) if v<=tend]+[v for v in [.01,.03,.08,.15,.25,.3,.8] if v<=tend]))
    ri=0;hist=[];snaps={}
    def record():
        D,phi,c,k,u,dp,alpha,beta=fields(y,p)
        hist.append([t,*c[:,-1],u,dp,*np.mean(y[1],axis=1),*np.mean(y.sum(axis=0),axis=1),*incoming,*outgoing,*maxerr])
        if any(abs(t-v)<1e-9 for v in [.01,.03,.08,.15,.25,.3,.8]):snaps[f'{t:.2f}']=y.copy()
    while t<tend-1e-13:
        if abs(t-rec[ri])<1e-10:record();ri+=1
        D,phi,c,k,u,dp,alpha,beta=fields(y,p)
        boost=np.maximum(1,np.maximum(np.asarray(p.charging),np.asarray(p.cap)/np.asarray(p.B2)))[:,None] if p.kinetics!='linear' else np.ones((2,1))
        dt=min(cfl*p.m0*phi.min()*dx/max(u,1e-300),cfl*np.min(p.m0*phi/(alpha*boost+1e-300)),dtmax,rec[ri]-t,tend-t)
        # Conservative upwind transport, shared hydraulic velocity at all faces.
        cin=inlet(t)
        faces=u*np.concatenate([cin[:,None],c],axis=1)
        yt=y.copy();yt[0]-=dt*(faces[:,1:]-faces[:,:-1])/(p.m0*dx)
        # Freeze post-transport concentration and apply exact branch/event source flow.
        _,_,cs,_,_,_,als,bes=fields(yt,p)
        newB,hc,ch=source_advance(yt[1],cs,als,bes,dt,p)
        transfer=newB-yt[1];yt[0]-=transfer;yt[1]=newB
        y=yt;incoming+=dt*faces[:,0]/p.m0;outgoing+=dt*faces[:,-1]/p.m0
        holds+=hc;caphits+=ch;t+=dt;steps+=1
        inv=np.mean(y.sum(axis=0),axis=1)
        maxerr=np.maximum(maxerr,abs(inv-incoming+outgoing))
        D,phi,c,*_=fields(y,p)
        S1=phi-y[0].sum(axis=0)
        min_carrier=min(min_carrier,float(S1.min()))
        minstate=min(minstate,float(y.min()),float(S1.min()))
        maxD=max(maxD,float(D.max()));maxc=max(maxc,float(c.sum(axis=0).max()))
        if minstate < -1e-12 or D.max()>=1:raise RuntimeError(('admissibility',t,minstate,D.max()))
    record();h=np.array(hist)
    H=np.asarray(p.eps)[:,None]/(1-np.asarray(p.eps)[:,None])*y[1]
    return dict(p=asdict(p),n=n,steps=steps,time=h[:,0],outlet=h[:,1:3],u=h[:,3],dp=h[:,4],retained=h[:,5:7],inventory=h[:,7:9],incoming=incoming,outgoing=outgoing,maxerr=maxerr,minstate=minstate,min_carrier=min_carrier,maxD=maxD,maxc=maxc,holds=holds,cap_hits=caphits,final=y,snapshots=snaps,hist=h,porosity_identity_error=float(abs((1-np.asarray(p.eps)[:,None])*H-np.asarray(p.eps)[:,None]*y[1]).max()))

def summary(r):
    dyn={}
    for t in [.03,.08,.15,.25]:
        j=int(np.argmin(abs(r['time']-t)));dyn[str(t)]=dict(flow=float(r['u'][j]),pressure=float(r['dp'][j]),retained=r['retained'][j].tolist(),outlet=r['outlet'][j].tolist())
    return dict(n=r['n'],steps=r['steps'],flow_final=float(r['u'][-1]),pressure_final=float(r['dp'][-1]),retained_final=r['retained'][-1].tolist(),input=r['incoming'].tolist(),output=r['outgoing'].tolist(),max_species_balance_error=r['maxerr'].tolist(),minimum_of_all_seven_saturations=r['minstate'],minimum_mobile_carrier_saturation=r['min_carrier'],maximum_occupancy=r['maxD'],minimum_mobile_fraction=1-r['maxD'],maximum_total_concentration=r['maxc'],B1_hold_events=r['holds'],cap_crossings=r['cap_hits'],fixed_porosity_identity_error=r['porosity_identity_error'],dynamic_samples=dyn,peak_normalized_outlet=(np.max(r['outlet'],axis=0)/np.array([.05,.01])).tolist())

def save(name,r):
    header='t_star,c1_out,c2_out,u_star,delta_p_star,B1_integral,B2_integral,inventory1,inventory2,input1,input2,output1,output2,max_balance1,max_balance2'
    np.savetxt(ROOT/(name+'.csv'),r['hist'],delimiter=',',header=header,comments='')
    np.savez_compressed(ROOT/(name+'_fields.npz'),final=r['final'],**{f't{t}':v for t,v in r['snapshots'].items()})

def event_benchmark():
    p=Params();C=np.array(p.alpha_s)*np.array([.05,.01]);be=np.array(p.beta_s);B1=np.array(p.B1);B2=np.array(p.B2);cap=np.array(p.cap);chi=np.array(p.charging)
    t1=p.m0*B1/(chi*C)
    t2=t1-p.m0/be*np.log((C/be-B2)/(C/be-B1))
    tc=t2-p.m0/(2*be)*np.log((C*cap/be-cap*cap)/(C*cap/be-B2*B2))
    errors=[]
    for t in [.005,.02,.08,.2]:
        expected=np.zeros(2)
        for i in range(2):
            if t<t1[i]:expected[i]=chi[i]*C[i]*t/p.m0
            elif t<t2[i]:expected[i]=C[i]/be[i]+(B1[i]-C[i]/be[i])*np.exp(-be[i]*(t-t1[i])/p.m0)
            elif t<tc[i]:expected[i]=np.sqrt(C[i]*cap[i]/be[i]+(B2[i]**2-C[i]*cap[i]/be[i])*np.exp(-2*be[i]*(t-t2[i])/p.m0))
            else:expected[i]=cap[i]
        actual,*_=source_advance(np.zeros((2,1)),np.array([[.05],[.01]]),np.ones((2,1)),np.ones((2,1)),t,p)
        errors.append(float(abs(actual[:,0]-expected).max()))
    # Original linear source, independent inlet formula from 2017.
    lp=replace(p,kinetics='linear');t=.03
    actual,*_=source_advance(np.zeros((2,1)),np.array([[.05],[.01]]),np.ones((2,1)),np.ones((2,1)),t,lp)
    exact=np.array([.05,.01])*(1-np.exp(-t/p.m0))
    # At the B1 discontinuity and insufficient concentration the charged layer holds.
    held,*_=source_advance(B1[:,None],np.zeros((2,1)),np.ones((2,1)),np.ones((2,1)),.1,p)
    arrested,*_=source_advance(cap[:,None],np.zeros((2,1)),np.ones((2,1)),np.ones((2,1)),.1,p)
    start=np.array([.006,.0015]);con=np.array([.0005,.00005]);totaldt=.3
    event=p.m0/2*np.log((start**2-con*cap)/(B2**2-con*cap))
    expected=con+(B2-con)*np.exp(-(totaldt-event)/p.m0)
    down,*_=source_advance(start[:,None],con[:,None],np.ones((2,1)),np.ones((2,1)),totaldt,p)
    return dict(charging_end_times=t1.tolist(),transition_end_times=t2.tolist(),cap_times=tc.tolist(),piecewise_max_absolute_error=max(errors),linear_max_absolute_error=float(abs(actual[:,0]-exact).max()),B1_holding_error=float(abs(held[:,0]-B1).max()),cap_terminal_arrest_error=float(abs(arrested[:,0]-cap).max()),downward_B2_event_error=float(abs(down[:,0]-expected).max()))

def convergence():
    p=Params();T=.2
    runs={n:simulate(p,n=n,tend=T,inlet=smooth_feed,record_dt=.01) for n in [80,160,320,640]}
    ref=runs[640]['final'];space=[]
    for n in [80,160,320]:
        projected=ref.reshape(2,2,n,640//n).mean(axis=-1)
        space.append(dict(n=n,L1=float(np.mean(abs(runs[n]['final']-projected)))))
    n=80;ref=simulate(p,n=n,tend=T,inlet=smooth_feed,dtmax=.0000125,record_dt=.01)['final'];temporal=[]
    for dt in [.0002,.0001,.00005]:
        y=simulate(p,n=n,tend=T,inlet=smooth_feed,dtmax=dt,record_dt=.01)['final']
        temporal.append(dict(dt=dt,L1=float(np.mean(abs(y-ref)))))
    return space,temporal

def overshoot_refinement():
    result=[];p=replace(Params(),alpha_s=(1.,3.))
    for n in [80,160,320,640]:
        r=simulate(p,n=n,tend=.25,record_dt=.0005)
        j=int(np.argmax(r['outlet'][:,1]))
        result.append(dict(n=n,peak_c2_normalized=float(r['outlet'][j,1]/.01),peak_time=float(r['time'][j]),maximum_total_concentration=r['maxc']))
    return result

def figures(runs,s):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300})
    styles={'linear':('#777777','--','Original linear kinetics'),'multistage':('#0072B2','-','Multistage kinetics')}
    fig,axs=plt.subplots(2,2,figsize=(9,6),layout='constrained')
    for name,(co,ls,label) in styles.items():
        r=runs[name];axs[0,0].plot(r['time'],r['outlet'][:,0]/.05,color=co,ls=ls,label=label);axs[0,1].plot(r['time'],r['outlet'][:,1]/.01,color=co,ls=ls)
        axs[1,0].plot(r['time'],r['u'],color=co,ls=ls);axs[1,1].plot(r['time'],r['retained'].sum(axis=1),color=co,ls=ls)
    axs[0,0].set_ylabel(r'$c_{1,out}/c_{1,in}$');axs[0,1].set_ylabel(r'$c_{2,out}/c_{2,in}$');axs[1,0].set_ylabel(r'Flow rate $u^*$');axs[1,1].set_ylabel('Retained solid inventory');axs[0,0].legend(fontsize=10)
    for ax,l in zip(axs.flat,'abcd'):ax.set_xlabel(r'Time $t^*=U_0t/L$');ax.grid(alpha=.15);ax.text(.02,.96,f'({l})',transform=ax.transAxes,va='top')
    fig.savefig(ROOT/'figure1_kinetics_comparison.png');fig.savefig(ROOT/'figure1_kinetics_comparison.pdf');plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(9,6),layout='constrained');r=runs['multistage'];x=(np.arange(r['n'])+.5)/r['n'];colors=['#0072B2','#D55E00','#009E73','#CC79A7']
    for t,col in zip(['0.01','0.03','0.08','0.15'],colors):
        y=r['snapshots'][t];p=Params();_,_,c,*_=fields(y,p)
        axs[0,0].plot(x,c[0]/.05,color=col,label=f't*={t}');axs[0,1].plot(x,c[1]/.01,color=col);axs[1,0].plot(x,y[1,0],color=col);axs[1,1].plot(x,y[1,1],color=col)
    axs[0,0].legend(fontsize=10);axs[0,0].set_ylabel(r'$c_1/c_{1,in}$');axs[0,1].set_ylabel(r'$c_2/c_{2,in}$');axs[1,0].set_ylabel(r'$B_1=S_{41}$');axs[1,1].set_ylabel(r'$B_2=S_{42}$')
    for ax,l in zip(axs.flat,'abcd'):ax.set_xlabel(r'$x^*=x/L$');ax.grid(alpha=.15);ax.text(.02,.96,f'({l})',transform=ax.transAxes,va='top')
    fig.savefig(ROOT/'figure2_multistage_profiles.png');fig.savefig(ROOT/'figure2_multistage_profiles.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(9,3.6),layout='constrained')
    for name,col,label in [('charging_slow','#D55E00',r'$\chi_1=0.25$'),('multistage','#0072B2',r'$\chi_1=0.5$'),('charging_fast','#009E73',r'$\chi_1=1$')]:
        r=runs[name];axs[0].plot(r['time'],r['outlet'][:,0]/.05,color=col,label=label);axs[1].plot(r['time'],r['retained'][:,0],color=col)
    for name,col,ls,label in [('multistage','#0072B2','-',r'$\alpha_{s2}^*=1$'),('species2_fast','#CC79A7','--',r'$\alpha_{s2}^*=3$')]:
        r=runs[name];axs[2].plot(r['time'],r['outlet'][:,1]/.01,color=col,ls=ls,label=label)
    axs[0].set_ylabel(r'$c_{1,out}/c_{1,in}$');axs[1].set_ylabel('Type 1 retained inventory');axs[2].set_ylabel(r'$c_{2,out}/c_{2,in}$');axs[0].legend(fontsize=10);axs[2].legend(fontsize=10)
    for ax,l in zip(axs,'abc'):ax.set_xlabel(r'Time $t^*=U_0t/L$');ax.grid(alpha=.15);ax.text(.02,.96,f'({l})',transform=ax.transAxes,va='top')
    fig.savefig(ROOT/'figure3_kinetic_sensitivity.png');fig.savefig(ROOT/'figure3_kinetic_sensitivity.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(9,3.6),layout='constrained')
    sp=s['spatial_convergence'];tc=s['temporal_convergence'];axs[0].loglog([1/v['n'] for v in sp],[v['L1'] for v in sp],'o-',color='#0072B2');axs[0].set_xticks([1/320,1/160,1/80],labels=['1/320','1/160','1/80']);axs[0].set_xlabel('Grid spacing');axs[0].set_ylabel('Mean absolute field error')
    axs[1].loglog([v['dt'] for v in tc],[v['L1'] for v in tc],'o-',color='#D55E00');axs[1].set_xticks([.00005,.0001,.0002],labels=['0.5','1','2']);axs[1].set_xlabel(r'Time step $\times10^{-4}$');axs[1].set_ylabel('Mean absolute field error')
    p=Params();times=np.linspace(0,.2,401);B=np.array([source_advance(np.zeros((2,1)),np.array([[.05],[.01]]),np.ones((2,1)),np.ones((2,1)),t,p)[0][:,0] for t in times]);axs[2].plot(times,B[:,0],label='Type 1');axs[2].plot(times,B[:,1],label='Type 2');axs[2].set_xlabel(r'Time $t^*$');axs[2].set_ylabel('Constant-concentration deposit');axs[2].legend(fontsize=10)
    from matplotlib.ticker import NullLocator
    for ax,l in zip(axs,'abc'):ax.xaxis.set_minor_locator(NullLocator());ax.grid(alpha=.15);ax.margins(x=.1);ax.tick_params(labelsize=11);ax.text(.02,.96,f'({l})',transform=ax.transAxes,va='top')
    fig.savefig(ROOT/'figure4_verification.png');fig.savefig(ROOT/'figure4_verification.pdf');plt.close(fig)

def main():
    p=Params();cases={'linear':replace(p,kinetics='linear'),'multistage':p,'charging_slow':replace(p,charging=(.25,.5)),'charging_fast':replace(p,charging=(1.,.5)),'species2_fast':replace(p,alpha_s=(1.,3.)),'linear_fixed_rate':replace(p,kinetics='linear',mode='rate'),'multistage_fixed_rate':replace(p,mode='rate')}
    runs={}
    for name,q in cases.items():
        print('running',name,flush=True);runs[name]=simulate(q,n=320);save(name,runs[name])
    print('verification',flush=True);space,temporal=convergence()
    s={k:summary(r) for k,r in runs.items()};s['parameters']=asdict(p);s['event_benchmark']=event_benchmark();s['spatial_convergence']=space;s['temporal_convergence']=temporal
    s['overshoot_refinement']=overshoot_refinement()
    Dcap=sum(np.asarray(p.cap)/(1-np.asarray(p.eps)));kcap=(1-np.sqrt(Dcap))**3;mu=1.15
    exactu=10*kcap/mu;exactdp=10*mu/kcap
    s['coupled_steady_benchmark']=dict(Dcap=float(Dcap),relative_permeability=float(kcap),exact_flow=exactu,numerical_flow=float(runs['multistage']['u'][-1]),flow_absolute_error=abs(float(runs['multistage']['u'][-1])-exactu),exact_pressure=exactdp,numerical_pressure=float(runs['multistage_fixed_rate']['dp'][-1]),pressure_absolute_error=abs(float(runs['multistage_fixed_rate']['dp'][-1])-exactdp))
    (ROOT/'summary.json').write_text(json.dumps(s,indent=2))
    with (ROOT/'case_metrics.csv').open('w') as f:
        w=csv.writer(f);w.writerow(['case','flow_final','pressure_final','retained_total_final','max_balance_error','maxD','max_total_c','B1_hold_events','cap_crossings'])
        for name in runs:
            v=s[name];w.writerow([name,v['flow_final'],v['pressure_final'],sum(v['retained_final']),max(v['max_species_balance_error']),v['maximum_occupancy'],v['maximum_total_concentration'],v['B1_hold_events'],v['cap_crossings']])
    figures(runs,s);print(json.dumps(s,indent=2),flush=True)

if __name__=='__main__':main()
