"""Publication plots from fixed Alem observations and N=240 saved predictions.

This script only reads observations/results and renders figures. It does not fit
models, alter data or sample/interpolate any published source-image curve.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

BASE = Path(__file__).resolve().parent / 'alem'
RESULTS = BASE / 'results'
STYLE = {'linear': {'color':'#2267A0','linestyle':(0,(5,2.3)),'label':'Linear kinetics'},
         'multistage': {'color':'#D16A24','linestyle':'-','label':'Multistage kinetics'}}
plt.rcParams.update({
    'font.family':'DejaVu Sans', 'font.size':12,
    'axes.labelsize':12, 'axes.titlesize':12, 'xtick.labelsize':11,
    'ytick.labelsize':11, 'legend.fontsize':11, 'axes.linewidth':0.85,
    'lines.linewidth':1.9, 'pdf.fonttype':42, 'ps.fonttype':42,
    'svg.fonttype':'none', 'savefig.facecolor':'white',
    'axes.spines.top':False, 'axes.spines.right':False,
    'xtick.direction':'out', 'ytick.direction':'out',
    'xtick.major.size':4, 'ytick.major.size':4,
})


def read_csv(path: Path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def filehash(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


observations = {
    'C_over_C0':read_csv(BASE / 'data_breakthrough.csv'),
    'K_over_K0':read_csv(BASE / 'data_hydraulics.csv'),
}
fields = {(model,mode):np.load(RESULTS/f'{model}_{mode}_fields.npz')
          for model in STYLE for mode in ('rate','pressure')}
for (model,mode), values in fields.items():
    if values['B_profile'].shape[1] != 240:
        raise ValueError(f'{model} {mode}: expected N=240 saved fields')
    if not np.all(np.diff(values['pv']) > 0):
        raise ValueError('Pore-volume coordinates must increase strictly')


def panel_style(ax, letter, title):
    ax.set_title(title, loc='left', pad=13)
    ax.text(-0.115, 1.075, f'({letter})', transform=ax.transAxes,
            ha='left', va='bottom', fontsize=12, fontweight='bold')
    ax.grid(axis='y', color='#D5DADF', linewidth=0.65, alpha=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(pad=5)


def plot_observations(ax, quantity, mode):
    selected=[r for r in observations[quantity] if r['mode']==mode]
    for ambiguous in (False,True):
        rows=[r for r in selected if (r['overlapping_series'].lower()=='true')==ambiguous]
        if not rows:
            continue
        ax.errorbar([float(r['pore_volumes']) for r in rows],
                    [float(r['value']) for r in rows],
                    xerr=[float(r['digitization_x_bound_PV']) for r in rows],
                    yerr=[float(r['digitization_value_bound']) for r in rows],
                    fmt='s' if ambiguous else 'o',
                    markersize=4.8 if ambiguous else 4.1,
                    markerfacecolor='white' if ambiguous else '#20252A',
                    markeredgecolor='#20252A', markeredgewidth=1.0,
                    ecolor='#737C85', elinewidth=0.8, capsize=2.1,
                    capthick=0.8, linestyle='none', zorder=5)


def save(fig, stem):
    for ext in ('pdf','svg','png'):
        fig.savefig(RESULTS / f'{stem}.{ext}', dpi=600 if ext=='png' else 300,
                    bbox_inches='tight', pad_inches=0.10)


fig, axs=plt.subplots(2,2,figsize=(11.4,8.35),sharex=True)
fig.subplots_adjust(left=0.105,right=0.98,bottom=0.215,top=0.91,
                    hspace=0.41,wspace=0.28)
# Full saved curves are shown, including early excursions.
max_c=max(float(np.max(v['c_normalized'])) for v in fields.values())
c_ylim_top=max(0.6, np.ceil((max_c+0.025)/0.2)*0.2)
for column,(field_mode,obs_mode) in enumerate((('rate','constant_flow'),('pressure','constant_head'))):
    for model,style in STYLE.items():
        d=fields[(model,field_mode)]
        axs[0,column].plot(d['pv'], d['c_normalized'][:,0],
                           color=style['color'], linestyle=style['linestyle'],
                           zorder=3,label=style['label'])
        # Use effective intrinsic permeability, not flow/pressure conductance.
        axs[1,column].plot(d['pv'], d['keff'],
                           color=style['color'], linestyle=style['linestyle'],zorder=3)
    plot_observations(axs[0,column],'C_over_C0',obs_mode)
    plot_observations(axs[1,column],'K_over_K0',obs_mode)
    axs[0,column].set_ylim(0,c_ylim_top)
    axs[0,column].set_yticks(np.arange(0,c_ylim_top+1e-8,0.2))
    axs[1,column].set_ylim(0.66,1.025)
    axs[1,column].set_yticks([0.7,0.8,0.9,1.0])
    axs[1,column].set_xlabel('Injected pore volumes, PV')
    for row in range(2):
        axs[row,column].set_xlim(0,85)
        axs[row,column].set_xticks([0,20,40,60,80])
        axs[row,column].set_ylabel(r'$C/C_0$' if row==0 else r'$K/K_0$')
panel_style(axs[0,0],'a','Constant flow: breakthrough calibration')
panel_style(axs[0,1],'b','Constant head: held-out breakthrough')
panel_style(axs[1,0],'c','Constant flow: held-out hydraulics')
panel_style(axs[1,1],'d','Constant head: held-out hydraulics')
legend_handles=[
    Line2D([],[],marker='o',color='none',markerfacecolor='#20252A',
           markeredgecolor='#20252A',markersize=5,label='Digitized observations'),
    Line2D([],[],marker='s',color='none',markerfacecolor='white',
           markeredgecolor='#20252A',markersize=5,label='Overlapping source symbols'),
    *[Line2D([],[],color=s['color'],linestyle=s['linestyle'],linewidth=1.9,label=s['label'])
      for s in STYLE.values()]
]
fig.legend(handles=legend_handles,loc='lower center',bbox_to_anchor=(0.53,0.09),
           ncol=2,frameon=False,columnspacing=2.7,handlelength=3.0)
# Measurement definitions and digitization bounds are stated in the article
# caption and Methods, avoiding small duplicate text below the panel legend.
save(fig,'alem_comparison')
plt.close(fig)

# The nine retained-mass bins are preserved exactly; model cell values are drawn
# against the physical column coordinate. Horizontal bin widths are intervals,
# not measurement error bars.
profile=read_csv(BASE / 'data_retention_profile_table1_flow.csv')
left=np.array([float(r['x_left_m']) for r in profile])
right=np.array([float(r['x_right_m']) for r in profile])
B=np.array([float(r['B_mean_in_bin_per_original_pore_volume']) for r in profile])
edges=np.r_[left,right[-1]]
if not np.allclose(right[:-1],left[1:]):
    raise ValueError('Measured retention bins are not contiguous')
fig,ax=plt.subplots(figsize=(8.9,5.35))
fig.subplots_adjust(left=0.13,right=0.98,bottom=0.27,top=0.88)
ax.bar((left+right)/2, B*1e3, width=right-left,
       align='center',facecolor='#E5E8EB',edgecolor='#737C85',
       linewidth=0.8,zorder=1,label='Measured bin means (Table 1)')
ax.stairs(B*1e3, edges,baseline=None,color='#323B43',linewidth=1.3,zorder=3)
for model,style in STYLE.items():
    d=fields[(model,'rate')]
    # x is the dimensionless cell-center position; physical length is 0.4 m.
    ax.plot(d['x']*0.4,d['B_profile'][-1]*1e3,
            color=style['color'],linestyle=style['linestyle'],
            linewidth=2.1,zorder=4,label=style['label'])
ax.set_xlabel('Depth into the column, m')
ax.set_ylabel(r'Retained solid volume, $10^3 B$')
ax.set_xlim(0,0.4)
ax.set_ylim(0,max(float(np.max(B*1e3)),float(np.max(fields[('linear','rate')]['B_profile'][-1]*1e3)),float(np.max(fields[('multistage','rate')]['B_profile'][-1]*1e3)))*1.10)
ax.set_xticks([0,.1,.2,.3,.4])
ax.set_yticks([0,4,8,12,16])
ax.set_title('Retained-solids profile after 83 PV: constant flow',loc='left',pad=14)
ax.grid(axis='y',color='#D5DADF',linewidth=0.65,alpha=0.7)
ax.set_axisbelow(True)
handles=[Patch(facecolor='#E5E8EB',edgecolor='#737C85',label='Measured bin means (Table 1)'),
         *[Line2D([],[],color=s['color'],linestyle=s['linestyle'],linewidth=2.1,label=s['label']) for s in STYLE.values()]]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.55,.102),
           ncol=3,frameon=False,columnspacing=1.7,handlelength=2.8,fontsize=10.5)
fig.text(.13,.045,
         r'$B$: solid volume per initial pore volume. Bars span the nine sampled bins; no uncertainty bars supplied.',
         fontsize=10,color='#4B535B')
fig.text(.13,.014,
         'Lines are saved model profiles at N=240; hydraulic consistency bounds are reported separately.',
         fontsize=10,color='#4B535B')
save(fig,'alem_retention_profile')
plt.close(fig)

inputs=[BASE/'data_breakthrough.csv',BASE/'data_hydraulics.csv',
        BASE/'data_retention_profile_table1_flow.csv',
        *[RESULTS/f'{m}_{k}_fields.npz' for m in STYLE for k in ('rate','pressure')]]
(RESULTS/'plot_provenance.json').write_text(json.dumps({
    'script':'work_validation/plot_alem.py',
    'source_doi':'10.1007/s11242-014-0402-8',
    'saved_model_grid_cells':240,'effective_permeability_field':'keff',
    'hydraulic_conductance_used':False,'data_modified':False,
    'bounds':'Digitization envelopes only, not experimental uncertainty or confidence intervals.',
    'overlap_flag':'Open squares mark observations with overlapping source symbols; each CSV-provided bound is preserved.',
    'profile':'Printed retained-mass fractions converted to B in an existing input table; visual axis shows 1000*B.',
    'inputs_sha256':{str(p.relative_to(BASE)):filehash(p) for p in inputs},
    'outputs':['alem_comparison.pdf','alem_comparison.svg','alem_comparison.png',
               'alem_retention_profile.pdf','alem_retention_profile.svg','alem_retention_profile.png'],
},indent=2)+'\n')
print('Rendered Alem assessment and retention-profile figures (PDF, SVG, PNG).')
