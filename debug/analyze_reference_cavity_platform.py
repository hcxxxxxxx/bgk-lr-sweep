"""Validate the original-cavity continuation and plot the measured residuals."""
from pathlib import Path
import csv
import hashlib
import json
import os
import re

os.environ.setdefault('MPLCONFIGDIR','/private/tmp/bgk_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'results/reference_cavity_residual_platform'
SAVED=ROOT/'reference_runs/nonoscillatory/original_config'
ORIGINAL=Path('/Users/hanchenxi/Desktop/NonOscillatory')


def load(folder):
    with (folder/'equation_residual.csv').open() as stream:
        return [{key:float(value) for key,value in row.items()} for row in csv.DictReader(stream)]


def main():
    settings=json.loads((FOLDER/'settings.json').read_text())
    manifest=json.loads((FOLDER/'solver_module_sha256.json').read_text())
    for name,digest in manifest.items():
        assert hashlib.sha256((FOLDER/'source'/name).read_bytes()).hexdigest()==digest
    initial=load(FOLDER/'zero_state')
    assert len(initial)==1 and initial[0]['step']==0
    scale=settings['common_zero_state_equation_residual']
    report=dict(settings=settings,solver_modules_unchanged=True,cases=[])
    def without_stopping(text):
        for key in ['TOL','TMAX']:
            text=re.sub(rf'(?m)^(\s*{key}\s*=)[^\n]*',key+'=AUDIT_VALUE',text)
        return text
    original_control=without_stopping((ORIGINAL/'control.in').read_text())
    mesh_digest=hashlib.sha256((ORIGINAL/'Rect_A1_Nx17_Ny17.msh').read_bytes()).hexdigest()
    for mode in ['zero_state','cis','gsis']:
        assert without_stopping((FOLDER/mode/'control.in').read_text())==original_control
        assert hashlib.sha256((FOLDER/mode/'Rect_A1_Nx17_Ny17.msh').read_bytes()).hexdigest()==mesh_digest
    report['original_control_unchanged_except_stopping']=True
    report['original_mesh_sha256']=mesh_digest
    collected={};rows=[]
    for mode,offset in [('cis',16),('gsis',50)]:
        assert (FOLDER/mode/'final_status.dat').exists()
        data=load(FOLDER/mode)
        oldstatus=(SAVED/mode/'final_status.dat').read_text().splitlines()[1].split()
        assert abs(data[0]['mixed_mass']-float(oldstatus[-1]))<1e-12
        for row in data:
            assert abs(row['relative_original']-row['absolute_original']/scale)<1e-15
            assert abs(row['absolute_original']-row['absolute_source'])<=row['absolute_coupled']+1e-12
            assert abs(row['absolute_original']-row['absolute_coupled'])<=row['absolute_source']+1e-12
            row['total_step']=offset+int(row['step'])
        # The stage-zero update was initialized for monitoring; replace it
        # with the saved update solely for plotting, not as a new iteration.
        data[0]['macro_update']=float(oldstatus[2])
        if mode=='cis':
            assert max(abs(row['absolute_original']-row['absolute_coupled']) for row in data)<1e-12
        final=data[-1]
        tail=data[-min(4,len(data)):]
        report['cases'].append(dict(mode=mode.upper(),restart_step=offset,
                                   total_step=final['total_step'],converged=final['macro_update']<1e-11,
                                   final=final,
                                   tail_original_relative_range=(max(r['relative_original'] for r in tail)-min(r['relative_original'] for r in tail))/final['relative_original'],
                                   kinetic_mass_change_during_continuation=final['kinetic_mass']-data[0]['kinetic_mass']))
        collected[mode]=data
        rows.extend(dict(mode=mode.upper(),**row) for row in data)
    with (FOLDER/'validated_summary.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    c,g=[np.loadtxt(FOLDER/mode/'final_cells.dat') for mode in ['cis','gsis']]
    assert np.allclose(c[:,:4],g[:,:4])
    area=c[:,1]
    norm=lambda value:np.sqrt(np.sum(area[:,None]*value**2))
    report['gsis_vs_cis_kinetic_velocity_percent']=100*norm(g[:,5:7]-c[:,5:7])/norm(c[:,5:7])
    previous=np.loadtxt(SAVED/'gsis/final_cells.dat')
    density_change=g[:,4]-previous[:,4]
    mean_density_change=np.sum(area*density_change)/np.sum(area)
    report['gsis_step50_to_final']=dict(
        velocity_vector_change_percent=100*norm(g[:,5:7]-previous[:,5:7])/norm(previous[:,5:7]),
        mean_density_change=mean_density_change,
        maximum_departure_from_uniform_density_shift=float(np.max(abs(density_change-mean_density_change))))
    report['note']='The GSIS state must not be called a converged solution unless its macro-update stopping test is satisfied. Residual plateaus can coexist with secular density drift.'

    fig,axes=plt.subplots(1,3,figsize=(13.5,4.1))
    colors={'cis':'#2878b5','gsis':'#d95319'}
    for mode in ['cis','gsis']:
        data=collected[mode];x=[r['total_step'] for r in data]
        axes[0].semilogy(x,[r['macro_update'] for r in data],'o-',ms=3,color=colors[mode],label=mode.upper())
        axes[1].semilogy(x,[r['relative_original'] for r in data],'o-',ms=3,color=colors[mode],label=f'{mode.upper()}: original')
        axes[2].plot(x,[r['kinetic_mass'] for r in data],'o-',ms=3,color=colors[mode],label=f'{mode.upper()}: kinetic mass')
    data=collected['gsis'];x=[r['total_step'] for r in data]
    axes[1].semilogy(x,[r['relative_coupled'] for r in data],'s--',ms=3,color='#7e2f8e',label='GSIS: coupled')
    axes[1].semilogy(x,[r['relative_source'] for r in data],':',lw=2,color='#252525',label='GSIS: source gap')
    axes[0].axhline(1e-5,color='gray',ls=':',lw=1,label='original tolerance 1e-5')
    axes[0].axhline(1e-11,color='black',ls=':',lw=1,label='requested tolerance 1e-11')
    axes[0].set_title('Original macro-update stopping metric');axes[0].set_ylabel('Relative update')
    axes[1].set_title('Equation residuals');axes[1].set_ylabel('Residual / common zero-state residual')
    axes[2].set_title('Domain-integrated density');axes[2].set_ylabel('Kinetic mass')
    for ax in axes:
        ax.set_xlabel('Total iteration (including saved run)');ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.suptitle('Original lid-driven cavity: P4, 578 triangles, 24x24x24 velocities, delta = 0.099999653',fontsize=11)
    fig.text(.5,.015,'Continuation restores both H and mixed collision sources; original solver modules and stabilization are unchanged.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,.94));fig.savefig(FOLDER/'cavity_residual_platform.png',dpi=180);plt.close(fig)
    (FOLDER/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
