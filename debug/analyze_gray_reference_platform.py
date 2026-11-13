"""Summarize the native Gray reference audit without conflating residual norms."""
from pathlib import Path
import csv
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'results/gray_reference_platform'


def read(mode):
    folder=FOLDER/mode
    with (folder/'equation_residual.csv').open() as f:
        rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
    log=(folder/'run.log').read_text()
    stop=re.search(r'NATIVE_FIRST_STOP\s+(\d+)\s+([\d.Ee+-]+)',log)
    macro=[]
    for line in log.splitlines():
        match=re.fullmatch(r'\s*(\d+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*',line)
        if match:macro.append((int(match[1]),float(match[2])))
    return rows,macro,stop,('**** RUN END' in log)


def main():
    gsis,gm,gstop,gdone=read('gsis');cis,cm,cstop,cdone=read('cis')
    assert gdone and cdone,'Reference processes must finish cleanly before final analysis.'
    assert abs(gsis[0]['common_initial_absolute']-cis[0]['common_initial_absolute'])<1e-14
    assert max(r['relative_source_gap'] for r in cis)<1e-12
    assert max(abs(r['absolute_original']-r['absolute_coupled']) for r in cis)<1e-12
    tail=[r for r in gsis if r['step']>=40]
    original=np.array([r['initial_relative_original'] for r in tail])
    span=float(np.ptp(original)/np.mean(original))
    last=gsis[-1]
    summary=dict(model='Gray phonon Callaway thermal cavity (not gas lid-driven cavity)',
        triangles=200,polynomial_degree=3,angular_points=[80,40],tau_R=.01,tau_N=.01,
        wall_temperatures=[0,1,0,0],stabilization=[1,1,1],
        tolerance=1e-5,gsis_native_first_stop=int(gstop[1]),
        gsis_update_at_native_stop=float(gstop[2]),gsis_extended_final_step=int(last['step']),
        gsis_final_macro_update=last['macro_update'],
        gsis_original_absolute=last['absolute_original'],gsis_coupled_absolute=last['absolute_coupled'],
        gsis_original_relative_initial=last['initial_relative_original'],
        gsis_coupled_relative_initial=last['initial_relative_coupled'],
        gsis_original_relative_rhs=last['rhs_relative_original'],
        gsis_coupled_relative_rhs=last['rhs_relative_coupled'],
        gsis_relative_macro_kinetic_gap=last['relative_source_gap'],
        plateau_relative_span_steps40to100=span,
        cis_final_step=int(cis[-1]['step']),cis_native_converged=cstop is not None,
        cis_final_original_relative_initial=cis[-1]['initial_relative_original'],
        cis_final_source_gap=cis[-1]['relative_source_gap'],
        unchanged_solver_equations=json.loads((FOLDER/'manifest.json').read_text())['unchanged_transport'],
        unchanged_macro_equations=json.loads((FOLDER/'manifest.json').read_text())['unchanged_synthetic_acceleration'])
    assert span<1e-6 and last['initial_relative_coupled']<1e-10
    (FOLDER/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    selected=[]
    for r in gsis:
        if int(r['step']) in [8,10,20,40,60,80,100]:selected.append(dict(mode='GSIS',**r))
    for r in cis:
        if int(r['step']) in [100,200,500]:selected.append(dict(mode='CIS',**r))
    with (FOLDER/'summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(selected[0]));writer.writeheader();writer.writerows(selected)
    fig,ax=plt.subplots(1,3,figsize=(14,4.2),constrained_layout=True)
    steps,updates=np.array(gm).T
    ax[0].semilogy(steps,updates,label='Macro update')
    ax[0].axhline(1e-5,color='gray',ls=':',label='Native tolerance')
    ax[0].axvline(int(gstop[1]),color='gray',lw=.8)
    ax[0].set_title('Native GSIS stopping metric')
    gx=[r['step'] for r in gsis if r['step']>0]
    data=[r for r in gsis if r['step']>0]
    ax[1].semilogy(gx,[r['initial_relative_original'] for r in data],'o-',label='Original kinetic equation')
    ax[1].semilogy(gx,[r['initial_relative_coupled'] for r in data],'s-',label='Returned-source equation')
    ax[1].set_title('Both divided by initial original residual')
    ax[2].semilogy(gx,[r['relative_source_gap'] for r in data],'o-',label='Macro vs VDF moments')
    ax[2].set_title('Relative L2 macro/kinetic moment gap')
    for a in ax:
        a.set_xlabel('Outer iteration');a.grid(True,which='both',alpha=.25);a.legend(fontsize=8)
    fig.suptitle('Uploaded Gray reference: original P3 thermal cavity, 80 x 40 angles')
    fig.savefig(FOLDER/'reference_platform.png',dpi=180);plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
