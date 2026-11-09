"""Reproducible PPT figures: solution differences, residuals, grids, literature.

Run debug/run_channel_presentation_checks.m in MATLAB first, then this script.
Numerical solver and reference stopping are unchanged by plotting.
"""
from pathlib import Path
import csv
import html
import json
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from matplotlib.transforms import Bbox
import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'results/channel_presentation_20261109'
OUT = ROOT / 'output/channel_figures_20261109'
OUT.mkdir(parents=True, exist_ok=True)
CSV = OUT / 'data'
CSV.mkdir(exist_ok=True)
PANELS=OUT/'pdf_panels'
PANELS.mkdir(exist_ok=True)
BLUE, RED, GOLD, BLACK = '#2568A6', '#BF4F39', '#AD811E', '#303C49'
FLOW = {'couette': 'Couette', 'poiseuille': 'Poiseuille'}
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 12, 'axes.titlesize': 13,
    'axes.labelsize': 12, 'legend.fontsize': 10, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.edgecolor': '#81909D',
    'axes.grid': True, 'grid.alpha': .18, 'grid.linewidth': .6,
    'lines.linewidth': 1.8, 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'svg.fonttype': 'none',
    'mathtext.fontset': 'dejavusans', 'axes.formatter.use_mathtext': True,
})
raw = loadmat(DATA / 'presentation_data.mat', simplify_cells=True)
cases = raw['cases']
base = sorted([c for c in cases if c['cells'] == 64],
              key=lambda c: (c['delta'], c['flow'] != 'couette'))
lookup = {(c['flow'], float(c['delta']), int(c['cells'])): c for c in cases}
figures = []


def delta_tag(delta):
    return f'{delta:g}'.replace('.', 'p')


def write_csv(name, rows):
    assert rows
    with (CSV / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def read_csv(path):
    with path.open() as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]


def finish(fig, name, title, caption, chinese):
    fig.suptitle(title, fontsize=17, fontweight='semibold', x=.06, ha='left', y=.98)
    fig.text(.06, .025, caption, fontsize=9.5, color='#526270', va='bottom')
    fig.tight_layout(rect=(.015, .07, .995, .92), pad=1.2)
    fig.savefig(OUT / f'{name}.png', dpi=240)
    fig.savefig(OUT / f'{name}.svg')
    fig.canvas.draw()
    groups=[]
    if name[:2] in ['05','06','07','08','09']:
        groups=[(delta_tag(delta),fig.axes[2*j:2*j+2]) for j,delta in enumerate([.1,1,10])]
    elif name.startswith('10_'):
        groups=[('part1',fig.axes[:2]),('part2',fig.axes[2:])]
    elif name.startswith('11_'):
        groups=[('flow',fig.axes[:1]),('differences',fig.axes[1:])]
    for suffix,axes in groups:
        bounds=Bbox.union([a.get_tightbbox(fig.canvas.get_renderer()) for a in axes])
        bounds=bounds.transformed(fig.dpi_scale_trans.inverted()).expanded(1.015,1.06)
        fig.savefig(PANELS/f'{name}_{suffix}.png',dpi=240,bbox_inches=bounds)
    figures.append(dict(name=name, title=chinese, caption=caption))
    plt.close(fig)


def metric(value, reference, h):
    d = value-reference
    absolute = float(np.sqrt(h*np.sum(d*d)))
    return absolute, 100*absolute/float(np.sqrt(h*np.sum(reference*reference)))


summary, profile_rows, residual_rows = [], [], []
for case in base:
    c, g = case['cis'], case['gsis']
    flow, delta = case['flow'], float(case['delta'])
    tag=delta_tag(delta)
    i={(1.,'couette'):1,(1.,'poiseuille'):2,(10.,'couette'):3,(10.,'poiseuille'):4,
       (.1,'couette'):14,(.1,'poiseuille'):15}[(delta,flow)]
    x = np.asarray(c['mesh']['center'])
    u, ug = c['cellMoments'][:, 2], g['cellMoments'][:, 2]
    ub = case['reference']['cellVelocity']
    h = c['config']['gap']/64
    absolute, relative = metric(ug, u, h)
    _, shear_relative = metric(g['shearStress'], c['shearStress'], h)
    _, c_ref_error = metric(u, ub, h)
    _, g_ref_error = metric(ug, ub, h)
    denominator = max(abs(u))
    rel_profile = 100*abs(ug-u)/denominator
    row = dict(flow=flow, delta=delta, cells=64, Nvx=32, Nvy=32,
               cisIterations=int(c['iterations']), gsisIterations=int(g['iterations']),
               iterationRatio=c['iterations']/g['iterations'],
               absoluteVelocityL2=absolute, maximumAbsoluteVelocityDifference=float(max(abs(ug-u))),
               relativeVelocityL2Percent=relative, relativeShearL2Percent=shear_relative,
               cisVsIndependentBGKPercent=c_ref_error, gsisVsIndependentBGKPercent=g_ref_error,
               independentBGKSelfRefinementPercent=100*case['reference']['selfRelativeL2'],
               cisFlowRate=float(c['flowRate']), gsisFlowRate=float(g['flowRate']),
               absoluteFlowRateDifference=abs(g['flowRate']-c['flowRate']),
               relativeFlowDifferencePercent=(100*abs(g['flowRate']-c['flowRate'])/abs(c['flowRate'])
                                              if flow == 'poiseuille' else float('nan')),
               cisOriginalResidual=float(c['diagnostics']['kinetic']),
               gsisOriginalResidual=float(g['diagnostics']['kinetic']),
               gsisCoupledResidual=float(g['diagnostics']['coupled']),
               cisRelativeMacroUpdate=float(c['convergence']['referenceMacroUpdate']),
               gsisRelativeMacroUpdate=float(g['convergence']['referenceMacroUpdate']))
    summary.append(row)
    rows = [dict(x=float(x[j]), cisUy=float(u[j]), gsisUy=float(ug[j]),
                 independentBGKUy=float(ub[j]), signedGSISMinusCIS=float(ug[j]-u[j]),
                 absoluteGSISMinusCIS=float(abs(ug[j]-u[j])),
                 globallyNormalizedAbsDifferencePercent=float(rel_profile[j]),
                 cisMinusIndependentBGK=float(u[j]-ub[j]),
                 gsisMinusIndependentBGK=float(ug[j]-ub[j])) for j in range(len(x))]
    write_csv(f'{flow}_delta{tag}_profiles.csv', rows)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.9))
    ax[0].plot(x, u, color=BLUE, label='CIS')
    ax[0].plot(x, ug, color=RED, linestyle='--', label='GSIS')
    ax[0].set(title='Velocity profile', xlabel='$x$', ylabel='$u_y$')
    ax[0].legend(loc='best')
    ax[1].plot(x, ug-u, color=RED)
    ax[1].axhline(0, color=BLACK, linewidth=.8)
    ax[1].set(title='Signed absolute difference', xlabel='$x$',
              ylabel=r'$u_y^{\rm GSIS}-u_y^{\rm CIS}$')
    ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax[2].plot(x, rel_profile, color=RED)
    ax[2].set(title=f'Relative $L^2$ difference: {relative:.4g}%', xlabel='$x$',
              ylabel=r'$100\,|\Delta u_y|/\max_x|u_y^{\rm CIS}|$ [%]')
    ax[2].set_ylim(bottom=0)
    finish(fig, f'{i:02d}_{flow}_delta{tag}_profiles',
           f'{FLOW[flow]} | $\delta={delta:g}$ | velocity and differences',
           'Nx=64, Nv=32x32, relative-update stop=1e-5. Pointwise differences use one global amplitude; no division by a zero crossing.',
           f'{FLOW[flow]} δ={delta:g}：速度剖面、绝对差与相对差')

write_csv('main_results.csv', summary)

# Save complete available histories, including NaNs at skipped diagnostics.
for case in base:
    for method in ['cis', 'gsis']:
        r = case[method]
        columns = list(r['historyColumns'])
        rows = [dict(zip(columns, (float(v) for v in line))) for line in r['history']]
        write_csv(f"{case['flow']}_delta{delta_tag(case['delta'])}_{method}_history.csv", rows)


def histories(kind):
    fig, axes = plt.subplots(len(base)//2, 2, figsize=(12.5, 3.5*(len(base)//2)+1))
    for ax, case in zip(axes.flat, base):
        c, g = case['cis'], case['gsis']
        hc, hg = c['history'], g['history']
        if kind == 'absolute_update':
            ax.semilogy(hc[:, 0], hc[:, 1], color=BLUE, label='CIS: kinetic update')
            ax.semilogy(hg[:, 0], hg[:, 1], color=RED, label='GSIS: kinetic update')
            ax.semilogy(hg[:, 0], hg[:, 7], color=GOLD, linestyle=':', label='GSIS: source update')
            ax.set_ylabel('Absolute moment update (RMS)')
        elif kind == 'relative_update':
            ax.semilogy(hc[:, 0], hc[:, -4], color=BLUE, label='CIS: kinetic monitor')
            ax.semilogy(hg[:, 0], hg[:, -4], color=RED, label='GSIS: mixed-source monitor')
            ax.axhline(1e-5, color=BLACK, linestyle='--', linewidth=1, label='Stopping tolerance')
            ax.set_ylabel('Relative cell-macro update')
        else:
            denominator = float(hc[0, 2]) if kind == 'relative_residual' else 1.
            for hist, col, color, style, label in [
                (hc, 2, BLUE, '-', 'CIS: original equation'),
                (hg, 2, RED, '-', 'GSIS: original CIS equation'),
                (hg, 8, GOLD, '--', 'GSIS: coupled equation'),
            ]:
                mask = np.isfinite(hist[:, col])
                ax.semilogy(hist[mask, 0], hist[mask, col]/denominator,
                            color=color, linestyle=style, label=label)
            ax.set_ylabel(r'$R/R^{\rm CIS}_1$' if denominator != 1 else 'Absolute equation residual $R$')
        ax.set_title(f"{FLOW[case['flow']]}  |  $\delta={case['delta']:g}$")
        ax.set_xlabel('Outer iteration')
        ax.legend(loc='best', fontsize=9)
    settings = {
        'absolute_update': ('05_absolute_updates', 'Absolute iteration updates',
                            'Max over four moment RMS updates. These are changes between iterations, not errors against a reference solution.',
                            '相邻迭代的绝对宏观量更新'),
        'relative_update': ('06_relative_updates', 'Relative macro updates used for stopping',
                            'Relative L2 changes of cell-integrated physical rho, T and speed magnitude; max of three. GSIS monitors the mixed source.',
                            '实际用于停机的相对宏观量更新'),
        'absolute_residual': ('07_absolute_residuals', 'Absolute equation residuals',
                              'Weighted kinetic-DG weak residual norm. Equation residuals are diagnostics, not the 1e-5 stopping criterion.',
                              '绝对方程残差：原 CIS 方程与耦合方程'),
        'relative_residual': ('08_relative_residuals', 'Relative equation residuals',
                              'All curves in each panel share the CIS first-iteration residual as denominator. This is not the stopping metric.',
                              '相对方程残差：统一除以各算例 CIS 首轮残差'),
    }
    finish(fig, *settings[kind])


for kind in ['absolute_update', 'relative_update', 'absolute_residual', 'relative_residual']:
    histories(kind)

spatial = []
fig, axes = plt.subplots(len(base)//2, 2, figsize=(12.5, 3.5*(len(base)//2)+1))
for ax, baseline in zip(axes.flat, base):
    group = [lookup[(baseline['flow'], float(baseline['delta']), n)] for n in [32, 64, 128]]
    nc = [int(c['cis']['iterations']) for c in group]
    ng = [int(c['gsis']['iterations']) for c in group]
    ax.plot([32,64,128], nc, 'o-', color=BLUE, label='CIS')
    ax.plot([32,64,128], ng, 's--', color=RED, markerfacecolor='white', label='GSIS')
    for x, y, z in zip([32,64,128], nc, ng):
        ax.annotate(str(y), (x,y), xytext=(0,8), textcoords='offset points', ha='center', color=BLUE)
        if y != z:
            ax.annotate(str(z), (x,z), xytext=(0,-17), textcoords='offset points', ha='center', color=RED)
    ax.set_xscale('log', base=2);ax.set_xticks([32,64,128], ['32','64','128'])
    ax.set(xlabel='Spatial cells $N_x$', ylabel='Iterations',
           title=f"{FLOW[baseline['flow']]}  |  $\delta={baseline['delta']:g}$")
    ax.set_ylim(0, max(nc+ng)*1.25);ax.legend(loc='best')
    for item in group:
        c,g=item['cis'],item['gsis'];h=c['config']['gap']/item['cells']
        abs_error,rel_error=metric(g['cellMoments'][:,2],c['cellMoments'][:,2],h)
        spatial.append(dict(flow=item['flow'],delta=float(item['delta']),cells=int(item['cells']),
                            cisIterations=int(c['iterations']),gsisIterations=int(g['iterations']),
                            iterationRatio=c['iterations']/g['iterations'],
                            absoluteVelocityL2=abs_error,relativeVelocityL2Percent=rel_error,
                            omega=float(np.min(g['synthetic']['omega'])),
                            stopCriterion='reference',stopTolerance=1e-5,Nvx=32,Nvy=32))
finish(fig,'09_spatial_iterations','Iteration counts under spatial refinement',
       'Nv=32x32, same 1e-5 stopping. The GSIS mixing coefficient changes with cell size; this is not a fixed-omega study.',
       '空间分辨率与迭代次数：统一阈值和速度网格')
write_csv('spatial_iterations_and_differences.csv',spatial)

fig, axes = plt.subplots(1,4,figsize=(18,4.9))
for ax,delta in zip(axes[:3],[.1,1,10]):
    for flow,color in [('couette',BLUE),('poiseuille',RED)]:
        group=[r for r in spatial if r['flow']==flow and r['delta']==delta]
        ax.loglog([r['cells'] for r in group],[r['relativeVelocityL2Percent'] for r in group],
                  'o-',color=color,label=FLOW[flow])
    ax.set(title=f'GSIS vs CIS | $\delta={delta}$',xlabel='Spatial cells $N_x$',ylabel='Relative velocity $L^2$ difference [%]')
    ax.set_xscale('log',base=2);ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticks([32,64,128],['32','64','128']);ax.legend(loc='best')
for delta,color in [(.1,GOLD),(1,BLUE),(10,RED)]:
    group=[r for r in spatial if r['flow']=='couette' and r['delta']==delta]
    axes[3].plot([r['cells'] for r in group],[r['omega'] for r in group], 'o-',color=color,label=f'$\delta={delta:g}$')
axes[3].set(title='Actual GSIS mixing',xlabel='Spatial cells $N_x$',ylabel=r'$\omega=\min(1,2\delta\Delta x/\sqrt{\pi})$')
axes[3].set_xscale('log',base=2);axes[3].set_xticks([32,64,128],['32','64','128']);axes[3].legend()
finish(fig,'10_spatial_differences_and_mixing','Spatial refinement: differences and mixing',
       'Changes in iteration counts and solution differences occur together with the prescribed cell-size-dependent mixing.',
       '空间加密后的解差异及混合系数')

literature_dir=ROOT/'results/channel_1d2v/literature_poiseuille'
old=read_csv(literature_dir/'comparison.csv')
velocity_check=read_csv(literature_dir/'velocity_check_summary.csv')
velocity_rows=read_csv(literature_dir/'velocity_check.csv')
assert len(old)==6 and len(velocity_rows)==30
expected_paper=[2.0327,1.6022,1.5392,1.5957,1.9922,2.7673]
assert [r['paperBGK'] for r in old]==expected_paper
paper={r['delta']:r['paperBGK'] for r in old}
for row in old:
    row['signedAbsoluteDifference']=row['fineQp']-row['paperBGK']
    row['absoluteDifference']=abs(row['signedAbsoluteDifference'])
    row['absoluteRefinementChange']=abs(row['fineQp']-row['coarseQp'])
write_csv('historical_literature_comparison.csv',old)
write_csv('historical_velocity_check.csv',velocity_rows)
write_csv('historical_velocity_check_summary.csv',velocity_check)
current=[]
for case in base:
    if case['flow']!='poiseuille':continue
    for method in ['cis','gsis']:
        r=case[method];cfg=r['config'];delta=case['delta']
        qp=r['flowRate']/(cfg['acceleration']*cfg['gap']**2)
        current.append(dict(delta=float(delta),method=method.upper(),cells=64,Nvx=32,Nvy=32,
                            stopCriterion='reference',stopTolerance=1e-5,Qp=float(qp),paperBGK=paper[delta],
                            signedAbsoluteDifference=float(qp-paper[delta]),
                            signedRelativeDifferencePercent=100*(qp/paper[delta]-1)))
write_csv('current_literature_comparison.csv',current)

fig,ax=plt.subplots(1,3,figsize=(16,5.5))
d=np.array([r['delta'] for r in old]);qp=np.array([r['fineQp'] for r in old]);p=np.array(expected_paper)
ax[0].semilogx(d,qp,'o-',color=BLUE,label='Historical fine CIS')
ax[0].semilogx(d,p,'s',color=BLACK,markerfacecolor='white',label='Published BGK [36]')
ax[1].semilogx(d,qp-p,'o-',color=BLUE,label='Historical fine CIS')
ax[2].semilogx(d,100*(qp/p-1),'o-',color=BLUE,label='Historical fine CIS')
rounding=100*.00005/p
ax[2].fill_between(d,-rounding,rounding,color='#D5DBE0',alpha=.7,label='4-decimal rounding only')
for method,color,marker in [('CIS',GOLD,'D'),('GSIS',RED,'^')]:
    group=[r for r in current if r['method']==method]
    ax[0].semilogx([r['delta'] for r in group],[r['Qp'] for r in group],marker,
                  markersize=8,color=color,markerfacecolor='white',label=f'Current {method}, 32x32 velocities')
    ax[1].semilogx([r['delta'] for r in group],[r['signedAbsoluteDifference'] for r in group],
                  marker,markersize=8,color=color,markerfacecolor='white',label=f'Current {method}')
    ax[2].semilogx([r['delta'] for r in group],[r['signedRelativeDifferencePercent'] for r in group],
                  marker,markersize=8,color=color,markerfacecolor='white',label=f'Current {method}')
ax[0].set(title='Poiseuille flow coefficient',xlabel='$\delta$',ylabel='$Q_p$')
ax[1].set(title='Signed absolute difference',xlabel='$\delta$',ylabel=r'$Q_p-Q_{p,\rm paper}$')
ax[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax[2].set(title='Signed relative difference',xlabel='$\delta$',ylabel=r'$100(Q_p/Q_{p,\rm paper}-1)$ [%]')
ax[1].axhline(0,color=BLACK,linestyle=':',linewidth=1)
ax[2].axhline(0,color=BLACK,linestyle=':',linewidth=1)
for a in ax:a.legend(fontsize=9,loc='best')
finish(fig,'11_literature_Qp','Published BGK benchmark | infinite parallel plates',
       'Sharipov (1999), Table 1, a/b=0. Historical: Nx256/Nv512x8, equation stop 1e-10. Current: Nx64/Nv32x32, update stop 1e-5.',
       'Poiseuille 文献流量对照：历史精细 CIS 与当前 CIS/GSIS')

fig,ax=plt.subplots(1,2,figsize=(13,5.5))
ax[0].loglog(d,[r['refinementPercent'] for r in old],'o-',color=BLUE)
ax[0].set(title='Historical combined resolution check',xlabel='$\delta$',ylabel='Relative change of $Q_p$ [%]')
for col,color,label,style in [
    ('Nvx256to512',BLUE,'Normal velocity: 256 to 512','o-'),
    ('Nvx512to1024',RED,'Normal velocity: 512 to 1024','s-'),
    ('cutoff8to10_atNvx1024',GOLD,'Cutoff: 8 to 10 (Nvx=1024)','^-'),
    ('Nvy8to16_atNvx512',BLACK,'Tangential velocity: 8 to 16','d--'),
]:
    ax[1].loglog(d,[max(100*r[col],1e-14) for r in velocity_check],style,color=color,label=label)
ax[1].set(title='One-setting-at-a-time velocity checks',xlabel='$\delta$',ylabel='Relative change of $Q_p$ [%]')
ax[1].legend(fontsize=9,loc='best')
finish(fig,'12_literature_resolution','Resolution checks for the historical literature benchmark',
       'Left: Nx128/Nv256x8 to Nx256/Nv512x8. Right: Nx256 fixed, equation stop 1e-10. Values below 1e-14% are drawn at 1e-14%.',
       '历史文献对照的空间/速度分辨率核查')

velocity_iterations=[]
fig,ax=plt.subplots(figsize=(10,5.5))
colors=[BLUE,RED,GOLD,BLACK,'#5E8150','#79588E']
for delta,color in zip(d,colors):
    group=sorted([r for r in velocity_rows if r['delta']==delta and r['Nvy']==8 and r['normalCutoff']==8],
                 key=lambda r:r['Nvx'])
    counts=[int(r['iterations']) for r in group]
    assert len(group)==3 and len(set(counts))==1
    ax.semilogy([r['Nvx'] for r in group],counts,'o-',color=color,
                label=f'$\delta={delta:g}$: {counts[0]} iterations')
    velocity_iterations.append(dict(delta=delta,Nvx256=counts[0],Nvx512=counts[1],Nvx1024=counts[2],
                                    cells=256,Nvy=8,stopCriterion='equation',stopTolerance=1e-10))
ax.set_xscale('log',base=2);ax.set_xticks([256,512,1024],['256','512','1024'])
ax.set(xlabel='Normal velocity points $N_{v_x}$',ylabel='CIS iterations')
ax.set_yticks([15,30,48,89,259,721],['15','30','48','89','259','721'])
ax.legend(ncol=3,loc='lower center',bbox_to_anchor=(.5,1.02),fontsize=10)
finish(fig,'13_historical_velocity_iterations','Historical CIS iterations under velocity refinement',
       'Historical CIS only: Nx256, Nvy8, cutoff8, unit drive, absolute equation stop 1e-10. All variants start from H=0.',
       '历史 CIS：速度网格加密后的迭代次数')
write_csv('historical_velocity_iterations.csv',velocity_iterations)
figures.sort(key=lambda item:int(item['name'].split('_',1)[0]))

provenance=dict(date='2026-11-09',data_file=str(DATA/'presentation_data.mat'),
    source_paper='/Users/hanchenxi/Desktop/Non-isothermal gas flow through rectangular microchannel.pdf',
    paper_table='Sharipov (1999), Table 1, printed page 398, a/b=0, BGK work [36]',
    rounding_band='Assumed rounding to four decimal places (+/-0.00005 in Qp); not uncertainty reported by the paper.',
    paper_delta='delta = gap*p/(mu*c0), c0=sqrt(2*k*T/m); Kn=sqrt(pi)/(2*delta)',
    code_delta='paper_delta=code_delta*gap; all presented cases have gap=1',
    flow_normalization='Qp=flowRate/(acceleration*gap^2); linear amplitudes are normalized before comparison',
    absolute_difference='Signed pointwise delta_u=GSIS-CIS; absolute L2=sqrt(dx*sum(delta_u^2))',
    relative_difference='Global L2=norm(delta_u)/norm(CIS). Profile plots use abs(delta_u)/max(abs(CIS)), avoiding zero-crossing division.',
    absolute_residual='Mass-matrix-inverse and Maxwell-velocity-weighted kinetic-DG weak residual norm',
    relative_residual='R_k/R_CIS_1: all curves share the CIS first-iteration residual for the same flow and delta',
    stopping='Max relative updates of cell-integrated physical density, temperature and speed magnitude <1e-5. CIS monitors kinetic moments; GSIS monitors mixed source.',
    caution='The historical literature table uses CIS, fine velocity grids and equation stopping 1e-10. Do not label it as a current 32x32 GSIS validation.',
    figures=figures)
(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2,ensure_ascii=False)+'\n')


def table(rows,columns):
    head=''.join(f'<th>{html.escape(label)}</th>' for _,label in columns)
    body=''
    for row in rows:
        cells=''
        for key,_ in columns:
            v=row[key]
            if isinstance(v,(float,np.floating)):
                v='—' if np.isnan(v) else f'{v:.7g}'
            cells+=f'<td>{html.escape(str(v))}</td>'
        body+=f'<tr>{cells}</tr>'
    return f'<div class="scroll"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>1D2V CIS / GSIS 算例图表</title><style>
body{font:16px/1.7 -apple-system,BlinkMacSystemFont,"PingFang SC",sans-serif;color:#243447;background:#f5f7fa;margin:0}
main{max-width:1180px;margin:auto;padding:34px 28px}h1{font-size:30px;margin:0 0 8px}h2{margin-top:36px}
p{max-width:1000px}section{background:white;border:1px solid #dde4eb;border-radius:10px;padding:22px;margin:22px 0}
img{width:100%;height:auto}table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:14px}
th,td{padding:9px 11px;text-align:right;border-bottom:1px solid #dde4eb;white-space:nowrap}th{background:#edf2f7}th:first-child,td:first-child{text-align:left}
a{color:#2568a6}.scroll{overflow-x:auto}.meta{color:#526270}nav a{display:inline-block;margin:0 18px 10px 0}.note{border-left:4px solid #ad811e;padding-left:16px}
</style><main><h1>1D2V CIS / GSIS：差异、残差与分辨率</h1><p class="meta">2026-11-09 · PNG 适合直接插入 PPT；SVG 可缩放。所有 CSV 保留原始精度。</p>
<nav><a href="#main">主结果</a><a href="#resolution">分辨率</a><a href="#paper">文献对照</a><a href="#figures">全部图片</a><a href="provenance.json">定义和数据来源</a></nav>
<section id="main"><h2>当前主算例</h2><p>空间 P1，64 单元；速度 32×32，法向半区 Gauss-Legendre、流向 Gauss-Hermite；gap=1，初值 H=0。
Couette：壁速 [-0.01,0.01]、加速度 0；Poiseuille：静止壁、加速度 0.01。两壁等温完全漫反射，稳定参数 [1,δ,δ,δ]。
CIS/GSIS 统一采用宏观量相对更新 &lt;10⁻⁵ 停机。</p>'''
page+=table(summary,[('flow','流动'),('delta','δ'),('cisIterations','CIS 步数'),('gsisIterations','GSIS 步数'),
                     ('iterationRatio','迭代数比'),('absoluteVelocityL2','速度绝对 L2 差'),('relativeVelocityL2Percent','速度相对差 [%]'),
                     ('relativeShearL2Percent','剪应力相对差 [%]'),('relativeFlowDifferencePercent','流量相对差 [%]')])
page+='<p>迭代数比不等于运行时间加速比。以上差异以同网格 CIS 为基准，不是相对精确解的误差。</p>'
page+='<p class="note">绝对差 Δu=GSIS−CIS；全局相对 L2 差 =‖Δu‖₂/‖CIS‖₂。剖面图中的局部曲线使用 |Δu(x)|/max|CIS|，避免 Couette 中心零点导致除零。'
page+='绝对方程残差按原 kinetic-DG 弱残差的加权范数计算；相对残差统一为 Rₖ/R₁ᶜᴵˢ。这些残差均不参与当前停机。</p>'
page+='</section>'
page+='<section id="resolution"><h2>不同空间分辨率的迭代次数</h2><p>速度固定 32×32，停机统一为 10⁻⁵。CIS 的次数不变；δ=1 的 GSIS 也不变，但 δ=10 的 GSIS 次数随加密增加。'
page+='同时 ω=min(1,2δΔx/√π) 随网格变化，因此这不是固定混合系数的研究，不能宣称当前 GSIS 迭代次数与网格无关。</p>'
page+=table(spatial,[('flow','流动'),('delta','δ'),('cells','空间单元'),('cisIterations','CIS'),('gsisIterations','GSIS'),
                     ('relativeVelocityL2Percent','速度相对差 [%]'),('omega','ω')])
page+='</section><section id="paper"><h2>上次的 Poiseuille 文献对照</h2><p>Sharipov, <i>Non-isothermal gas flow through rectangular microchannels</i>, '
page+='J. Micromech. Microeng. 9 (1999), 394–401，Table 1（印刷页 398）。这里只用平行板 a/b=0 的 BGK work [36] 列。'
page+='δ=gap·p/(μc₀)，Kn=√π/(2δ)，Qₚ=flowRate/(aᵧ·gap²)。</p>'
page+='<p>以下是历史精细 CIS：Nx=256、速度 512×8、单位驱动、绝对方程判据 10⁻¹⁰；与当前 32×32 GSIS 结果分别报告。</p>'
page+=table(old,[('delta','δ'),('fineQp','历史 CIS Qp'),('paperBGK','文献 BGK Qp'),('signedAbsoluteDifference','有符号绝对差'),('bgkDifferencePercent','有符号相对差 [%]'),
                ('refinementPercent','联合加密变化 [%]'),('fineIterations','历史 CIS 步数')])
page+='<p>历史联合加密为 Nx128/Nv256×8 → Nx256/Nv512×8；速度求积单独加密检查见 data/historical_velocity_check.csv。'
page+='文献 BGK 数字保留四位小数，其舍入区间不足以解释所有 δ 下的差异，不能把剩余偏差统称为打印精度误差。</p>'
page+='<p>图中灰带仅表示假设 Qp 四舍五入到四位小数的区间（±0.00005），不是论文报告的误差条。</p>'
page+=table(velocity_iterations,[('delta','δ'),('Nvx256','Nvx=256 步数'),('Nvx512','Nvx=512 步数'),('Nvx1024','Nvx=1024 步数')])
page+=table(current,[('delta','δ'),('method','当前方法'),('Qp','当前 Qp'),('paperBGK','文献 BGK Qp'),('signedRelativeDifferencePercent','有符号相对差 [%]')])
page+='</section><h2 id="figures">可用于 PPT 的图片</h2>'
for fig in figures:
    name=fig['name']
    page+=f'<section><h2>{html.escape(fig["title"])}</h2><p><a href="{name}.png">PNG 原图</a> · <a href="{name}.svg">SVG 矢量图</a></p>'
    page+=f'<a href="{name}.png"><img loading="lazy" src="{name}.png" alt="{html.escape(fig["title"])}"></a>'
    page+=f'<p class="meta">{html.escape(fig["caption"])}</p></section>'
page+='</main></html>'
(OUT/'index.html').write_text(page)

# A contact sheet is a navigation aid; use individual images in the slides.
from PIL import Image, ImageOps, ImageDraw
tiles=[]
for fig in figures:
    with Image.open(OUT/f"{fig['name']}.png") as im:
        tile=Image.new('RGB',(840,560),'white')
        thumb=ImageOps.contain(im.convert('RGB'),(820,510))
        tile.paste(thumb,((840-thumb.width)//2,35+(510-thumb.height)//2))
        ImageDraw.Draw(tile).text((20,12),fig['name'],fill=BLACK)
        tiles.append(tile)
sheet=Image.new('RGB',(840*3,560*((len(tiles)+2)//3)),'#E9EEF3')
for i,tile in enumerate(tiles):sheet.paste(tile,((i%3)*840,(i//3)*560))
sheet.save(OUT/'overview.png')

archive=OUT.parent/'channel_figures_20261109.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
    for path in sorted(OUT.rglob('*')):
        if path.is_file():z.write(path,path.relative_to(OUT.parent))
print('Saved',len(figures),'figures as PNG + SVG, CSV tables, index.html and ZIP.')
print('Output:',OUT)
print('Main results:')
for row in summary:print(row['flow'],row['delta'],row['absoluteVelocityL2'],row['relativeVelocityL2Percent'])
print('Current literature differences:',current)
