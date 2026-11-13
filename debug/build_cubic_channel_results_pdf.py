"""Six-page retest of every brief-report metric after cubic velocity mapping."""
from pathlib import Path

# Reuse layout helpers, fonts and readers, without running the detailed report.
source=Path(__file__).with_name('build_channel_results_pdf.py').read_text()
source=source.split('# Supplementary diagnostics, recomputed only from saved fields.')[0]
source=source.replace("FIG=ROOT/'output/channel_figures_20261109'","FIG=ROOT/'output/channel_cubic_figures_20261109'")
source=source.replace("TMP=ROOT/'tmp/pdfs/channel_results_20261109'","TMP=ROOT/'tmp/pdfs/channel_cubic_results_20261109'")
source=source.replace("TARGET=OUT/'1D2V_CIS_GSIS_指标与逐图说明_20261109.pdf'",
                      "TARGET=OUT/'1D2V_CIS_GSIS_全指标与停机敏感性_20261109.pdf'")
source=source.replace("ROOT/'results/channel_presentation_20261109/presentation_data.mat'",
                      "ROOT/'results/channel_cubic_presentation_20261109/presentation_data.mat'")
source=source.replace('全部指标、算例结果与逐图说明','三次映射全指标复测')
source=source.replace('所有速度与流量均为无量纲量；当前主算例的停机标准为相对宏观更新 < 1e-5。',
                      '主结果停机阈值 1e-5；敏感性另测 1e-7、1e-9。所有速度与流量均为无量纲量。')
exec(compile(source,str(Path(__file__)),'exec'))

import matplotlib.pyplot as plt
import zipfile
from shutil import copyfile

oldmain=csv_rows(ROOT/'output/channel_figures_20261109/data/main_results.csv')
oldspatial=csv_rows(ROOT/'output/channel_figures_20261109/data/spatial_iterations_and_differences.csv')
old={(r['flow'],float(r['delta'])):r for r in oldmain}
oldsaved=loadmat(ROOT/'results/channel_presentation_20261109/presentation_data.mat',simplify_cells=True)
oldcases={(r['flow'],float(r['delta'])):r for r in oldsaved['cases'] if r['cells']==64}
flowrows=csv_rows(ROOT/'results/cubic_poiseuille_flowrates32/cubic_results.csv')
oldflows=csv_rows(ROOT/'results/cubic_poiseuille_flowrates32/summary.csv')
oldflows={float(r['delta']):r for r in oldflows if r['grid']=='original'}
velocity=csv_rows(ROOT/'results/channel_cubic_presentation_20261109/velocity_refinement.csv')
stopping=csv_rows(ROOT/'results/cubic_stopping_residuals32/summary.csv')
stop={(r['flow'],float(r['delta']),float(r['tolerance'])):r for r in stopping}
assert len(stop)==18
copyfile(ROOT/'results/channel_cubic_presentation_20261109/velocity_refinement.csv',DATA/'current_velocity_refinement.csv')
copyfile(ROOT/'results/cubic_poiseuille_flowrates32/cubic_results.csv',DATA/'all_six_flowrates.csv')
copyfile(ROOT/'results/cubic_poiseuille_flowrates32/strict_cis_1e9.csv',DATA/'strict_cis_flowrates_1e9.csv')
copyfile(ROOT/'results/cubic_stopping_residuals32/summary.csv',DATA/'stopping_sensitivity_all_six_cases.csv')


def flowlabel(flow,delta):return f'{"C" if flow=="couette" else "P"} / {delta:g}'
def pair(a,b,digits=4):return f'{f(a,digits)} → {f(b,digits)}'
def write_csv(name,rows):
    with (DATA/name).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


# Keep all diagnostics and old/new comparisons as CSV, without lengthening PDF.
changes=[];diagnostics=[]
for case in cases:
    flow,delta=key(case);r=summary[(flow,delta)];previous=old[(flow,delta)]
    row={'flow':flow,'delta':delta}
    for name in ['cisIterations','gsisIterations','absoluteVelocityL2','maximumAbsoluteVelocityDifference',
                 'relativeVelocityL2Percent','relativeShearL2Percent','relativeFlowDifferencePercent',
                 'cisOriginalResidual','gsisOriginalResidual','gsisCoupledResidual',
                 'cisRelativeMacroUpdate','gsisRelativeMacroUpdate']:
        row['old_'+name]=float(previous[name]);row['new_'+name]=float(r[name])
    changes.append(row)
for case in saved['cases']:
    flow,delta=key(case);den=case['cis']['history'][0,2]
    for method in ['cis','gsis']:
        r=case[method];d=r['diagnostics'];q=r['cellMoments'];stress=r['shearStress']
        row=dict(flow=flow,delta=delta,cells=int(case['cells']),method=method,
            iterations=int(r['iterations']),absoluteKineticUpdate=float(r['history'][-1,1]),
            absoluteSourceUpdate=float(r['history'][-1,7]) if method=='gsis' else float(r['history'][-1,1]),
            relativeMacroUpdate=float(r['convergence']['referenceMacroUpdate']),
            originalResidual=float(d['kinetic']),coupledResidual=float(d['coupled']),
            relativeOriginalResidual=float(d['kinetic']/den),relativeCoupledResidual=float(d['coupled']/den),
            wallResidual=float(d['wall']),wallMassFlux=float(d['wallMassFlux']),massError=float(d['massError']),
            collisionInvariantDefect=float(d['collisionInvariantDefect']),sourceMomentGap=float(d['sourceMomentGap']),
            maxRhoPerturbation=float(np.max(abs(q[:,0]))),maxNormalVelocity=float(np.max(abs(q[:,1]))),
            maxTemperaturePerturbation=float(np.max(abs(q[:,3]))),meanShear=float(np.mean(stress)),
            maxAbsShear=float(np.max(abs(stress))),
            relativeShearSpanPercent=float(100*np.ptp(stress)/abs(np.mean(stress))) if flow=='couette' else float('nan'))
        for j,name in enumerate(['massBalance','streamwiseMomentumBalance','normalMomentumBalance','energyBalance']):
            row[name]=float(d['balance'][j])
        if method=='gsis':
            v=r['velocity'];hw=r['H'][[0,-1],:];w=v['weight']
            rho=hw@w;temp=hw@(w*(v['vx']**2+v['vy']**2-1));uy=hw@(w*v['vy'])
            target=np.column_stack((rho+temp,np.zeros(2),uy,temp))
            hdg=r['synthetic']['hdg']
            row['macroBoundaryTraceError']=float(np.max(abs(target-hdg['Uhat'][[0,-1],:])))
            for name in ['localResidual','traceResidual','fluxContinuityResidual','conservationResidual']:
                row['hdg_'+name]=float(hdg[name])
        else:
            row.update(macroBoundaryTraceError=float('nan'))
            for name in ['localResidual','traceResidual','fluxContinuityResidual','conservationResidual']:
                row['hdg_'+name]=float('nan')
        diagnostics.append(row)
write_csv('old_vs_cubic_main_metrics.csv',changes)
write_csv('all_spatial_solver_diagnostics.csv',diagnostics)

# New compact scientific figures are also available separately for slides.
plt.rcParams.update({'font.size':8.5,'axes.titlesize':9.5,'axes.labelsize':8.5,
                     'legend.fontsize':7.2,'lines.linewidth':1.5})
blue='#2568A6';red='#BF4F39';gold='#AD811E';grey='#7D8993'


def savefig(fig,name):
    path=FIG/f'{name}.png';fig.savefig(path,dpi=240,facecolor='white');fig.savefig(FIG/f'{name}.svg')
    plt.close(fig);return path


fig,axes=plt.subplots(3,2,figsize=(11.8,5.5),constrained_layout=True)
for i,delta in enumerate([.1,1,10]):
    for j,flow in enumerate(['couette','poiseuille']):
        case=lookup[(flow,delta)];ax=axes[i,j];x=case['cis']['mesh']['center']
        ax.plot(x,oldcases[(flow,delta)]['cis']['cellMoments'][:,2],':',color=grey,label='Old-grid CIS')
        ax.plot(x,case['cis']['cellMoments'][:,2],color=blue,label='Cubic CIS')
        ax.plot(x,case['gsis']['cellMoments'][:,2],'--',color=red,label='Cubic GSIS')
        ax.set(title=f'{FLOW[flow]} | delta={delta:g}',xlabel='$x$',ylabel='$u_y$');ax.grid(alpha=.15)
        if i==j==0:ax.legend(ncol=3,loc='upper left')
profiles=savefig(fig,'16_cubic_profiles_with_old_grid')

fig,axes=plt.subplots(2,4,figsize=(12,5.4),constrained_layout=True)
for i,flow in enumerate(['couette','poiseuille']):
    case=lookup[(flow,10.)];hc=case['cis']['history'];hg=case['gsis']['history']
    for j,kind in enumerate(['absolute update','relative update','absolute residual','relative residual']):
        ax=axes[i,j]
        if j==0:
            ax.semilogy(hc[:,0],hc[:,1],color=blue,label='CIS kinetic')
            ax.semilogy(hg[:,0],hg[:,1],color=red,label='GSIS kinetic')
            ax.semilogy(hg[:,0],hg[:,7],':',color=gold,label='GSIS source')
            ylabel='RMS update'
        elif j==1:
            ax.semilogy(hc[:,0],hc[:,-4],color=blue,label='CIS monitor')
            ax.semilogy(hg[:,0],hg[:,-4],color=red,label='GSIS source monitor')
            ax.axhline(1e-5,color=grey,ls='--',lw=.8);ylabel='Relative macro update'
        else:
            denominator=hc[0,2] if j==3 else 1.
            for h,col,color,style,lab in [(hc,2,blue,'-','CIS original'),(hg,2,red,'-','GSIS original'),(hg,8,gold,'--','GSIS coupled')]:
                mask=np.isfinite(h[:,col]);ax.semilogy(h[mask,0],h[mask,col]/denominator,style,color=color,label=lab)
            ylabel='$R/R_{CIS,1}$' if j==3 else '$R$'
        ax.set(title=f'{FLOW[flow]} | {kind}',xlabel='Iteration',ylabel=ylabel)
        ax.grid(alpha=.15);ax.legend(loc='best',fontsize=6.3)
histories=savefig(fig,'17_cubic_all_four_convergence_metrics_delta10')

fig,axes=plt.subplots(1,2,figsize=(11.8,3),constrained_layout=True)
d=np.array([float(r['delta']) for r in flowrows]);paper=np.array([float(r['paperBGKQp']) for r in flowrows])
for method,color,marker in [('cis',blue,'o'),('gsis',red,'^')]:
    q=np.array([float(r[method+'Qp']) for r in flowrows]);oq=np.array([float(oldflows[x][method+'Qp']) for x in d])
    axes[0].semilogx(d,q,marker+'-',color=color,label=f'Cubic {method.upper()}')
    axes[1].semilogx(d,100*(q/paper-1),marker+'-',color=color,label=f'Cubic {method.upper()}')
    axes[1].semilogx(d,100*(oq/paper-1),':',color=color,alpha=.45,label=f'Old {method.upper()}')
axes[0].semilogx(d,paper,'ks',mfc='white',label='Published BGK [36]')
axes[1].axhline(0,color=grey,ls=':',lw=.7)
for ax in axes:ax.grid(alpha=.15);ax.legend(fontsize=7)
axes[0].set(xlabel='$\delta$',ylabel='$Q_p$',title='Flow coefficient, all six deltas')
axes[1].set(xlabel='$\delta$',ylabel='Signed difference to paper [%]',title='Old grid and cubic grid')
flowfig=savefig(fig,'18_cubic_all_six_flowrates_old_vs_new')

fig,axes=plt.subplots(1,2,figsize=(11.8,2.8),constrained_layout=True)
for flow,color in [('couette',blue),('poiseuille',red)]:
    newgroup=sorted([r for r in spatial if r['flow']==flow and float(r['delta'])==10],key=lambda r:float(r['cells']))
    oldgroup=sorted([r for r in oldspatial if r['flow']==flow and float(r['delta'])==10],key=lambda r:float(r['cells']))
    x=[int(r['cells']) for r in newgroup]
    axes[0].loglog(x,[float(r['relativeVelocityL2Percent']) for r in newgroup],'o-',color=color,label=f'Cubic {FLOW[flow]}')
    axes[0].loglog(x,[float(r['relativeVelocityL2Percent']) for r in oldgroup],':',color=color,label=f'Old {FLOW[flow]}')
    axes[1].plot(x,[int(r['gsisIterations']) for r in newgroup],'o-',color=color,label=f'GSIS {FLOW[flow]}')
axes[0].set(title='delta=10 | GSIS-CIS velocity difference',ylabel='Relative L2 difference [%]')
axes[1].set(title='delta=10 | GSIS counts: old = cubic',ylabel='Iterations')
for ax in axes:
    ax.set_xscale('log',base=2);ax.set_xticks([32,64,128],['32','64','128']);ax.set_xlabel('Spatial cells');ax.legend(fontsize=7);ax.grid(alpha=.15)
spacefig=savefig(fig,'19_cubic_spatial_refinement_old_vs_new')

fig,axes=plt.subplots(1,2,figsize=(11.8,3),constrained_layout=True)
for delta,color in [(.1,blue),(1,gold),(10,red)]:
    for ax,flow in zip(axes,['couette','poiseuille']):
        group=sorted([r for r in velocity if r['flow']==flow and float(r['delta'])==delta and int(r['Nvy'])==32],key=lambda r:int(r['Nvx']))
        ax.loglog([int(r['Nvx']) for r in group],[float(r['cisVelocityVs128x64Percent']) for r in group],'o-',color=color,label=f'delta={delta:g}')
        ax.set(title=FLOW[flow],xlabel='Normal velocity points',ylabel='CIS difference to 128x64 [%]')
        ax.set_xticks([32,64,128],['32','64','128']);ax.legend();ax.grid(alpha=.15)
savefig(fig,'20_cubic_current_velocity_refinement')

fig,axes=plt.subplots(3,2,figsize=(11.8,8),constrained_layout=True)
for i,delta in enumerate([.1,1,10]):
    for j,flow in enumerate(['couette','poiseuille']):
        ax=axes[i,j];group=[stop[(flow,delta,t)] for t in [1e-5,1e-7,1e-9]]
        x=[float(r['tolerance']) for r in group]
        for column,color,marker,label in [
            ('cisRelativeOriginalResidual',blue,'o','CIS original'),
            ('gsisRelativeOriginalResidual',red,'s','GSIS original'),
            ('gsisRelativeCoupledResidual',gold,'^','GSIS coupled')]:
            ax.loglog(x,[float(r[column]) for r in group],marker+'-',color=color,label=label)
        ax.invert_xaxis();ax.set(title=f'{FLOW[flow]} | delta={delta:g}',xlabel='Relative macro stopping tolerance',ylabel='Relative equation residual')
        ax.set_xticks([1e-5,1e-7,1e-9]);ax.legend(fontsize=7);ax.grid(alpha=.15)
savefig(fig,'21_cubic_stopping_sensitivity_all_six_cases')

vrows=[]
for delta in [.1,1,10]:
    c=next(r for r in velocity if r['flow']=='couette' and float(r['delta'])==delta and int(r['Nvx'])==32)
    q=next(r for r in velocity if r['flow']=='poiseuille' and float(r['delta'])==delta and int(r['Nvx'])==32)
    vrows.append([f(delta),f(c['cisVelocityVs128x64Percent'],5),f(q['cisVelocityVs128x64Percent'],5),e(q['cisFlowVs128x64Percent'])])

new('1. 三次映射复测：设置与主结果','六页对照版')
p('1D2V 线性 BGK，x 为壁面法向，间距 L=1，P1-DG。Couette 壁速 [-0.01,0.01]、无体力；Poiseuille 静止壁、a_y=0.01。等温完全漫反射，H 初值为零。主网格 Nx64，相对宏观更新 1e-5 停机。')
table(['设置','本次复测'],[
    ['新速度网格','v_x、v_y 各 32 点：ξ_i=(2i-1-N)/(N-1)，v_i=4ξ_i³；裸权重 Δv_i=24ξ_i²/(N-1)，二维权重 w_ij=Δv_x,i Δv_y,j M0，M0=exp(-v_x²-v_y²)/π。两方向均截断 [-4,4]。'],
    ['零附近分辨率','仍为 1024 点；min|v|≈0.000134269，|v|<0.1 内每方向 10 点。旧网格为法向半区 Legendre 截断 8、流向 Hermite。'],
    ['保持一致的设置','宏观稳定参数 [1,δ,δ,δ]；ω=min(1,2δΔx/√π)。只混合下一轮碰撞源，不直接改变 H；CIS 监测 kinetic 矩，GSIS 监测混合源矩。'],
],[1.05,5.8],size=8.7)
table(['算例 / δ','CIS / GSIS 步数','新绝对 L2 差','速度相对差 [%]：旧 → 新','应力相对差 [%]：旧 → 新'],[
    [flowlabel(r['flow'],float(r['delta'])),f"{r['cisIterations']} / {r['gsisIterations']}",e(r['absoluteVelocityL2']),
     pair(old[(r['flow'],float(r['delta']))]['relativeVelocityL2Percent'],r['relativeVelocityL2Percent']),
     pair(old[(r['flow'],float(r['delta']))]['relativeShearL2Percent'],r['relativeShearL2Percent'])] for r in main
],[1,1.2,1.2,2,2],size=8.7)
p('<b>结果：</b>C=Couette，P=Poiseuille。六组迭代次数与旧网格完全相同；δ=1 的 GSIS/CIS 速度差减小，δ=0.1、10 略增。这里比较的是算法间一致性，不能把差异增减直接解释为相对精确解的精度增减。',small=True)
p('δ=10 的迭代数比仍为 Couette 1.74、Poiseuille 3.80。第 6 页严格停机检查发现原方程残差平台；以保持 CIS 离散解为目标，当前耦合一致性检查未通过，不能据迭代数减少认定同一离散解的加速。',small=True)
table(['δ','Couette：CIS / GSIS 原方程 R','Couette：GSIS 耦合 R','Poiseuille：CIS / GSIS 原方程 R','Poiseuille：GSIS 耦合 R'],[
    [f(delta),e(summary[('couette',delta)]['cisOriginalResidual'])+' / '+e(summary[('couette',delta)]['gsisOriginalResidual']),
     e(summary[('couette',delta)]['gsisCoupledResidual']),
     e(summary[('poiseuille',delta)]['cisOriginalResidual'])+' / '+e(summary[('poiseuille',delta)]['gsisOriginalResidual']),
     e(summary[('poiseuille',delta)]['gsisCoupledResidual'])] for delta in [.1,1,10]
],[.35,2.1,1.2,2.1,1.2],size=8)

new('2. 速度剖面：旧网格与新网格','剖面与解读')
picture(profiles,maxheight=349)
p('<b>坐标与曲线：</b>横轴为法向位置 x，纵轴为单元平均流向速度 u_y。由上到下为 δ=0.1、1、10，左 Couette、右 Poiseuille。灰点线为旧网格 CIS，蓝线为三次映射 CIS，红虚线为三次映射 GSIS。',small=True)
p('<b>观察：</b>δ=0.1 的剖面变化最明显，δ=1、10 的新旧曲线接近；Couette 反对称、Poiseuille 对称。新网格下 CIS/GSIS 曲线接近，仍需用第 1 页的小差异指标作定量区分。',small=True)
p('<b>近壁说明：</b>所画点为壁内单元平均速度；气体壁面迹与规定壁速之间的滑移另外记录。法向速度、密度和温度扰动接近零，详见配套诊断 CSV。',small=True)

new('3. δ=10：四类迭代指标一起检查','更新与残差')
picture(histories,maxheight=344)
p('<b>从左到右：</b>绝对宏观量 RMS 更新、实际相对更新、绝对方程残差、相对方程残差。上排 Couette，下排 Poiseuille；横轴均为外层迭代数，纵轴用对数坐标。',small=True)
p('<b>曲线：</b>更新图的蓝/红为 CIS/GSIS kinetic，黄色为 GSIS 混合源；相对更新监测 CIS kinetic 与 GSIS 混合源，虚线 1e-5 为停机阈值。残差图蓝色为 CIS 原方程，红色为 GSIS 原 CIS 方程，黄色为 GSIS 耦合方程；相对残差统一除以同算例 CIS 首轮残差。',small=True)
p('<b>结论：</b>GSIS 的更新与残差均显著下降，但最终原 CIS 方程残差仍高于 CIS。δ=10 新网格 GSIS 的原方程残差比旧网格略高；宏观混合反馈造成的解差异仍需单独处理。全部 δ 的四类曲线见配套图 05-08。',small=True)
p('指标定义：绝对解差 Eabs=√(Δx∑(uG-uC)²)，相对差=100‖uG-uC‖₂/‖uC‖₂，应力同理。绝对更新为四个节点宏观量 RMS 变化的最大值；绝对方程残差采用 Maxwell 权重和 DG 质量矩阵逆范数。',small=True)

new('4. 六个 δ 的 Poiseuille 流量对照','Q 与 Qp')
p('归一化 Qp=Q/(a_yL²)=100Q；实际原始流量为 Q=Qp/100。下表均为新网格、Nx64、相对更新阈值 1e-5；偏差=100(Qp/Qp,paper-1)，保留正负号。',small=True)
table(['δ','文献 BGK Qp','CIS Qp','GSIS Qp','CIS 偏差 [%]','GSIS 偏差 [%]','CIS / GSIS 步数'],[
    [f(r['delta']),r['paperBGKQp'],f(r['cisQp'],9),f(r['gsisQp'],9),f(r['cisVsPaperPercent'],5),f(r['gsisVsPaperPercent'],5),
     f"{r['cisIterations']} / {r['gsisIterations']}"] for r in flowrows
],[.4,1.2,1.3,1.3,1.15,1.15,1.2],size=8.7)
picture(flowfig,maxheight=175)
p('<b>图：</b>左为 Qp 随 δ 的变化及文献表值；右为有符号相对文献偏差，蓝 CIS、红 GSIS，淡点线为旧速度网格。δ=0.1 的改善最大；δ=5、10 的 GSIS 仍约低 0.29%、0.31%。',small=True)
p('文献为 Sharipov (1999), Table 1, printed p.398，a/b=0 的 BGK work [36] 列。δ=Lp/(μc0)，c0=√(2kT/m)，Kn=√π/(2δ)，L 为完整间距、本次为 1。另存的 CIS 1e-9 补查最大文献偏差约 0.0726%，没有混入本表。',small=True)

new('5. 空间与速度网格检查','Nx 与速度分辨率')
gridrows=[]
for flow in ['couette','poiseuille']:
    for delta in [.1,1,10]:
        group=sorted([r for r in spatial if r['flow']==flow and float(r['delta'])==delta],key=lambda r:float(r['cells']))
        og=sorted([r for r in oldspatial if r['flow']==flow and float(r['delta'])==delta],key=lambda r:float(r['cells']))
        assert [r['cisIterations'] for r in group]==[r['cisIterations'] for r in og]
        assert [r['gsisIterations'] for r in group]==[r['gsisIterations'] for r in og]
        gridrows.append([flowlabel(flow,delta),' / '.join(r['cisIterations'] for r in group),
                         ' / '.join(r['gsisIterations'] for r in group),' / '.join(f(r['relativeVelocityL2Percent'],4) for r in group)])
table(['算例 / δ','CIS：Nx32 / 64 / 128','GSIS：Nx32 / 64 / 128','新网格速度相对差 [%]：Nx32 / 64 / 128'],gridrows,[1,1.5,1.5,3.2],size=8.7)
picture(spacefig,maxheight=118)
p('上图左为 δ=10 的新旧网格速度差，右为 GSIS 次数。18 组空间对照的次数均与旧网格相同；δ=10 的 GSIS 次数随 Nx 增加，未体现网格无关加速。ω 与 Δx 成正比，Nx64 时三组 δ 对应 0.001763/0.01763/0.1763，因此不能仅据三点趋势判定空间收敛阶。',small=True)
table(['δ','CIS 速度加密变化：C [%]','CIS 速度加密变化：P [%]','CIS 流量加密变化：P [%]'],vrows,[.5,1.8,1.8,1.8],size=8.2)
p('速度加密表比较当前 32×32 与同三次映射 128×64，固定 Nx64、相同停机阈值；中间还测了 64×32、128×32。24 组速度配置中的 CIS/GSIS 次数均保持原值。这是有限网格差异，不是精确解误差。',small=True)
maxinv=max(r['collisionInvariantDefect'] for r in diagnostics);maxflux=max(r['wallMassFlux'] for r in diagnostics)
maxtrace=max(r['macroBoundaryTraceError'] for r in diagnostics if r['method']=='gsis')
p(f'空间对照均达到相对更新停机条件；碰撞不变量缺口最大 {e(maxinv)}，壁面质量通量最大 {e(maxflux)}，宏观边界 trace 施加误差最大 {e(maxtrace)}。这些检查通过仍不足以证明宏观/kinetic 离散闭合一致。',small=True)

new('6. 收紧停机后：大小 δ 均存在残差平台','代码一致性诊断')
p('固定 Nx64、32×32 三次映射、相同边界和稳定参数；阈值 1e-5/1e-7/1e-9 各自从 H=0 开始。停机监测的是单元完整密度、温度和速度模的相对更新；下表为 Rrel=Rabs/Rabs,CIS(1)，同算例、不同阈值使用相同分母，数值不是百分数。',small=True)
table(['Poiseuille δ=10：阈值','CIS 原方程 Rrel','GSIS 原 CIS 方程 Rrel','GSIS 耦合方程 Rrel','CIS / GSIS 步数'],[
    [f(t),e(stop[('poiseuille',10.,t)]['cisRelativeOriginalResidual']),
     e(stop[('poiseuille',10.,t)]['gsisRelativeOriginalResidual']),e(stop[('poiseuille',10.,t)]['gsisRelativeCoupledResidual']),
     stop[('poiseuille',10.,t)]['cisIterations']+' / '+stop[('poiseuille',10.,t)]['gsisIterations']] for t in [1e-5,1e-7,1e-9]
],[1.35,1.25,1.65,1.5,1.05],size=8.5)
p('<b>收紧阈值有作用：</b>CIS 原方程、GSIS 耦合方程残差都下降约四个数量级；GSIS 对原 CIS 方程的残差却仍约为 4.35e-3。问题不能归结为停机阈值不够小。',small=True)
table(['算例 / δ','GSIS 原方程 Rrel：1e-5 → 1e-9','CIS 原方程 Rrel @1e-9','GSIS 耦合 Rrel @1e-9','混合源缺口范数 @1e-9'],[
    [flowlabel(flow,delta),e(stop[(flow,delta,1e-5)]['gsisRelativeOriginalResidual'])+' → '+e(stop[(flow,delta,1e-9)]['gsisRelativeOriginalResidual']),
     e(stop[(flow,delta,1e-9)]['cisRelativeOriginalResidual']),e(stop[(flow,delta,1e-9)]['gsisRelativeCoupledResidual']),
     e(stop[(flow,delta,1e-9)]['relativeSourceDefect'])] for delta in [.1,1,10] for flow in ['couette','poiseuille']
],[1,2.1,1.4,1.4,1.45],size=8.1)
p('<b>小 δ 也存在：</b>δ=0.1 的原方程残差平台约为 Couette 1.92e-6、Poiseuille 3.43e-6；量级很小，但继续收紧迭代更新并不能使其归零，个别组还会略增。',small=True)
p('<b>已证实：</b>严格停机后，混合源宏观量与 H 自身矩仍有非零缺口；该源项缺口经与方程残差相同的范数及分母归一化后，与原方程残差平台相符。以保持原 CIS 离散解为目标，当前 GSIS 耦合一致性检查未通过。',small=True)
p('<b>尚待定位：</b>需要继续核查宏观/kinetic 离散通量、边界闭合及宏观反馈投影等环节。现有测试尚不能确定具体哪一处实现有误。仅满足宏观更新停机条件，不足以证明已解到同一个原 CIS 离散稳态。',small=True)
p('完整 18 行停机对照及诊断数据：output/channel_cubic_figures_20261109/data/；复现停机测试：addpath(\'debug\'); check_cubic_stopping_residuals。',small=True)

assert page==6
C.save();(TMP/'page_index.json').write_text(json.dumps(toc,ensure_ascii=False,indent=2)+'\n')
extra_titles={
    '16_cubic_profiles_with_old_grid':'新旧速度网格的六组速度剖面',
    '17_cubic_all_four_convergence_metrics_delta10':'δ=10 的绝对/相对更新与绝对/相对方程残差',
    '18_cubic_all_six_flowrates_old_vs_new':'六个 δ 的流量与文献偏差：新旧网格对照',
    '19_cubic_spatial_refinement_old_vs_new':'空间加密：新旧网格的解差异和迭代次数',
    '20_cubic_current_velocity_refinement':'当前三次映射速度网格的加密检查',
    '21_cubic_stopping_sensitivity_all_six_cases':'六组算例的停机阈值与相对方程残差平台',
}
index=(FIG/'index.html').read_text()
if '<h2 id="cubic-retest">' in index:index=index.split('<h2 id="cubic-retest">')[0]+'</main></html>'
sections='<h2 id="cubic-retest">六页复测报告的组合图片</h2>'
for name,title in extra_titles.items():
    sections+=f'<section><h2>{title}</h2><p><a href="{name}.png">PNG</a> · <a href="{name}.svg">SVG</a></p><img loading="lazy" src="{name}.png" alt="{title}"></section>'
index=index.replace('</main></html>',sections+'</main></html>');(FIG/'index.html').write_text(index)
prov=json.loads((FIG/'provenance.json').read_text())
prov['figures']=[r for r in prov['figures'] if int(r['name'][:2])<16]
prov['figures'] += [dict(name=name,title=title,caption='Relative macro stops 1e-5, 1e-7, 1e-9; common CIS first residual.' if name.startswith('21_') else 'Cubic-grid retest, reference macro stop 1e-5.') for name,title in extra_titles.items()]
(FIG/'provenance.json').write_text(json.dumps(prov,ensure_ascii=False,indent=2)+'\n')
with zipfile.ZipFile(FIG.with_suffix('.zip'),'w',zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(FIG.rglob('*')):
        if path.is_file():archive.write(path,path.relative_to(FIG.parent))
print('PDF:',TARGET);print('Pages:',page);print('All retest metrics:',DATA/'old_vs_cubic_main_metrics.csv')
