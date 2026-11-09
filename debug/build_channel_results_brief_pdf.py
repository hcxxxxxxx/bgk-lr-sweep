"""Six-page presentation summary, using the verified report's layout helpers."""
from pathlib import Path

# Load only definitions and saved data; exclude diagnostic writes/report pages.
source = Path(__file__).with_name('build_channel_results_pdf.py').read_text()
source = source.split('# Supplementary diagnostics, recomputed only from saved fields.')[0]
source = source.replace("TMP=ROOT/'tmp/pdfs/channel_results_20261109'",
                        "TMP=ROOT/'tmp/pdfs/channel_results_brief_20261109'")
source = source.replace("TARGET=OUT/'1D2V_CIS_GSIS_指标与逐图说明_20261109.pdf'",
                        "TARGET=OUT/'1D2V_CIS_GSIS_六页汇报版_无独立BGK_20261109.pdf'")
source = source.replace('全部指标、算例结果与逐图说明', '六页汇报版')
exec(compile(source, str(Path(__file__)), 'exec'))

import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10,
                     'axes.labelsize': 9, 'legend.fontsize': 8})
fig, axes = plt.subplots(3, 2, figsize=(11.8, 5.5), constrained_layout=True)
for row, delta in enumerate([.1, 1, 10]):
    for col, flow in enumerate(['couette', 'poiseuille']):
        case = lookup[(flow, delta)]
        ax = axes[row, col]
        x = case['cis']['mesh']['center']
        ax.plot(x, case['cis']['cellMoments'][:, 2], color='#2568A6', label='CIS')
        ax.plot(x, case['gsis']['cellMoments'][:, 2], '--', color='#D34A32', label='GSIS')
        ax.set(title=f'{FLOW[flow]} | delta={delta:g}', xlabel='$x$', ylabel='$u_y$')
        ax.grid(alpha=.15)
        if row == col == 0:
            ax.legend(ncol=2, loc='upper left')
profile_figure = TMP/'six_velocity_profiles.png'
fig.savefig(profile_figure, dpi=230, facecolor='white')
plt.close(fig)

new('1. 核心结果与统一设置', '六页汇报版')
p('测试：1D2V 线性 BGK，x 为壁面法向、y 为流向，间距 L=1。主网格 Nx=64、P1-DG，速度 32×32；完全漫反射等温壁面，H 初值为零。Couette 壁速 [-0.01,0.01]、无体力；Poiseuille 两壁静止、a_y=0.01。')
p('CIS 和 GSIS 均以相对宏观更新小于 1e-5 停机。GSIS 只将混合后的宏观量回传到下一轮碰撞源，不直接改 H；宏观稳定参数为 [1,δ,δ,δ]，ω=min(1,2δΔx/√π)。', small=True)
table(['速度离散','节点与权重（当前全部主算例相同）'],[
    ['法向 v_x：32 点','[-8,0]、[0,8] 各 16 点 Gauss-Legendre。ξ_i∈(-1,1)，v_x=±4(1+ξ_i)；半区权重为 4ω_i exp(-v_x²)。min|v_x|≈0.042396，无零节点。'],
    ['流向 v_y：32 点','32 点 Gauss-Hermite，定义于全实轴，权函数 exp(-v_y²)；节点范围约 [-7.12581,7.12581]，不另设截断。'],
    ['二维组合：1024 点','张量积权重 w_ij=w_x,i w_y,j/π，已包含 M0=exp(-v_x²-v_y²)/π。当前未额外做零附近的坐标映射加密。'],
],[1.2,5.7],size=9)
table(['流动','δ','CIS / GSIS 步数','迭代数比','速度绝对 L2 差','速度相对差 [%]','应力相对差 [%]'], [
    [FLOW[r['flow']],f(r['delta']),f"{r['cisIterations']} / {r['gsisIterations']}",
     f(r['iterationRatio'],4),e(r['absoluteVelocityL2']),f(r['relativeVelocityL2Percent'],5),
     f(r['relativeShearL2Percent'],5)] for r in main
], [1.3,.45,1.25,.8,1.35,1.4,1.4])
sub('可以用于汇报的结论')
p('<b>δ=0.1、1：</b>两种方法的解接近，迭代次数相同。<b>δ=10：</b>Couette 从 82 步降至 47 步，Poiseuille 从 251 步降至 66 步；同时，GSIS 与 CIS 的速度差异分别约 0.0774% 和 0.300%。')
p('迭代数比=N_CIS/N_GSIS，并非运行时间加速比。表中的解差异均以同网格 CIS 为基准，反映两算法的一致性。')

new('2. 速度剖面：六组算例合并看', '图与解读')
picture(profile_figure, maxheight=349)
p('<b>看什么：</b>横轴为法向位置 x，纵轴为单元平均流向速度 u_y；从上到下为 δ=0.1、1、10。蓝实线为 CIS，红虚线为 GSIS；左列 Couette，右列 Poiseuille。',small=True)
p('<b>说明什么：</b>Couette 关于中心反对称，Poiseuille 对称；气体侧近壁速度体现滑移。两种方法的曲线在当前显示尺度下接近，细微差异用第 1 页的绝对及相对 L2 指标衡量。',small=True)
p('<b>定量比较：</b>Couette 的 GSIS-CIS 速度相对差依次为 6.28e-7%、0.001515%、0.07736%；Poiseuille 依次为 2.08e-5%、0.009321%、0.3003%。',small=True)

new('3. 迭代下降：更新量与方程残差', '图与解读')
picture(FIG/'pdf_panels'/'06_relative_updates_10.png',maxheight=162)
picture(FIG/'pdf_panels'/'07_absolute_residuals_10.png',maxheight=162)
p('<b>上排：δ=10 的相对更新。</b>横轴为迭代数；蓝线监测 CIS kinetic 矩，红线监测 GSIS 混合源矩。水平虚线 1e-5 是实际停机阈值。左 Couette、右 Poiseuille，GSIS 更早达到阈值。',small=True)
p('<b>下排：同两组的绝对方程残差。</b>蓝线为 CIS 原方程，红线为 GSIS 代入原 CIS 方程，黄虚线为 GSIS 代入混合源耦合方程；红黄使用相同 H、不同碰撞源。红线确实下降，但最终高于 CIS。',small=True)
p('<b>结论：</b>“更新足够小”“方程残差下降”“与 CIS 解相同”是三个不同判断。δ=0.1、1 没有迭代数收益，其计数见第 1 页。',small=True)

new('4. Poiseuille 流量与文献对照', '流量 Qp')
p('<b>为什么乘 100：</b>Q=∫u_y dx 是实际计算的原始流量；Qp=Q/(a_yL²) 消除驱动力幅值。当前 a_y=0.01、L=1，所以 Qp=100Q。这是归一化，不是把实际流量改大。')
table(['δ','方法','原始 Q','归一化 Qp','文献 BGK Qp','相对文献差 [%]'], [
    [f(r['delta']),r['method'],f(float(r['Qp'])/100,9),f(r['Qp'],9),r['paperBGK'],
     f(r['signedRelativeDifferencePercent'],5)] for r in current
], [.5,.8,1.5,1.5,1.3,1.6])
picture(FIG/'pdf_panels'/'11_literature_Qp_differences.png',maxheight=166)
p('<b>图：</b>横轴为 δ；左为有符号差 Qp-Qp,paper，右为 100(Qp/Qp,paper-1)。金色当前 CIS、红色当前 GSIS、蓝色历史精细 CIS；灰带仅表示表格四位小数的舍入区间。',small=True)
p('<b>来源与口径：</b>Sharipov (1999), Table 1, p.398，a/b=0 的 BGK work [36] 列。历史精细 CIS 为 Nx256、Nv512×8、方程阈值 1e-10，最大文献差约 0.0719%；与当前设置分别展示。采用黏性平均自由程约定时 Kn=√π/(2δ)，本次间距为 1、δ 定义已对齐。',small=True)

new('5. 网格检查：次数是否稳定，差异是否下降', '空间与速度分辨率')
gridrows=[]
for flow in ['couette','poiseuille']:
    for delta in [.1,1,10]:
        grp=sorted([r for r in spatial if r['flow']==flow and float(r['delta'])==delta],key=lambda r:float(r['cells']))
        gridrows.append([label(flow,delta),' / '.join(r['cisIterations'] for r in grp),
                         ' / '.join(r['gsisIterations'] for r in grp),
                         ' / '.join(f(r['relativeVelocityL2Percent'],4) for r in grp)])
table(['算例 / δ','CIS：Nx32 / 64 / 128','GSIS：Nx32 / 64 / 128','速度相对差 [%]：Nx32 / 64 / 128'],
      gridrows,[1.6,1.5,1.5,3])
picture(FIG/'pdf_panels'/'10_spatial_differences_and_mixing_part2.png',maxheight=177)
p('<b>图：</b>左为 δ=10 时 GSIS-CIS 的速度相对 L2 差，蓝 Couette、红 Poiseuille；右为随网格改变的 ω。CIS 次数不变，δ=10 的 GSIS 次数随加密增加，因此当前没有证明网格无关加速。',small=True)
p('<b>需要一起看：</b>加密时差异减小，但 ω 也变小，不能仅由该趋势断言二阶精度。历史精细流量检查中，联合空间/速度加密的最大变化为 0.000241%；这不能替代当前 32×32 网格的精度研究。',small=True)

new('6. 指标速查与最终残差', '定义与数据来源')
table(['指标','定义 / 读法'],[
    ['解的绝对 / 相对差','绝对 L2=√(Δx∑(uG-uC)²)；相对 L2=100‖uG-uC‖₂/‖uC‖₂。应力 Pxy 用同一形式。'],
    ['相对更新 D_ref','单元积分的完整密度、温度、速度模，相邻轮相对变化的最大值；CIS 监测 kinetic 矩，GSIS 监测混合源矩。'],
    ['绝对残差 R','离散 kinetic-DG 稳态方程缺口的速度加权、质量矩阵逆范数，包含体项与面项；不是相对参考解误差。'],
    ['相对残差','R(k)/R_CIS(1)，同算例各曲线共用 CIS 首轮分母；它不是实际停机指标。'],
    ['流量与近零量','Qp=Q/(a_yL²)；文献差=100(Qp/Qp,paper-1)。Couette 净 Q≈0，不报告相对流量差。'],
],[1.1,5.6],size=9)
sub('最终绝对残差与 GSIS 停机值（Nx64）')
table(['算例 / δ','CIS 原方程 R','GSIS 原方程 R','GSIS 耦合方程 R','GSIS 相对更新'],[
    [label(r['flow'],float(r['delta'])),e(r['cisOriginalResidual']),e(r['gsisOriginalResidual']),
     e(r['gsisCoupledResidual']),e(r['gsisRelativeMacroUpdate'])] for r in main
],[1.7,1.3,1.3,1.3,1.3])
p('所有六组均达到相对更新 1e-5。δ=10 的原方程残差与 CIS 仍有差距；宏观线性系统解到舍入误差附近，也不能证明整个 GSIS 保持 CIS 的同一个离散稳态。',small=True)
p('数据：results/channel_presentation_20261109/presentation_data.mat；表格与原图：output/channel_figures_20261109/data/ 及上级目录。详细版保留逐图解释、壁面 trace、守恒及参考程序诊断。本版仅重新整理已保存结果，未改动求解器。',small=True)

assert page == 6
C.save()
(TMP/'page_index.json').write_text(json.dumps(toc,ensure_ascii=False,indent=2)+'\n')
print(TARGET)
print('Pages:',page)
