"""Chinese PDF report with metric definitions, all current results and figure explanations.

Dependencies can be installed in /private/tmp/gsis_pdf_deps. This script uses
only saved results; it never changes the kinetic or macro solvers.
"""
import sys
sys.path.insert(0, '/private/tmp/gsis_pdf_deps')
from pathlib import Path
import csv
import html
import io
import json
import math

import numpy as np
from scipy.io import loadmat
from PIL import Image as PILImage
import matplotlib
matplotlib.use('Agg')
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.colors import HexColor, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, Table, TableStyle, Image
from reportlab.lib.styles import ParagraphStyle

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'output/channel_figures_20261109'
DATA=FIG/'data'
TMP=ROOT/'tmp/pdfs/channel_results_20261109'
OUT=ROOT/'output/pdf'
TMP.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
TARGET=OUT/'1D2V_CIS_GSIS_指标与逐图说明_20261109.pdf'
FLOW={'couette':'Couette','poiseuille':'Poiseuille'}


def csv_rows(path):
    with path.open() as f:return list(csv.DictReader(f))


main=csv_rows(DATA/'main_results.csv')
spatial=csv_rows(DATA/'spatial_iterations_and_differences.csv')
historical=csv_rows(DATA/'historical_literature_comparison.csv')
current=csv_rows(DATA/'current_literature_comparison.csv')
vcheck=csv_rows(DATA/'historical_velocity_check_summary.csv')
viter=csv_rows(DATA/'historical_velocity_iterations.csv')
saved=loadmat(ROOT/'results/channel_presentation_20261109/presentation_data.mat',simplify_cells=True)
cases=sorted([c for c in saved['cases'] if c['cells']==64],key=lambda c:(c['delta'],c['flow']!='couette'))
lookup={(c['flow'],float(c['delta'])):c for c in cases}
summary={(r['flow'],float(r['delta'])):r for r in main}
assert len(cases)==6 and len(main)==6 and len(spatial)==18


def f(v,digits=5):
    v=float(v)
    return '不适用' if math.isnan(v) else f'{v:.{digits}g}'


def e(v):return f'{float(v):.3e}'
def tag(delta):return f'{delta:g}'.replace('.','p')
def label(flow,delta):return f'{FLOW[flow]} / δ={delta:g}'
def esc(s):return html.escape(str(s))


pdfmetrics.registerFont(TTFont('CN','/System/Library/Fonts/STHeiti Light.ttc',subfontIndex=1))
pdfmetrics.registerFont(TTFont('CNB','/System/Library/Fonts/STHeiti Medium.ttc',subfontIndex=1))
pdfmetrics.registerFontFamily('CN',normal='CN',bold='CNB',italic='CN',boldItalic='CNB')
pdfmetrics.registerFont(TTFont('Latin','/System/Library/Fonts/Supplemental/Arial.ttf'))
W,H=landscape(A4);MARGIN=38;CW=W-2*MARGIN
INK=HexColor('#203447');BLUE=HexColor('#2568A6');MUTED=HexColor('#536878')
PALE=HexColor('#EEF3F8');RULE=HexColor('#D6E0E8')
BODY=ParagraphStyle('body',fontName='CN',fontSize=10.5,leading=16,textColor=INK,wordWrap='CJK',spaceAfter=7)
SMALL=ParagraphStyle('small',parent=BODY,fontSize=9,leading=13,spaceAfter=5)
CELL=ParagraphStyle('cell',parent=BODY,fontSize=9,leading=12,spaceAfter=0)
CELLHEAD=ParagraphStyle('cellhead',parent=CELL,fontName='CNB',textColor=BLUE)
C=canvas.Canvas(str(TARGET),pagesize=(W,H),pageCompression=1)
C.setTitle('1D2V CIS / GSIS：全部指标、算例结果与逐图说明')
C.setAuthor('Numerical results report')
C.setSubject('Couette / Poiseuille, delta=0.1/1/10, common stopping, figure interpretation')
page=0;y=0;toc=[];formula_id=0


def new(title,section='结果与指标'):
    global page,y
    if page:C.showPage()
    page+=1;y=H-72
    C.setFillColor(BLUE);C.rect(0,H-7,W,7,fill=1,stroke=0)
    C.setFont('CN',9);C.setFillColor(MUTED)
    C.drawString(MARGIN,H-25,'1D2V BGK  |  CIS / GSIS')
    C.drawRightString(W-MARGIN,H-25,'2026-11-09  /  '+section)
    C.setFont('CNB',19);C.setFillColor(INK);C.drawString(MARGIN,H-53,title)
    C.setStrokeColor(RULE);C.line(MARGIN,32,W-MARGIN,32)
    C.setFont('CN',8);C.setFillColor(MUTED)
    C.drawString(MARGIN,19,'所有速度与流量均为无量纲量；当前主算例的停机标准为相对宏观更新 < 1e-5。')
    C.drawRightString(W-MARGIN,19,f'{page:02d}')
    bookmark=f'p{page}';C.bookmarkPage(bookmark);C.addOutlineEntry(title,bookmark,level=0)
    toc.append({'page':page,'title':title,'section':section})


def p(text,small=False,gap=5):
    global y
    paragraph=Paragraph(text,SMALL if small else BODY)
    _,height=paragraph.wrap(CW,1000)
    if y-height<43:raise RuntimeError(f'Page {page}: paragraph overflow at {y-height:.1f}: {text[:60]}')
    paragraph.drawOn(C,MARGIN,y-height);y-=height+gap


def sub(text):
    global y
    y-=4;C.setFont('CNB',12);C.setFillColor(BLUE);C.drawString(MARGIN,y-13,text);y-=24


def table(headers,rows,widths=None,size=9):
    global y
    if widths is None:widths=[CW/len(headers)]*len(headers)
    else:widths=[v/sum(widths)*CW for v in widths]
    style=ParagraphStyle('tablecell',parent=CELL,fontSize=size,leading=size+3)
    content=[[Paragraph(esc(v),CELLHEAD) for v in headers]]
    content.extend([[Paragraph(esc(v),style) for v in row] for row in rows])
    t=Table(content,colWidths=widths,hAlign='LEFT')
    t.setStyle(TableStyle([
        ('BACKGROUND',(0,0),(-1,0),PALE),('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#FAFCFE')]),
        ('LINEBELOW',(0,0),(-1,0),.8,BLUE),('LINEBELOW',(0,1),(-1,-1),.3,RULE),
        ('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),5),
        ('LEFTPADDING',(0,0),(-1,-1),7),('RIGHTPADDING',(0,0),(-1,-1),7),('VALIGN',(0,0),(-1,-1),'MIDDLE'),
    ]))
    _,height=t.wrap(CW,1000)
    if y-height<43:raise RuntimeError(f'Page {page}: table overflow at {y-height:.1f}')
    t.drawOn(C,MARGIN,y-height);y-=height+10


def picture(path,maxheight=255,width=None):
    global y
    width=width or CW
    with PILImage.open(path) as im:iw,ih=im.size
    scale=min(width/iw,maxheight/ih);dw,dh=iw*scale,ih*scale
    if y-dh<43:raise RuntimeError(f'Page {page}: image overflow')
    Image(str(path),width=dw,height=dh).drawOn(C,MARGIN+(CW-dw)/2,y-dh)
    y-=dh+7


def formula(expression):
    global y,formula_id
    formula_id+=1;path=TMP/f'formula_{formula_id:02d}.png'
    math_to_image(expression,str(path),prop=FontProperties(size=13),dpi=210,color='#203447')
    with PILImage.open(path) as im:iw,ih=im.size
    dw,dh=iw*72/210,ih*72/210
    if dw>CW-24:dh*=((CW-24)/dw);dw=CW-24
    C.setFillColor(PALE);C.roundRect(MARGIN,y-dh-15,CW,dh+16,5,fill=1,stroke=0)
    Image(str(path),width=dw,height=dh).drawOn(C,MARGIN+12,y-dh-7)
    y-=dh+24


def key(c):return (c['flow'],float(c['delta']))


# Supplementary diagnostics, recomputed only from saved fields.
diagnostics=[];wall=[]
for case in cases:
    c,g=case['cis'],case['gsis'];v=g['velocity'];w=v['weight'];vx=v['vx'];vy=v['vy']
    hw=g['H'][[0,-1],:]
    rho=hw@w;temp=hw@(w*(vx*vx+vy*vy-1));uy=hw@(w*vy)
    target=np.column_stack((rho+temp,np.zeros(2),uy,temp))
    hdg=g['synthetic']['hdg'];hat=hdg['Uhat'][[0,-1],:]
    applied=float(np.max(abs(target-hat)))
    extracted=float(np.max(abs(target-g['synthetic']['hot']['boundary']['trace'])))
    assert max(applied,extracted)<1e-12
    kinwall=g['kineticTrace'][[0,-1],:]
    kin_shear=kinwall@(w*vx*vy)
    macro_flux=np.vstack((hdg['outwardFlux'][:4,0],hdg['outwardFlux'][4:,-1]))
    macro_shear=np.array([-1,1])*macro_flux[:,2]/2
    wall.append(dict(flow=case['flow'],delta=float(case['delta']),
                     macroBoundaryDataError=extracted,macroAppliedTraceError=applied,
                     leftGasUy=float(uy[0]),rightGasUy=float(uy[1]),
                     leftSlip=float(g['slipVelocity'][0]),rightSlip=float(g['slipVelocity'][1]),
                     leftKineticPxy=float(kin_shear[0]),rightKineticPxy=float(kin_shear[1]),
                     leftMacroPxy=float(macro_shear[0]),rightMacroPxy=float(macro_shear[1])))
    for method,r in [('CIS',c),('GSIS',g)]:
        d=r['diagnostics'];qm=r['cellMoments'];heat=r['heatFlux'];stress=r['shearStress']
        row=dict(flow=case['flow'],delta=float(case['delta']),method=method,iterations=int(r['iterations']),
                 relativeMacroUpdate=float(r['convergence']['referenceMacroUpdate']),
                 kineticResidual=float(d['kinetic']),coupledResidual=float(d['coupled']),
                 wallResidual=float(d['wall']),wallMassFlux=float(d['wallMassFlux']),
                 meanDensityError=float(d['massError']),collisionInvariantDefect=float(d['collisionInvariantDefect']),
                 sourceMomentGap=float(d['sourceMomentGap']),
                 maximumRhoPerturbation=float(np.max(abs(qm[:,0]))),maximumUx=float(np.max(abs(qm[:,1]))),
                 maximumTPerturbation=float(np.max(abs(qm[:,3]))),
                 maximumAbsHeatX=float(np.max(abs(heat[:,0]))),maximumAbsHeatY=float(np.max(abs(heat[:,1]))),
                 meanPxy=float(np.mean(stress)),maximumAbsPxy=float(np.max(abs(stress))),
                 relativePxySpan=(float(np.ptp(stress)/abs(np.mean(stress))) if case['flow']=='couette' else float('nan')),
                 massBalanceDefect=float(d['balance'][0]),streamwiseMomentumBalanceDefect=float(d['balance'][1]),
                 normalMomentumBalanceDefect=float(d['balance'][2]),energyBalanceDefect=float(d['balance'][3]),
                 recordedLoopSeconds=float(r['wallTime']))
        if method=='GSIS':
            row.update(hdgLocalResidual=float(hdg['localResidual']),hdgTraceResidual=float(hdg['traceResidual']),
                       hdgFluxContinuityResidual=float(hdg['fluxContinuityResidual']),
                       hdgConservationResidual=float(hdg['conservationResidual']))
        else:row.update(hdgLocalResidual=float('nan'),hdgTraceResidual=float('nan'),
                        hdgFluxContinuityResidual=float('nan'),hdgConservationResidual=float('nan'))
        diagnostics.append(row)
for name,rows in [('all_solver_diagnostics.csv',diagnostics),('all_wall_diagnostics.csv',wall)]:
    with (DATA/name).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


new('1D2V CIS / GSIS：指标、算例结果与逐图说明','报告总览')
p('覆盖 Couette 与 Poiseuille 两种流动，新增 δ=0.1，并汇总 δ=0.1、1、10 的主结果、空间加密检查、文献流量对照以及 15 组图片。所有结论均来自已保存的数值结果。')
sub('本次最重要的四个结果')
p('<b>δ=0.1：</b>64 单元下 Couette 为 CIS/GSIS=5/5 步，Poiseuille 为 8/8 步。速度相对 CIS 的差异分别为 6.276e-7% 与 2.084e-5%，未减少迭代次数。')
p('<b>δ=1：</b>Couette 为 11/11 步，Poiseuille 为 23/23 步。两算法接近，但仍没有迭代数收益。')
p('<b>δ=10：</b>Couette 为 82/47 步，Poiseuille 为 251/66 步；迭代数比为 1.74 与 3.80，速度相对 CIS 的差异为 0.07736% 与 0.3003%。')
p('<b>接近 CIS 不等于接近参考解：</b>δ=0.1 时，两算法虽然几乎一致，但当前 32×32 速度网格下，相对独立 BGK 数值参考的速度 L2 偏差约为 Couette 4.34%、Poiseuille 0.640%。')
sub('如何阅读')
table(['内容','用途'],[
    ['设置与指标定义','先弄清当前解、比较基准、停机指标和残差归一化。'],
    ['数值汇总','直接取用于 PPT 的步数、速度/剪应力差异、流量、残差与网格结果。'],
    ['逐图解读','逐一解释图 01-15 的目的、坐标、曲线和可支持的结论。'],
    ['附录诊断与来源','查看壁面、守恒、宏观线性系统、单次计时和参考 Fortran 程序的边界条件。'],
],[1,5])

new('1. 模型、坐标、边界与迭代设置','计算设置')
formula(r'$v_x\partial_x H=\delta[\mathcal{P}(Q)-H]+2a_yv_y,\quad f=M_0(1+H),\quad M_0=e^{-(v_x^2+v_y^2)}/\pi$')
p('只离散空间 x；x 是壁面法向，y 是流向。Q=[ρ′,u_x,u_y,T′]，密度和温度输出是相对平衡态 1 的扰动。实际碰撞增益使用离散 Gram 投影，守恒四个碰撞不变量。')
table(['设置','当前主算例'],[
    ['区域与空间离散','x∈[-0.5,0.5]，间距 L=1；主算例 64 单元，P1-DG；加密为 32/64/128 单元。'],
    ['速度求积','32×32=1024 点；法向正、负半区各 16 点 Gauss-Legendre，范围 [-8,8]；流向 32 点 Gauss-Hermite。'],
    ['Couette','左、右壁切向速度 [-0.01,0.01]；a_y=0。'],
    ['Poiseuille','两壁静止；a_y=0.01。'],
    ['kinetic 壁面','等温完全漫反射；左壁 v_x>0、右壁 v_x<0 为入射；H_in=n_w+2v_yU_w。n_w 从出射分布按零质量通量求得。'],
    ['宏观壁面 trace','由当前内部 H([1,end],:) 的全速度矩给出 [p,0,u_y,T′]，p=ρ′+T′；切向 trace 是气体速度，含滑移。'],
    ['初值与密度规范','H=0；保留原 kinetic 的平均密度扰动为零约束。'],
    ['宏观稳定参数','[1,δ,δ,δ]，新增 δ=0.1 时为 [1,0.1,0.1,0.1]；没有调整默认算法。'],
],[1.2,5])
formula(r'$Q_s=Q_{\rm kin}+\omega(Q_{\rm macro}-Q_{\rm kin}),\qquad\omega=\min(1,2\delta\Delta x/\sqrt{\pi})$')
p('宏观求解只改变下一轮 kinetic 的碰撞源，不直接提升或替换 H。64 单元时，δ=0.1/1/10 的 ω 分别为 0.001763/0.01763/0.1763。',small=True)

new('2. 解的差异、剪切应力与流量指标','指标定义')
p('以下主结果以同网格 CIS 的 kinetic 矩为比较基准，GSIS 也取最终 kinetic H 的矩。宏观求解器的原始输出和混合源并不替代最终物理剖面。')
formula(r'$E_{\rm abs}=\sqrt{\Delta x\sum_i(u_{G,i}-u_{C,i})^2},\quad E_{\rm rel}=100\frac{E_{\rm abs}}{\sqrt{\Delta x\sum_i u_{C,i}^2}},\quad E_{\infty}=\max_i|u_{G,i}-u_{C,i}|$')
p('<b>速度指标：</b>u_C、u_G 是各 DG 单元平均 u_y。E_abs 是绝对 L2 差；E_rel 是全局相对 L2 百分比；E∞ 是最大绝对差。它们反映两算法的离散解一致性，不能直接称为相对精确解的误差。')
formula(r'$d(x)=u_y^{G}(x)-u_y^{C}(x),\qquad d_{\rm plot}(x)=100\,|d(x)|/\max_x|u_y^{C}(x)|$')
p('剖面图的中间面板保留差异的正负号；右侧曲线用整条 CIS 剖面的最大幅值归一化，避免 Couette 中心速度过零时除零。右图标题给出的 E_rel 与这条局部曲线的定义不同。')
formula(r'$P_{xy}(x)=\int v_xv_yH M_0\,d\boldsymbol{v},\quad Q=\int_{-L/2}^{L/2}u_y\,dx,\quad Q_p=Q/(a_yL^2)$')
p('<b>剪切应力：</b>代码输出 P_xy，宏观辅助量 σ_xy=2P_xy；本报告统一比较 P_xy。剪切应力相对 L2 差按速度的同一形式计算。Couette 稳态 P_xy 应近似常数；Poiseuille 满足 dP_xy/dx=a_y。')
p('<b>流量：</b>Q 是原始无量纲截面积分；本模型参考密度为 1，因此它也是线性化质量流率的对应量。当前 Poiseuille 的 a_y=0.01、L=1，故 Q_p=100Q。流量相对差为 100|Q_G-Q_C|/|Q_C|。对反对称 Couette，净 Q 约为零，不计算这个相对百分比。')

new('3. 方程残差与迭代更新：必须区分','指标定义')
formula(r'$\mathcal{R}(H;Q)=A_{\rm DG}H-M[\delta\mathcal{P}(Q)+G]-b_w(H),\qquad G=2a_yv_y$')
p('A_DG 包含离散输运及碰撞损失；M 是 DG 质量矩阵；b_w 按当前 H 重算漫反射边界。残差包含体项、内部面和壁面项，表示当前 H 对离散稳态方程的缺口。')
formula(r'$R_{\rm abs}=\sqrt{\frac{1}{L}\sum_jw_j\mathcal{R}_j^TM^{-1}\mathcal{R}_j},\qquad R_{\rm rel}^{k}=R_{\rm abs}^{k}/R_{\rm abs,CIS}^{1}$')
p('速度权重 w_j 已含 Maxwell 分布。相对残差统一除以同算例 CIS 第一次迭代的绝对残差；它不是除以右端范数，也不是除以最终残差。')
table(['图中曲线','代入的 H','碰撞增益所用宏观量','含义'],[
    ['蓝线：CIS 原方程','CIS 当前 H','Q_kin','检查 CIS 解对原 kinetic-DG 方程的满足程度。'],
    ['红线：GSIS 原方程','GSIS 当前 H','同一 H 自身的 Q_kin','检查 GSIS 解对原 CIS 方程的满足程度。'],
    ['黄色虚线：GSIS 耦合方程','与红线完全相同的 H','混合后的 Q_s','检查当前混合源耦合系统。'],
],[1.5,1,1.25,3])
p('红线和黄线的差别来自碰撞源。残差显著下降说明方程缺口减小；但某条残差下降，不单独证明两算法具有同一个最终离散解。')
formula(r'$D_{\rm kin}^{k}=\max_{\ell=1,\ldots,4}\sqrt{\frac{1}{2N}\sum_{a=1}^{2N}(Q_{\ell,a}^{k}-Q_{\ell,a}^{k-1})^2}$')
p('绝对更新是相邻两轮四个节点宏观量 RMS 更新的最大值。GSIS 另记录 Q_s 的绝对更新。它与方程残差、相对 CIS 的解差异是三个不同指标。',small=True)

new('4. 真正的停机指标与收敛状态','指标定义')
formula(r'$Z_{\rho,i}=\int_i(1+\rho^\prime)dx,\quad Z_{T,i}=\int_i(1+T^\prime)dx,\quad Z_{U,i}=\sqrt{(\int_i u_xdx)^2+(\int_i u_ydx)^2}$')
formula(r'$D_{\rm ref}^{k}=\max_{a\in\{\rho,T,U\}}\frac{\|Z_a^k-Z_a^{k-1}\|_2}{\|Z_a^k\|_2}<10^{-5}$')
p('CIS 用 kinetic 矩构造 Z；GSIS 用混合后的源项宏观量构造 Z。这与参考 Fortran 的监测阶段对应。密度和温度取完整量，不用接近零的扰动作为分母；速度用单元速度积分的模。首次旧值为静止参考气体。')
p('若分母和变化都恰为零，指标记为零；分母为零而变化不为零，记为无穷。当前解满足上述条件就停止，方程残差和壁面诊断只记录，不再额外否决停机。')
table(['流动 / δ','CIS 最终相对更新','GSIS 最终相对更新','CIS / GSIS 步数'],[
    [label(r['flow'],float(r['delta'])),e(r['cisRelativeMacroUpdate']),e(r['gsisRelativeMacroUpdate']),f"{r['cisIterations']} / {r['gsisIterations']}"] for r in main
],[2.2,1.5,1.5,1.5])
p('<b>迭代数比：</b>N_CIS/N_GSIS。它只表示外层迭代次数的比，不是运行时间加速比。每轮 GSIS 有额外的宏观求解，并且其最终离散解尚与 CIS 有差异。')
p('<b>判据的含义：</b>停机阈值 1e-5 是相邻迭代的相对更新要求，不能推出相对 CIS 的误差、相对文献的误差或方程残差都小于 1e-5。')

new('5. 主结果：迭代次数与解的一致性','数值汇总')
table(['算例','δ','CIS','GSIS','迭代数比','速度相对 L2 [%]','剪应力相对 L2 [%]'],[
    [FLOW[r['flow']],f(r['delta']),r['cisIterations'],r['gsisIterations'],f(r['iterationRatio'],4),f(r['relativeVelocityL2Percent'],6),f(r['relativeShearL2Percent'],6)] for r in main
],[1.6,.5,.6,.6,1,1.7,1.7])
sub('绝对差异：保留量纲尺度信息')
table(['算例 / δ','速度绝对 L2 差','最大绝对速度差','原始流量绝对差'],[
    [label(r['flow'],float(r['delta'])),e(r['absoluteVelocityL2']),e(r['maximumAbsoluteVelocityDifference']),e(r['absoluteFlowRateDifference'])] for r in main
],[2.2,1.5,1.5,1.5])
p('δ=0.1 的 Couette/Poiseuille 速度相对 CIS 差异分别为约 6.28e-7% / 2.08e-5%，但两者迭代数均无改善。δ=10 显示减少迭代次数的潜力，同时伴随可测的离散解偏差。')
p('Couette 净流量的绝对差约为机器舍入量级；由于净流量本身为零，其相对流量差不具有稳定意义。不能用这个接近零的积分代替速度剖面的验证。',small=True)

new('6. 原始流量、归一化流量与独立 BGK 参考','数值汇总')
po=[r for r in main if r['flow']=='poiseuille']
table(['δ','CIS 原始 Q','GSIS 原始 Q','CIS Qp','GSIS Qp','流量相对差 [%]'],[
    [f(r['delta']),f(r['cisFlowRate'],11),f(r['gsisFlowRate'],11),f(float(r['cisFlowRate'])*100,10),f(float(r['gsisFlowRate'])*100,10),f(r['relativeFlowDifferencePercent'],6)] for r in po
],[.5,1.6,1.6,1.5,1.5,1.4])
p('上表是当前 64 单元、32×32 速度点、相对更新 1e-5 停机的完整 CIS/GSIS 结果；Qp=100Q。')
sub('速度剖面相对独立 BGK 数值参考的偏差')
table(['算例 / δ','CIS 相对参考 [%]','GSIS 相对参考 [%]','参考自身加密变化 [%]'],[
    [label(r['flow'],float(r['delta'])),f(r['cisVsIndependentBGKPercent'],6),f(r['gsisVsIndependentBGKPercent'],6),e(r['independentBGKSelfRefinementPercent'])] for r in main
],[2.2,1.5,1.5,1.7])
p('独立参考由特征线积分方程求得，使用 2048 空间单元及正速度半区 320 点。自身变化比较 1024/160 与 2048/320 对齐到 DG 单元后的速度。该参考仍有数值误差，不是解析精确解，也不是文献表值。')
p('<b>新增 δ=0.1 的重要含义：</b>GSIS 与 CIS 几乎相同，并不代表当前网格精度充分。其速度相对独立参考偏差约 4.34%（Couette）和 0.640%（Poiseuille），明显大于参考自身加密变化；尚需单独的空间与速度精度研究。')

new('7. 停机时方程残差：下降与最终水平','数值汇总')
table(['算例 / δ','CIS 原方程 R','GSIS 原 CIS 方程 R','GSIS 耦合方程 R'],[
    [label(r['flow'],float(r['delta'])),e(r['cisOriginalResidual']),e(r['gsisOriginalResidual']),e(r['gsisCoupledResidual'])] for r in main
],[2.2,1.5,1.7,1.7])
p('这些是绝对加权残差，不是百分数。GSIS 的原 CIS 残差与耦合残差使用同一份 H，仅碰撞源宏观量不同；两者不必单调保持大小关系。')
sub('相对残差的最终值：共同除以各算例 CIS 首轮 R')
rows=[]
for case in cases:
    r=summary[key(case)];den=case['cis']['history'][0,2]
    rows.append([label(*key(case)),e(float(r['cisOriginalResidual'])/den),e(float(r['gsisOriginalResidual'])/den),e(float(r['gsisCoupledResidual'])/den)])
table(['算例 / δ','CIS R/R1','GSIS 原方程 R/R1','GSIS 耦合方程 R/R1'],rows,[2.2,1.5,1.7,1.7])
p('原 CIS 方程残差确实被 GSIS 显著降低；δ=10 时其下降在后期放缓，Couette 还出现非单调变化。当前数据支持“减少迭代并降低残差”，尚不足以支持“保持同一 CIS 离散稳态且残差达到同样水平”。')

new('8. 空间分辨率：迭代次数与差异','网格对照')
gridrows=[]
for flow in ['couette','poiseuille']:
    for delta in [.1,1,10]:
        grp=sorted([r for r in spatial if r['flow']==flow and float(r['delta'])==delta],key=lambda r:float(r['cells']))
        gridrows.append([label(flow,delta),' / '.join(r['cisIterations'] for r in grp),' / '.join(r['gsisIterations'] for r in grp),
                         ' / '.join(f(r['relativeVelocityL2Percent'],4) for r in grp)])
table(['算例 / δ','CIS：32 / 64 / 128','GSIS：32 / 64 / 128','速度相对差 [%]：32 / 64 / 128'],gridrows,[1.8,1.6,1.6,3])
sub('加密时混合系数也在变')
table(['δ','ω：32 单元','ω：64 单元','ω：128 单元'],[
    [f(delta)]+[f(min(1,2*delta/n/np.sqrt(np.pi)),6) for n in [32,64,128]] for delta in [.1,1,10]
],[.6,1.5,1.5,1.5])
p('CIS 迭代数在本次空间网格范围内不变；δ=0.1、1 的 GSIS 也不变。δ=10 的 GSIS 随网格加密需要更多轮，因此当前实现未表现出网格无关的迭代次数。')
p('与此同时，GSIS 与 CIS 的速度差异下降。由于 ω 与空间离散一起改变，这组结果不是固定混合系数研究，不能把差异下降全部归因于空间精度，也不能仅由三点趋势断言二阶收敛。')

new('9. 文献流量对照：历史数据与当前结果分开','文献对照')
p('文献：F. Sharipov, Non-isothermal gas flow through rectangular microchannels, J. Micromech. Microeng. 9 (1999), 394-401。采用 Table 1（印刷页 398）平行板 a/b=0 的 BGK work [36] 列；不使用方形通道或 S-model 列作为本 BGK 精度基准。',small=True)
table(['δ','历史精细 CIS Qp','文献 BGK Qp','有符号相对差 [%]','联合加密变化 [%]'],[
    [f(r['delta']),f(r['fineQp'],10),r['paperBGK'],f(r['bgkDifferencePercent'],6),e(r['refinementPercent'])] for r in historical
],[.5,1.8,1.4,1.7,1.6])
p('历史精细 CIS：Nx=256、速度 512×8、a_y=1、方程判据 1e-10。联合加密比较 Nx128/Nv256×8 与 Nx256/Nv512×8，不是只加密空间。',small=True)
table(['δ','当前方法','当前 Qp','文献 BGK Qp','有符号相对差 [%]'],[
    [f(r['delta']),r['method'],f(r['Qp'],10),r['paperBGK'],f(r['signedRelativeDifferencePercent'],6)] for r in current
],[.5,1.1,1.8,1.4,1.7])
p('当前结果使用 64 单元、32×32 速度点和相对更新 1e-5。δ=0.1 的流量偏差约 -0.119%；δ=10 的当前 CIS 虽更接近表格值，也不能据此认定比精细网格更准确，需结合独立参考和误差抵消。',small=True)

new('10. δ、Kn 定义与剪切应力的物理检查','归一化与应力')
formula(r'$\delta=\frac{L_{\rm phys}p_0}{\mu_0c_0},\quad c_0=\sqrt{2RT_0},\quad \mathrm{Kn}_{\mu}=\frac{\sqrt{\pi}}{2\delta}$')
p('这里 L_phys 是两平行壁之间的完整间距，采用文献的黏性平均自由程定义。当前无量纲间距 L=1，因此文献 δ 与代码 δ 相同；一般情况下文献 δ=代码碰撞系数×无量纲间距。δ=0.1、1、10 分别对应 Kn_μ≈8.8623、0.88623、0.088623。其他平均自由程约定下的 Kn 不能直接混用。')
sub('Couette：稳态剪切应力近似为常数，但随 δ 改变')
cr=[r for r in diagnostics if r['flow']=='couette']
table(['δ','方法','平均 Pxy','(最大值-最小值)/|平均值| [%]'],[
    [f(r['delta']),r['method'],f(r['meanPxy'],10),e(100*r['relativePxySpan'])] for r in cr
],[.5,.8,1.6,3])
p('常数性是“固定 δ 时沿 x 近似不变”，不是不同 δ 之间数值相同。改变 δ 同时改变碰撞与动量传递，剪切应力随之变化。负号来自 Pxy 为动量通量的约定；对当前正的速度梯度，Pxy 为负。',small=True)
sub('Poiseuille：应力随 x 近似线性')
p('原 kinetic 稳态动量平衡为 dPxy/dx=a_y=0.01，对称解在中心 Pxy≈0。左右真实壁面的通量应接近 -0.005 与 +0.005；单元平均值位于壁内，不能直接要求其最大绝对值等于 0.005。壁面通量及局部动量缺口见附录 B、C。',small=True)

# Figures 01-04 and 14-15. Each page explains all three panels.
for flow,delta,index in [('couette',1,1),('poiseuille',1,2),('couette',10,3),('poiseuille',10,4),('couette',.1,14),('poiseuille',.1,15)]:
    r=summary[(flow,float(delta))];name=f'{index:02d}_{flow}_delta{tag(delta)}_profiles'
    new(f'图 {index:02d}：{label(flow,delta)} 的剖面与差异','逐图说明')
    picture(FIG/f'{name}.png',maxheight=242)
    p('<b>这张图回答：</b>完整 GSIS 与同网格 CIS 的速度剖面是否相同；差异出现在哪些位置；相对独立 BGK 数值参考的准确程度如何。')
    p('<b>左图：</b>横轴 x 为壁面法向坐标，纵轴为单元平均 u_y。蓝实线为 CIS，红虚线为 GSIS，深色点线为独立 BGK 数值参考。仅看重合外观可能漏掉小差异。')
    p('<b>中图：</b>纵轴是有符号的 u_GSIS-u_CIS，保留正负号，表示 GSIS 在哪里更大或更小。<b>右图：</b>纵轴为 100|差值|/max|CIS|；标题是另一项全局相对 L2 差，不能把右图某一点当作全局 L2。')
    p(f'<b>本图数值：</b>速度相对 CIS 差 {f(r["relativeVelocityL2Percent"],6)}%；绝对 L2 差 {e(r["absoluteVelocityL2"])}；最大绝对差 {e(r["maximumAbsoluteVelocityDifference"])}。GSIS 相对独立 BGK 的速度偏差为 {f(r["gsisVsIndependentBGKPercent"],6)}%。')
    if delta==.1:
        p('<b>解读：</b>两算法的差异极小，且没有减少迭代次数；独立参考曲线仍显示当前网格精度的限制。不能用“GSIS 与 CIS 重合”代替准确性验证。',small=True)
    elif flow=='couette':
        p('<b>物理形态：</b>剖面关于中心反对称，气体侧壁面速度与规定壁速不同，体现滑移。图中采样是单元平均值，不是直接把壁面速度接到曲线上。',small=True)
    else:
        p('<b>物理形态：</b>剖面关于中心对称，中心流速较高，静止壁旁气体速度非零体现滑移。δ 较大时趋向抛物线形状，但当前参数下仍不能等同无滑移 NS 极限。',small=True)

# Figures 05/06, grouped by delta for legible row panels.
for delta in [.1,1,10]:
    new(f'图 05 / 06：δ={delta:g} 的绝对更新与相对更新','逐图说明')
    picture(FIG/'pdf_panels'/f'05_absolute_updates_{tag(delta)}.png',maxheight=165)
    picture(FIG/'pdf_panels'/f'06_relative_updates_{tag(delta)}.png',maxheight=165)
    p('<b>上排图 05：</b>横轴为外层迭代数；纵轴为四个节点宏观量绝对 RMS 更新的最大值。蓝线是 CIS kinetic 更新，红线是 GSIS kinetic 更新，黄色点线是 GSIS 混合源更新。它描述“这一轮变化多大”，不描述对参考解的误差。',small=True)
    p('<b>下排图 06：</b>纵轴是实际停机指标 D_ref，蓝线监测 CIS kinetic 矩，红线监测 GSIS 混合源矩；水平虚线为 1e-5。两幅图的左列均为 Couette、右列为 Poiseuille。',small=True)
    if delta in [.1,1]:
        nc=5 if delta==.1 else 11;np_=8 if delta==.1 else 23
        p(f'<b>读图结论：</b>Couette 为 {nc}/{nc} 步，Poiseuille 为 {np_}/{np_} 步，两方法基本重合，未观察到迭代数收益。较小的绝对更新不等同于较高的解精度。',small=True)
    else:
        p('<b>读图结论：</b>GSIS 更早跨过相对更新阈值，Couette 为 82→47 步，Poiseuille 为 251→66 步；但这只是当前耦合算法下的停机改善，还须结合方程残差和最终解差异。',small=True)

# Figures 07/08. Same H, different source definitions are explained each time.
for delta in [.1,1,10]:
    new(f'图 07 / 08：δ={delta:g} 的绝对残差与相对残差','逐图说明')
    picture(FIG/'pdf_panels'/f'07_absolute_residuals_{tag(delta)}.png',maxheight=165)
    picture(FIG/'pdf_panels'/f'08_relative_residuals_{tag(delta)}.png',maxheight=165)
    p('<b>上排图 07：</b>绝对加权 kinetic-DG 方程残差。蓝线为 CIS 的原方程；红线为 GSIS 的原 CIS 方程；黄色虚线为 GSIS 的混合源耦合方程。<b>下排图 08：</b>同样三条曲线，统一除以各算例 CIS 首轮的绝对残差。左右列分别为 Couette、Poiseuille。',small=True)
    p('<b>如何解读：</b>红线与黄线用相同 H、不同碰撞源。残差下降表示方程缺口减小；红黄分离显示两个源项定义的差异。相对残差没有直接采用 1e-5 停机条件，图中也不把该阈值作为其通过线。',small=True)
    r1=summary[('couette',delta)];r2=summary[('poiseuille',delta)]
    p(f'<b>本次终点：</b>Couette 的 GSIS 原 CIS / 耦合残差为 {e(r1["gsisOriginalResidual"])} / {e(r1["gsisCoupledResidual"])}；Poiseuille 为 {e(r2["gsisOriginalResidual"])} / {e(r2["gsisCoupledResidual"])}。绝对值较小不等于相对误差也相同，幅值和残差归一化需要一起看。',small=True)

new('图 09：δ=0.1、1 的空间分辨率与迭代次数','逐图说明')
picture(FIG/'pdf_panels'/'09_spatial_iterations_0p1.png',maxheight=165)
picture(FIG/'pdf_panels'/'09_spatial_iterations_1.png',maxheight=165)
p('<b>目的与坐标：</b>检验空间单元数 32→64→128 时，达到同一停机标准需要多少外层迭代；横轴是空间单元数，纵轴是迭代次数。蓝色圆点为 CIS，红色空心方块为 GSIS；速度始终为 32×32。',small=True)
p('<b>结果：</b>δ=0.1 时 Couette 始终 5/5、Poiseuille 始终 8/8；δ=1 时分别为 11/11 和 23/23。当前这四组的迭代次数在所测网格范围内不变，两方法没有迭代数差别。',small=True)
p('<b>限制：</b>“本次次数不变”只适用于已测范围和当前驱动、速度网格、停机标准；不构成任意网格或任意 δ 下的普遍证明。',small=True)

new('图 09：δ=10 时 GSIS 对空间网格仍敏感','逐图说明')
picture(FIG/'pdf_panels'/'09_spatial_iterations_10.png',maxheight=230)
table(['流动','CIS：32 / 64 / 128','GSIS：32 / 64 / 128','迭代数比：32 / 64 / 128'],[
    ['Couette','82 / 82 / 82','34 / 47 / 59','2.41 / 1.74 / 1.39'],
    ['Poiseuille','251 / 251 / 251','41 / 66 / 102','6.12 / 3.80 / 2.46'],
],[1.2,1.8,1.8,2.2])
p('<b>图在检验什么：</b>若算法具备强的网格无关加速性质，GSIS 达到同样标准的次数应基本稳定。当前 CIS 次数稳定，GSIS 却随加密增加，因此本组没有证明网格无关加速。')
p('<b>同时发生的变化：</b>ω 从 0.3526 降到 0.1763，再降到 0.08815。宏观反馈逐渐变弱，而 GSIS 与 CIS 的解差异也在下降。需要固定 ω 的独立对照，才能分离网格和混合因素。')

new('图 10：δ=0.1、1 的空间加密与解差异','逐图说明')
picture(FIG/'pdf_panels'/'10_spatial_differences_and_mixing_part1.png',maxheight=255)
p('<b>目的：</b>观察网格加密时，GSIS 与 CIS 的离散解差异是否减小。横轴是空间单元数，纵轴是速度剖面相对 CIS 的 L2 百分比，使用对数坐标。蓝线为 Couette，红线为 Poiseuille。')
p('<b>左图 δ=0.1：</b>Couette 的相对差从 1.64e-6% 降到 1.86e-7%；Poiseuille 从 7.73e-5% 降到 5.41e-6%。两算法之间已非常接近。')
p('<b>右图 δ=1：</b>Couette 约 0.003985%→0.001515%→0.000468%；Poiseuille 约 0.02668%→0.009321%→0.002763%。')
p('<b>不能据此推出：</b>这是两算法间差异的下降，不是相对精确解的空间收敛阶；而且 ω 同时变化。δ=0.1 当前剖面仍与独立 BGK 参考有明显差异。')

new('图 10：δ=10 的差异与实际混合系数','逐图说明')
picture(FIG/'pdf_panels'/'10_spatial_differences_and_mixing_part2.png',maxheight=255)
p('<b>左图：</b>δ=10 的速度差随加密减小：Couette 为 0.1226%→0.07736%→0.04630%；Poiseuille 为 0.4047%→0.3003%→0.2069%。但当前范围内仍有可测偏差。')
p('<b>右图：</b>横轴为空间单元数，纵轴为实际 ω。三条线分别是 δ=0.1、1、10，均按 ω=min(1,2δΔx/√π) 随加密减小。它是迭代设置的变化，不是误差。')
p('<b>综合读法：</b>图 09 说明加密时 GSIS 次数增加；图 10 说明同时解差异下降、宏观混合减弱。因此目前观察到的是速度、一致性和混合强度之间的共同变化，尚不能仅靠调小 ω 宣称离散一致性已解决。')

new('图 11：流量在哪里，三组数据分别是什么','逐图说明')
picture(FIG/'pdf_panels'/'11_literature_Qp_flow.png',maxheight=270,width=360)
p('<b>目的与坐标：</b>横轴为 δ（对数坐标），纵轴为归一化 Poiseuille 流量 Qp。它展示文献 BGK 表格与本程序结果的整体对应关系，不是速度剖面或方程残差。')
p('<b>曲线和符号：</b>蓝线为历史精细 CIS（Nx256、Nv512×8、方程阈值 1e-10）；深色空心方块为文献 Table 1 的 BGK work [36]；金色菱形为当前 CIS，红色空心三角为当前 GSIS（均 Nx64、Nv32×32、相对更新阈值 1e-5）。')
p('<b>换算：</b>当前 Qp=100Q，所以 δ=0.1 的 GSIS 原始 Q≈0.02030282 对应 Qp≈2.030282；δ=1 对应约 1.538561；δ=10 对应约 2.759361。不同线性驱动幅值先归一化后才作比较。')

new('图 11：相对文献的绝对差与相对差','逐图说明')
picture(FIG/'pdf_panels'/'11_literature_Qp_differences.png',maxheight=255)
p('<b>左图：</b>有符号绝对差 Qp-Qp,paper；负值表示计算流量小于表格。<b>右图：</b>同一差值除以文献值再乘 100%，用于比较不同 δ 下的相对偏差。颜色和符号与上一页一致。')
p('<b>灰带的含义：</b>仅表示假设文献 Qp 四舍五入到四位小数时的 ±0.00005 区间，经过相对归一化后随 δ 改变。它不是论文报告的误差条，也不能解释所有已有偏差。')
p('<b>当前 GSIS：</b>δ=0.1、1、10 的有符号相对文献差分别约 -0.11895%、-0.04152%、-0.28687%。历史精细 CIS 的最大绝对相对差约 0.07193%。')
p('<b>如何使用：</b>文献对照检验模型和数值结果的外部一致性；GSIS-CIS 对照检验算法间一致性。两种百分比的分母不同，不能互相替换；历史精细 CIS 与当前 GSIS 的网格和停机定义也不同。')

new('图 12：历史文献对照的分辨率检查','逐图说明')
picture(FIG/'12_literature_resolution.png',maxheight=260)
p('<b>左图：</b>纵轴是 100|Qp,fine-Qp,coarse|/|Qp,fine|；联合加密为 Nx128/Nv256×8→Nx256/Nv512×8。因此它验证这两项共同改变后的稳定性，不能当成纯空间收敛阶。')
p('<b>右图：</b>固定 Nx=256、方程停机 1e-10，分别改变法向速度 256→512、512→1024，法向截断 8→10，或流向速度 8→16。每条线仅改一项；百分比与 CSV 中的相对小数相差因子 100。低于 1e-14% 的点在图中画到 1e-14%。')
p('<b>结论：</b>历史联合加密的最大 Qp 相对变化约 0.000241%；单独把法向速度 256 加到 512 的最大相对变化约 5.47e-7%。这些变化远小于多数 δ 下的文献差异。')
p('<b>限制：</b>这说明历史精细网格下，继续上述加密对流量影响很小；不能把它推广成当前 Nv=32 的速度剖面也已收敛，更不能把剩余文献差异直接归因于某一个尚未验证的因素。')

new('图 13：历史 CIS 的速度网格与迭代次数','逐图说明')
picture(FIG/'13_historical_velocity_iterations.png',maxheight=280)
p('<b>目的与坐标：</b>固定 Nx=256、流向速度 8 点、截断 8、单位驱动与方程阈值 1e-10，横轴只改变法向速度点数 256/512/1024；纵轴是 CIS 迭代次数。六条线对应六个 δ。')
p('<b>结果：</b>δ=0.1、0.5、1、2、5、10 的迭代数分别为 15、30、48、89、259、721，均不随上述法向速度加密而改变。')
p('<b>可得结论：</b>这组历史 CIS 对法向速度分辨率没有显示额外的迭代次数增长。<b>不能得到：</b>这不是 GSIS 的速度网格研究，也不是当前相对更新 1e-5 的迭代次数，不能与当前 5、8、11、23 等数值直接横向比较。')

new('附录 A：壁面、质量和碰撞守恒诊断','附加指标')
gs=[r for r in diagnostics if r['method']=='GSIS']
table(['算例 / δ','壁面关系残差','壁面质量通量最大值','平均密度误差','碰撞不变量缺口'],[
    [label(r['flow'],r['delta']),e(r['wallResidual']),e(r['wallMassFlux']),e(r['meanDensityError']),e(r['collisionInvariantDefect'])] for r in gs
],[2.2,1.4,1.8,1.4,1.5])
p('壁面关系残差比较本轮实际使用的入射分布与当前出射分布所要求的入射分布；壁面质量通量按实际 kinetic 数值迹计算，包含平衡背景的舍入误差。平均密度误差是域平均 ρ′ 的绝对值。碰撞不变量缺口检查 H 与离散投影的四个矩是否相等。',small=True)
table(['算例 / δ','max |ρ′|','max |u_x|','max |T′|','max |q_x|','max |q_y|'],[
    [label(r['flow'],r['delta'])]+[e(r[k]) for k in ['maximumRhoPerturbation','maximumUx','maximumTPerturbation','maximumAbsHeatX','maximumAbsHeatY']] for r in gs
],[2,1.1,1.1,1.1,1.1,1.1])
p('这些是最终 kinetic 单元平均量的最大绝对值。接近零的密度、温度、法向速度和法向热流不使用相对零基准误差；q_y 可以是稀薄非平衡热流，不能因等温壁面就强制设零。热流定义 q=∫v(v²-2)H M0 dv，为本 2V 模型的约定。',small=True)

new('附录 B：边界 trace、滑移与壁面应力通量','附加指标')
table(['算例 / δ','左气体 u_y','右气体 u_y','左滑移','右滑移','trace 施加误差'],[
    [label(r['flow'],r['delta']),f(r['leftGasUy'],7),f(r['rightGasUy'],7),f(r['leftSlip'],7),f(r['rightSlip'],7),e(r['macroAppliedTraceError'])] for r in wall
],[2,1.25,1.25,1.25,1.25,1.4])
p('滑移=气体侧速度-规定壁速。边界 trace 施加误差独立重算目标 [ρ′+T′,0,u_y,T′] 后，与 HDG 左右壁 Uhat 比较；本次六组均在舍入误差内。它验证边界取值和施加，不证明两套离散通量完全一致。',small=True)
table(['算例 / δ','左 kinetic Pxy','右 kinetic Pxy','左宏观数值 Pxy','右宏观数值 Pxy'],[
    [label(r['flow'],r['delta']),f(r['leftKineticPxy'],8),f(r['rightKineticPxy'],8),f(r['leftMacroPxy'],8),f(r['rightMacroPxy'],8)] for r in wall
],[2,1.6,1.6,1.7,1.7])
p('kinetic Pxy 用实际壁面分布迹积分；宏观值由带稳定项的外向切向动量数值通量，乘法向符号并除以 2 后换算得到。通量差异可能与体内耦合及离散方式共同相关，不能单凭这一张表判定壁面条件写错。',small=True)

new('附录 C：原 kinetic 通量平衡与宏观线性系统','附加指标')
table(['算例 / δ','质量平衡缺口','流向动量缺口','法向动量缺口','能量平衡缺口'],[
    [label(r['flow'],r['delta'])]+[e(r[k]) for k in ['massBalanceDefect','streamwiseMomentumBalanceDefect','normalMomentumBalanceDefect','energyBalanceDefect']] for r in gs
],[2,1.5,1.5,1.5,1.5])
p('上述值为相邻数值面通量差除以 Δx 后，相对原 kinetic 不变量方程物理右端的最大绝对缺口。流向动量减去 a_y，其他三项的目标导数为零。它与原 CIS 方程检查相关，不等同于 HDG 线性系统残差。',small=True)
table(['算例 / δ','HDG 局部残差','HDG trace 残差','内部面连续残差','HDG 单元守恒残差'],[
    [label(r['flow'],r['delta'])]+[e(r[k]) for k in ['hdgLocalResidual','hdgTraceResidual','hdgFluxContinuityResidual','hdgConservationResidual']] for r in gs
],[2,1.5,1.5,1.5,1.5])
p('HDG 的局部方程、凝聚 trace 线性系统、内部面通量连续与单元平衡均解至舍入误差附近。这说明宏观线性子问题被充分求解；不能据此推出其解与 kinetic-DG 离散稳态必然一致。',small=True)

new('附录 D：单次循环耗时与反馈强度','附加指标')
rows=[]
for case in cases:
    c,g=case['cis'],case['gsis']
    rows.append([label(*key(case)),f(c['wallTime'],5),f(g['wallTime'],5),f(c['wallTime']/g['wallTime'],4),
                 f(np.min(g['synthetic']['omega']),6),e(g['diagnostics']['sourceMomentGap'])])
table(['算例 / δ','CIS 循环秒数','GSIS 循环秒数','单次循环时间比','ω','源矩与 kinetic 矩 RMS 缺口'],rows,[2,1.2,1.2,1.5,1,2])
p('<b>耗时口径：</b>这是各次保存日志中求解循环的 wallTime，未包含全部前处理、速度/空间矩阵构造和所有输出；新 δ=0.1 与旧 δ=1/10 还来自不同批次。未控制预热、缓存或重复测量。')
p('因此不能把这些单次时间比直接称为可信的端到端加速比。迭代数减少是本次已经验证的现象，运行效率还需要统一环境、预热、重复计时和同等精度要求的独立基准。')
p('<b>反馈强度：</b>sourceMomentGap 是四个节点宏观量中，Q_s-Q_kin 的 RMS 最大值。它衡量最终碰撞源偏离 kinetic 自身宏观量的程度，不是速度相对误差，也不是方程残差。')
sub('当前可用于 PPT 的结论')
p('δ=0.1 和 1 时，GSIS 与 CIS 接近，但未减少迭代数；δ=10 时，GSIS 明显减少迭代，同时存在解差异和网格敏感性。小的算法间差异不能替代独立参考精度验证。')

new('附录 E：独立参考 Fortran 程序的实测记录','参考程序')
ref=json.loads((ROOT/'reference_runs/nonoscillatory/original_config/report.json').read_text())
p('这里是原参考代码的 2D3V Shakhov 方腔流，不是新增 δ=0.1 的 1D2V BGK 通道。保留原 control.in：δ≈0.1、578 三角形、P4、24³ 速度点，z 方向只存半空间；默认停机 1e-5、最多 50 步。')
table(['方法','步数','最终更新指标','是否满足 1e-5','程序打印的混合宏观质量'],[
    ['CIS','16','4.880e-6','是','0.99999856'],
    ['GSIS','50','7.605e-5','否，达到步数上限','0.99524733'],
],[1,1,1.5,2.2,2.3])
p('第 50 步 GSIS 的 kinetic 速度向量相对第 16 步 CIS 的单元平均、面积加权 L2 差为 0.02997%；这不是收敛后的 GSIS 稳态误差。')
p('另将参考 CIS 收紧到 1e-9 后启动参考 GSIS，推进两步，kinetic 速度场改变 0.02306%。这个检查说明该参考实现没有严格保持所测 CIS 状态不变；它不等于完整 GSIS 的最终差异。')
p('报告中的当前 δ=0.1 通道则是 Couette 5/5、Poiseuille 8/8 步，二者均已满足各自的当前停机指标。两类算例、模型、网格与迭代状态必须分别表述。')
p('参考程序副本只在测试驱动中增加方法切换、重启和高精度诊断，并跳过昂贵的原插值输出；求解模块保持原文件校验值。该记录用于审查参考程序表现，不用于替代当前通道的文献流量验证。',small=True)

new('附录 F：数据来源、复现方式与适用范围','来源与复现')
sources=[
    ('D1','output/channel_figures_20261109/data/main_results.csv','当前主结果、解差异、独立参考偏差和最终残差。'),
    ('D2','results/channel_presentation_20261109/presentation_data.mat','六组 Nx64 主结果及全部 18 组空间网格病例，含 H、矩、history、宏观量、通量及配置。'),
    ('D3','output/channel_figures_20261109/data/spatial_iterations_and_differences.csv','空间网格、迭代次数、速度差异和 ω。'),
    ('D4','output/channel_figures_20261109/data/historical_literature_comparison.csv','历史精细 CIS 的文献流量对照；不可重标为当前 GSIS 结果。'),
    ('D5','output/channel_figures_20261109/data/historical_velocity_check.csv','单项速度求积加密的 30 组历史结果。'),
    ('D6','output/channel_figures_20261109/data/all_solver_diagnostics.csv','本报告新增汇总的所有求解器诊断；每个主算例含 CIS 与 GSIS。'),
    ('D7','output/channel_figures_20261109/data/all_wall_diagnostics.csv','边界 trace、滑移和两种壁面应力通量。'),
    ('D8','reference_runs/nonoscillatory/original_config/report.json','参考 Fortran 方腔试跑及严格 CIS 重启检查。'),
]
table(['编号','相对项目根目录的文件路径','用途'],sources,[.45,4,2.7],size=8)
p('文献来源：F. Sharipov (1999), Table 1, printed p.398，平行板 a/b=0 的 BGK work [36]；文献中该列来源为 S. K. Loyalka, Phys. Fluids 17 (1974), 1053。本文已核对 δ、Kn 及 Qp 归一化。',small=True)
p('复现：MATLAB 运行 debug/run_channel_presentation_checks.m（δ=1/10），再运行 debug/extend_channel_delta01_presentation.m；Python 运行 debug/plot_channel_presentation.py 生成图片，运行 debug/build_channel_results_pdf.py 生成本报告。当前数值解未因制图或 PDF 汇总而修改。',small=True)
p('适用范围：线性、等温、切向剪切/体力驱动 BGK 通道。尚未完成统一的重复时间基准、固定 ω 的网格研究或证明 GSIS 严格保留 CIS 固定点。剩余文献偏差的来源仍需进一步验证，不能仅凭舍入、停机或单一网格因素解释。',small=True)

C.save()
(TMP/'page_index.json').write_text(json.dumps(toc,ensure_ascii=False,indent=2)+'\n')
print('PDF:',TARGET)
print('Pages:',page)
print('Supplementary CSVs:',DATA/'all_solver_diagnostics.csv',DATA/'all_wall_diagnostics.csv')
