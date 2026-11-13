"""Three-page inventory of completed 2D2V BGK / GSIS tests, from saved data.

No solver runs, no reinterpretation of old stopping criteria as current ones.
Run with reportlab available (or PYTHONPATH=/private/tmp/gsis_pdf_report_deps).
"""
from pathlib import Path
import csv
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/gsis_inventory_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.colors import HexColor, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, Table, TableStyle, Image
from reportlab.lib.styles import ParagraphStyle
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/pdf'
TMP = ROOT / 'tmp/pdfs/2d2v_gsis_inventory_20261112'
OUT.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)
TARGET = OUT / '2D2V_BGK_GSIS_已完成测试清单_20261112.pdf'


def read(relative):
    with (ROOT / relative).open() as f:
        return list(csv.DictReader(f))


sources = {
    'S1': 'results/cavity_gsis_gray_defaults_2d2v/summary.csv',
    'S2': 'results/cavity_gsis_comparison_delta100/summary.csv',
    'S3': 'results/cavity_gsis_delta1000_reference/summary.csv',
    'S4': 'results/cavity_2d2v/summary.csv',
    'S5': 'results/cavity_gsis_comparison_delta100/one_macro.csv',
    'S6': 'results/cavity_delta100_cis_tol1e5/delta_100_N32x32_stretch2/corner_summary.csv',
    'S7': 'results/gray_reference_platform/summary.csv',
    'S8': 'results/su2019_bgk_couette/refinement_summary.csv',
}
latest = read(sources['S1'])
old100 = read(sources['S2'])
old1000 = read(sources['S3'])
cis = read(sources['S4'])
one = read(sources['S5'])
cis100 = read(sources['S6'])[0]
gray = read(sources['S7'])
assert len(latest) == 2 and len(old100) == 2 and len(old1000) == 2
assert all(float(r['relativeOriginalResidual']) > 1e-5 for r in latest)

# Plot only the latest configuration: the three curves have distinct meanings.
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.15), constrained_layout=True)
for ax, row in zip(axes, latest):
    d = int(float(row['delta']))
    path = f'results/cavity_gsis_gray_defaults_2d2v/delta_{d}_N32x32_stretch2_GSIS/history.csv'
    h = read(path)
    x = [int(float(v['iteration'])) for v in h]
    for field, label, color, style in [
        ('macroRelativeUpdate', 'Macro update', '#1867A8', 'o-'),
        ('relativeKineticResidual', 'Original equation', '#B34436', 's-'),
        ('relativeCoupledResidual', 'Coupled equation', '#4B8170', '^--'),
    ]:
        ax.semilogy(x, [float(v[field]) for v in h], style, color=color,
                    markersize=2.5, linewidth=1.2, label=label)
    ax.axhline(1e-5, color='#8895A2', linestyle=':', linewidth=1)
    ax.set(title=f'delta = {d}', xlabel='Outer iteration', ylim=(1e-6, 100))
    ax.grid(alpha=.16)
    ax.legend(fontsize=6.5, loc='upper right')
    ax.tick_params(labelsize=7)
    ax.set_ylabel('Update / relative residual', fontsize=7)
chart = TMP / 'latest_histories.png'
fig.savefig(chart, dpi=210, facecolor='white')
plt.close(fig)

pdfmetrics.registerFont(TTFont('CN', '/System/Library/Fonts/STHeiti Light.ttc', subfontIndex=1))
pdfmetrics.registerFont(TTFont('CNB', '/System/Library/Fonts/STHeiti Medium.ttc', subfontIndex=1))
pdfmetrics.registerFontFamily('CN', normal='CN', bold='CNB', italic='CN', boldItalic='CNB')
INK = HexColor('#17324A')
MUTED = HexColor('#566B7E')
BLUE = HexColor('#1867A8')
PALE = HexColor('#EDF4F9')
AMBER = HexColor('#95601E')
W, H = A4
LEFT = 38
WIDTH = W - 2 * LEFT
BODY = ParagraphStyle('body', fontName='CN', fontSize=9.4, leading=14.1,
                      textColor=INK, wordWrap='CJK')
SMALL = ParagraphStyle('small', parent=BODY, fontSize=8.1, leading=11.7, textColor=MUTED)
CELL = ParagraphStyle('cell', parent=BODY, fontSize=8.5, leading=11.6)
HEAD = ParagraphStyle('head', parent=CELL, fontName='CNB', textColor=BLUE)
C = canvas.Canvas(str(TARGET), pagesize=A4)
C.setTitle('2D2V BGK / GSIS 已完成测试清单')
C.setAuthor('CIS / GSIS numerical test inventory')
page = 0
y = 0
page_index = []


def new(title, subtitle):
    global page, y
    if page:
        C.showPage()
    page += 1
    page_index.append({'page': page, 'title': title})
    C.setFillColor(BLUE)
    C.rect(0, H - 9, W, 9, fill=1, stroke=0)
    C.setFont('CNB', 17)
    C.setFillColor(INK)
    C.drawString(LEFT, H - 42, title)
    C.setFont('CN', 8.8)
    C.setFillColor(MUTED)
    C.drawString(LEFT, H - 61, subtitle)
    C.setStrokeColor(HexColor('#D9E3EB'))
    C.line(LEFT, 31, W - LEFT, 31)
    C.setFont('CN', 7.6)
    C.drawString(LEFT, 19, '2026-11-12  |  依据已保存结果与检查脚本整理；本报告未重跑流动算例')
    C.drawRightString(W - LEFT, 19, f'{page} / 3')
    y = H - 83


def p(text, small=False, gap=6):
    global y
    obj = Paragraph(text, SMALL if small else BODY)
    _, height = obj.wrap(WIDTH, 1000)
    assert y - height > 39, (page, 'paragraph overflow', y, height, text[:60])
    obj.drawOn(C, LEFT, y - height)
    y -= height + gap


def sub(text):
    global y
    y -= 3
    C.setFont('CNB', 10.8)
    C.setFillColor(BLUE)
    C.drawString(LEFT, y - 10, text)
    y -= 22


def table(headers, rows, fractions, size=8.5):
    global y
    cs = ParagraphStyle('cs', parent=CELL, fontSize=size, leading=size+3)
    hs = ParagraphStyle('hs', parent=HEAD, fontSize=size, leading=size+3)
    content = [[Paragraph(str(v), hs) for v in headers]]
    content += [[Paragraph(str(v), cs) for v in row] for row in rows]
    widths = [WIDTH * f / sum(fractions) for f in fractions]
    t = Table(content, colWidths=widths, hAlign='LEFT')
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), PALE),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [white, HexColor('#F8FAFC')]),
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('LEFTPADDING', (0, 0), (-1, -1), 6),
        ('RIGHTPADDING', (0, 0), (-1, -1), 6),
        ('TOPPADDING', (0, 0), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
        ('LINEBELOW', (0, 0), (-1, 0), .55, HexColor('#B9CCDB')),
        ('LINEBELOW', (0, -1), (-1, -1), .4, HexColor('#DCE5EC')),
    ]))
    _, height = t.wrap(WIDTH, 1000)
    assert y - height > 39, (page, 'table overflow', y, height)
    t.drawOn(C, LEFT, y - height)
    y -= height + 8


def image(path, height):
    global y
    im = Image(str(path))
    # Resolve lazy intrinsic dimensions before assigning display dimensions.
    iw, ih = im.imageWidth, im.imageHeight
    im.drawWidth = WIDTH
    im.drawHeight = WIDTH * ih / iw
    if im.drawHeight > height:
        ratio = height / im.drawHeight
        im.drawWidth *= ratio
        im.drawHeight = height
    assert y - im.drawHeight > 39, (page, 'image overflow')
    im.drawOn(C, LEFT + (WIDTH-im.drawWidth)/2, y-im.drawHeight)
    y -= im.drawHeight + 6


def sci(value, digits=3):
    return f'{float(value):.{digits}e}'


new('2D2V BGK / GSIS：已完成测试清单', '1  当前配置与正式复测  |  二维空间 (x,y) + 二维速度 (vx,vy)')
p('<b>范围：</b>正文列二维顶盖方腔 GSIS，以及支撑它的 CIS 和离散检查。1D2V Couette/Poiseuille、Su 论文及其他物理模型的参考程序单独注明，不作为二维精度验证。')
table(['项目', '实际设置'], [
    ['物理问题', '单位正方形；顶壁 ux=0.01，其余三壁静止；四壁等温、完全漫反射；无体力。线性扰动模型，不包含有限雷诺数的非线性惯性效应。'],
    ['空间与速度', '矩形 Q1-DG，每单元 4 节点；宏观 Q1-HDG。正式方腔速度均为 32×32 三次映射点，两个速度截断 ±4，零附近加密。'],
    ['当前宏观耦合', '四个宏观量回传碰撞源，不改写 H；稳定参数 [1,1,1,1]，无额外内部弱残差载荷，回传系数 ω=1。高阶项沿用本次扫描的壁面入射量，仿射壁面通量另由扫描后 H 构造。'],
    ['当前停机与网格', '单元积分宏观量逐分量相对更新 D&lt;1e-5；32×32 空间、meshStretch=2。原方程与耦合残差持续记录，默认不参与停机。'],
], [1, 5.4])
sub('A. 最新 Gray 结构适配版：两组正式复测 [S1]')
table(['δ', '步数', '宏观更新 D', '原方程 r_o', '耦合 r_c', '速度相对 CIS'], [
    [r['delta'], r['iterations'], sci(r['referenceMacroUpdate']), sci(r['relativeOriginalResidual']),
     sci(r['relativeCoupledResidual']), f"{float(r['velocityDifferenceVsCISPercent']):.5f}%" if r['velocityDifferenceVsCISPercent'] != 'NaN' else '无同配置 CIS 对照']
    for r in latest
], [.45, .48, 1.15, 1.15, 1.15, 1.5])
p('<b>判定：</b>两组达到宏观更新停止条件；原 kinetic 方程均未达到 1e-5。δ=100 的速度差异 0.15622% 是对已保存、有限阈值 CIS 的差异，不是文献误差。宏观壁面通量关系误差约 1e-19，仅反映子问题求解。')
image(chart, 166)
p('<b>图的含义：</b>蓝线为停机用的宏观更新；红线把 H 自身速度矩代回原 BGK-DG 方程；绿线代入宏观回传源。虚水平线为 1e-5。三者并非同一指标，更新停止不代表红线收敛。', small=True)
p('<b>指标：</b>r_o、r_c 均为速度加权、空间逆质量矩阵 L2 残差，除以各自算例共同的 CIS 首轮原方程残差。速度差异采用单元平均速度向量的相对 L2 范数。', small=True)
p('[S1] results/cavity_gsis_gray_defaults_2d2v/summary.csv；同目录各算例 history.csv 和 cavity.mat 保存逐步指标与完整分布。', small=True)

new('历史流动测试：设置、结果与局限', '2  历史版本分别列出  |  不与当前默认版本混算加速比')
p('<b>历史参数：</b>下列 GSIS 均采用 S=[1,δ,δ,δ]。δ=100 的两版还使用旧壁面高阶项取值时序及旧停机规则；与第 1 页同时改变了稳定参数、时序和停机，不能将步数变化归因于单个因素。')
sub('B. δ=100：同网格 CIS / 两版 GSIS 对照 [S2]')
table(['方法 / 宏观回传', '步数', '原方程相对残差', '速度相对 CIS', '原方程 &lt;1e-5'], [
    ['CIS 基准', cis100['cisIterations'], sci(cis100['relativeKineticResidual']), '基准', '是'],
    ['GSIS：无内部修正，ω=1', old100[0]['gsisIterations'], sci(old100[0]['relativeOriginalResidual']), f"{float(old100[0]['velocityDifferencePercent']):.5f}%", '否'],
    ['GSIS：八矩残差修正，ω=0.1', old100[1]['gsisIterations'], sci(old100[1]['relativeOriginalResidual']), f"{float(old100[1]['velocityDifferencePercent']):.5f}%", '是'],
], [2.35, .55, 1.5, 1.2, 1.0])
p('共用 32×32 空间、stretch=2、32×32 速度、tol=1e-5。无修正版按更新和耦合残差停机；修正版还检查原方程。三者均观察到一对底角结构；只有修正版额外匹配了原 kinetic 的八个弱矩残差，它不属于当前参考默认算法。', small=True)
p(f'<b>一次宏观解检查 [S5]：</b>已停机 CIS → 一次 HDG 宏观解，速度变化：无修正版 {float(one[0]["oneMacroVelocityDifferencePercent"]):.4f}%，修正版 {float(one[1]["oneMacroVelocityDifferencePercent"]):.5f}%。CIS 本身仍有有限残差，不能将这组值当作精确固定点测试。', small=True)
p('<b>探索性试跑：</b>δ=100、旧 S 的八矩修正版在 ω=1 时曾发散；ω=0.1 的有效结果见上表。布局修正前的试跑与未收敛探索不列为正式精度结果。', small=True)
sub('C. δ=1000：仅 GSIS 的空间 / 角涡探索 [S3]')
table(['空间 / stretch', '最小网格宽度', '步数', '原方程 r_o', '耦合 r_c', '底角候选数'], [
    [f"{r['cellsPerDirection']}×{r['cellsPerDirection']} / {r['meshStretch']}", sci(r['minimumCellWidth']),
     r['iterations'], sci(r['relativeOriginalResidual']), sci(r['relativeCoupledResidual']), r['bottomCornerCandidates']]
    for r in old1000
], [1.6, 1.25, .5, 1.1, 1.1, .9])
p('无内部修正、ω=1、32×32 速度、tol=1e-5；两组原方程均未通过阈值。64×64 网格出现每侧两个候选角涡，但第二级流函数量级约 6.1e-11，小于全域通量投影误差约 6.1e-8，尚未确认其物理精度。此处没有 δ=1000 CIS 基准。', small=True)
sub('D. 支撑用的纯 CIS 测试，不计作 GSIS 精度验证 [S4/S6]')
table(['δ', '空间', '步数', '原方程相对残差', '观察 / 网格检查'], [
    ['1', '均匀 24² / 32²', '22 / 22', '约 1.64e-8', '顺时针主涡；24→32 中心线变化 &lt;1%'],
    ['10', '均匀 24² / 32²', '160 / 160', '约 9.19e-8', '顺时针主涡；24→32 中心线变化 &lt;0.18%'],
    ['100', '32² / stretch=2', '2803', sci(cis100['relativeKineticResidual']), '主涡及两底角小涡；tol=1e-5'],
], [.45, 1.45, 1.05, 1.45, 2.65])
p('δ=1/10 使用 tol=1e-7；δ=100 使用 tol=1e-5。以上均未完成二维速度网格独立性研究。', small=True, gap=3)
p('[S2] cavity_gsis_comparison_delta100/summary.csv；[S3] cavity_gsis_delta1000_reference/summary.csv；[S4] cavity_2d2v/summary.csv；[S5] 同 S2 的 one_macro.csv；[S6] cavity_delta100_cis_tol1e5/.../corner_summary.csv。路径均相对 results/。', small=True)

new('离散检查与验证边界', '3  通过了什么检查  |  哪些结论尚未建立')
p('下表来自已执行的三个 debug 检查脚本。小网格检查用于排查实现错误；它们验证的是所选离散算子，不能替代连续解精度、外部文献或任意 δ 下的稳定性证明。', small=True)
table(['已完成检查', '测试设置 / 结果', '脚本'], [
    ['整体 kinetic 方程与四方向扫描', '2×2，4×4 GH，δ=2：直接整体解与 CIS 的 H 最大差约 4.4e-13；四象限显式扫描与逐速度 LU 最大差约 7.1e-15。', 'V1'],
    ['守恒、壁面与对称性', '零法向质量通量、碰撞不变量、反向驱动线性、静止态、方腔对称性通过；均匀及加密网格的 GMRES/CIS 对照通过。', 'V1'],
    ['二维 HDG 独立弱式组装', '3×2 非均匀网格，8×8 GH，δ=3；独立 Gauss 积分对照，当前组装最大差约 5.6e-17。', 'V2'],
    ['四壁 Jacobian 与凝聚求解', '四个法向、四个宏观扰动的漫反射 Jacobian 有限差分检查；局部/全局残差及凝聚/整体解对照通过。', 'V2'],
    ['宏观数组布局与源项回传', '逐单元/分量核对节点顺序；首轮 H 与 CIS 一致；宏观量完整回传，不直接提升或覆盖 H。', 'V2/V3'],
    ['高阶项的壁面取值时序', '3×2、δ=1000：单独扰动所存扫描入射量，HOT 面矩变化与独立积分一致，而仿射壁面数据不变；δ=3/1000 均检查复用扫描入射量。', 'V3'],
    ['当前默认参数与停机分母', 'S 全 1、无内部修正、ω=1；逐分量单元积分及近零分母分支通过。4×4、δ=10、8×8 GH 默认试算 9 步停止，原残差约 0.0608，属功能试算。', 'V3'],
    ['可选八矩残差修正与固定点', '独立 kinetic 弱残差取八矩对照通过；4×4、δ=10：一次宏观变化 &lt;1e-7；GSIS 与严格 CIS 的宏观相对差 &lt;1e-6。ω=1/0.1 分别 15/96 步，均检查原/耦合残差。', 'V2'],
], [1.45, 4.55, .55], size=8.1)
p('<b>复现脚本：</b>V1=debug/check_cavity_2d2v.m；V2=debug/check_gsis_cavity_2d2v.m；V3=debug/check_gray_defaults_cavity_2d2v.m。V2 的固定点及 15/96 步结果仅属于显式开启八矩修正的测试版本。', small=True)
sub('关联研究：记录过，但不能算作本模型的二维验证')
p('<b>参考程序：</b>NonOscillatory 原方腔是不同的 3V 气体配置，延续计算见 results/reference_cavity_residual_platform/，GSIS 有残差平台及密度漂移，未称充分收敛。最新 Gray ZIP 是 Callaway 声子模型：200 三角形、P3、80×40 角度点；默认第 8 步停止，延长到 100 步仍有原方程残差。初态归一化约 4.51e-4，参考自身 RHS 归一化约 8.22e-6，二者不能混用。[S7]', small=True)
p('<b>1D2V 文献对照：</b>Su 等 2019 年 BGK Couette 表 1 已做 29 次 CIS/GSIS、阈值和分辨率测试；加密 GSIS 的三项剪切量最大偏差约 0.0204%。它是特定 1D 剪切子问题的对照，原文为 3V、算法为 SIS；不是同模型/同算法复现，也不验证二维方腔。[S8]', small=True)
p('数据索引：[S7] results/gray_reference_platform/summary.csv；[S8] results/su2019_bgk_couette/refinement_summary.csv。', small=True)
sub('当前尚未完成')
p('① 当前全 1 稳定参数版的大 δ 原 kinetic 方程充分收敛；② 本模型二维方腔与匹配文献/高精度独立解的对照；③ 二维 GSIS 的系统空间、速度及阈值独立性与阶数验证；④ δ=1000 更高阶角涡的可靠确认。', small=True)
p('<b>结论边界：</b>已验证多项组装、边界和耦合实现；已有历史修正版在 δ=100 上达到原方程阈值。但当前参考默认版只达到宏观更新停机，不能把“图有主涡/角涡、宏观线性残差小、迭代快”合并成“整个 GSIS 已验证正确”。', small=True)

assert page == 3
C.save()
reader = PdfReader(TARGET)
assert len(reader.pages) == 3
for i, pg in enumerate(reader.pages):
    assert len(pg.extract_text()) > 200, i
manifest = {'date': '2026-11-12', 'pages': page_index,
            'source_sha256': {k: hashlib.sha256((ROOT/v).read_bytes()).hexdigest() for k, v in sources.items()},
            'scope': '2D2V BGK cavity GSIS; related 1D and other-model references explicitly separated',
            'solver_rerun': False}
(TMP/'inventory_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+'\n')
print(TARGET)
print('Pages:', len(reader.pages))
