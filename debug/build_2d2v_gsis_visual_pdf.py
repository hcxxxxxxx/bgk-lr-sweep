"""Two-page visual replacement: counts, residuals, velocity and vortex plots.

Reads saved current runs; never extrapolates velocity to impose no-slip walls.
Vortex markers are the elliptic zeros exported by plotCavityCorners2D2V.
"""
from pathlib import Path
import csv
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/gsis_inventory_matplotlib')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import landscape, A4
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'results/cavity_gsis_gray_defaults_2d2v'
TMP = ROOT/'tmp/pdfs/2d2v_gsis_visual_20261112'
OUT = ROOT/'output/pdf/2D2V_GSIS_迭代残差与速度涡旋_两页图版_20261112.pdf'
TMP.mkdir(parents=True, exist_ok=True)
OUT.parent.mkdir(parents=True, exist_ok=True)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


summary = read(DATA/'summary.csv')
plt.rcParams.update({'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 8,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 7, 'axes.spines.top': False,
                     'axes.spines.right': False})
assets = []
for row in summary:
    d = int(float(row['delta']))
    folder = DATA/f'delta_{d}_N32x32_stretch2_GSIS'
    h = read(folder/'history.csv')
    p = np.genfromtxt(folder/'profile.csv', delimiter=',', names=True)
    roots = np.genfromtxt(folder/'vortex_candidates.csv', delimiter=',', names=True)
    x = np.unique(p['x']); y = np.unique(p['y'])
    nx, ny = len(x), len(y)
    assert (nx, ny) == (32, 32)
    U = p['ux'].reshape(ny, nx)/.01
    V = p['uy'].reshape(ny, nx)/.01
    iu = RegularGridInterpolator((y, x), U)
    iv = RegularGridInterpolator((y, x), V)
    paths = {}

    def figure(name):
        f, ax = plt.subplots(figsize=(5.5, 2.95), constrained_layout=True)
        paths[name] = TMP/f'delta_{d}_{name}.png'
        return f, ax

    def save(f, name):
        f.savefig(paths[name], dpi=220, facecolor='white')
        plt.close(f)

    f, ax = figure('residual')
    steps = [int(float(r['iteration'])) for r in h]
    for field, name, color, style in [
        ('macroRelativeUpdate', 'Macro update', '#1867A8', 'o-'),
        ('relativeKineticResidual', 'Original kinetic residual', '#B34638', 's-'),
        ('relativeCoupledResidual', 'Coupled residual', '#3D816F', '^--'),
    ]:
        ax.semilogy(steps, [float(r[field]) for r in h], style, color=color,
                    markersize=3.5, linewidth=1.5, label=name)
    ax.axhline(1e-5, color='#8796A3', linestyle=':', linewidth=1, label='1e-5')
    ax.set(xlabel='Outer iteration', ylabel='Update / relative residual', ylim=(1e-6, 1e2))
    ax.grid(alpha=.18); ax.legend(loc='upper right')
    save(f, 'residual')

    f, ax = figure('velocity')
    speed = np.hypot(U, V)
    im = ax.contourf(x, y, speed, levels=np.linspace(0, 1, 25), cmap='viridis', extend='max')
    ax.quiver(x[::3], y[::3], U[::3, ::3], V[::3, ::3],
              color='white', scale=7, width=.0035, headwidth=3.5)
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel='x', ylabel='y', aspect='equal')
    f.colorbar(im, ax=ax, label='|u| / U_lid', ticks=[0, .25, .5, .75, 1])
    save(f, 'velocity')

    def streams(ax, bounds, density=1.4):
        # Uniform interpolation only WITHIN sampled cell-centre coordinates.
        a, b, c, e = bounds
        xx = np.linspace(max(a, x[0]), min(b, x[-1]), 260)
        yy = np.linspace(max(c, y[0]), min(e, y[-1]), 260)
        X, Y = np.meshgrid(xx, yy)
        uv = iu(np.stack([Y, X], axis=-1)); vv = iv(np.stack([Y, X], axis=-1))
        ax.streamplot(xx, yy, uv, vv, color='#2167A0', density=density,
                      linewidth=.7, arrowsize=.8)
        ax.set(xlim=(a, b), ylim=(c, e), aspect='equal', xlabel='x', ylabel='y')

    f, ax = figure('vortex')
    streams(ax, (0, 1, 0, 1), density=1.5)
    ax.annotate('', xy=(.82, 1), xytext=(.18, 1),
                arrowprops={'arrowstyle': '->', 'color': '#BA493A', 'lw': 2})
    main = roots[np.argmin(roots['vorticity'])]
    ax.plot(main['x'], main['y'], 'o', color='#B34638', ms=4)
    ax.annotate(f"CW\n({main['x']:.3f}, {main['y']:.3f})",
                xy=(main['x'], main['y']), xytext=(.07, .53), fontsize=8,
                arrowprops={'arrowstyle': '-', 'color': '#B34638'}, color='#973D32',
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': .9, 'pad': 2})
    save(f, 'vortex')

    f, axes = plt.subplots(1, 2, figsize=(5.5, 2.95), constrained_layout=True)
    paths['corners'] = TMP/f'delta_{d}_corners.png'
    bottom = roots[(roots['y'] < .12) & ((roots['x'] < .12) | (roots['x'] > .88))]
    for ax, bounds, label in zip(axes, [(0, .12, 0, .12), (.88, 1, 0, .12)],
                                  ['Bottom left', 'Bottom right']):
        streams(ax, bounds, density=1.6)
        selected = bottom[(bottom['x'] > bounds[0]) & (bottom['x'] < bounds[1])]
        for r in selected:
            ax.plot(r['x'], r['y'], 'o', mfc='white', mec='#B34638', ms=5, mew=1.3)
            ax.text(.5, .91, f"{'CCW' if r['vorticity'] > 0 else 'CW'} candidate",
                    ha='center', va='top', transform=ax.transAxes, color='#973D32', fontsize=7.5)
        ax.set_title(label)
        ax.set_xticks([bounds[0], (bounds[0]+bounds[1])/2, bounds[1]])
        ax.set_yticks([0, .06, .12])
    save(f, 'corners')
    assets.append({'row': row, 'paths': paths, 'bottom_candidates': len(bottom),
                   'primary_center': [float(main['x']), float(main['y'])]})

pdfmetrics.registerFont(TTFont('CN', '/System/Library/Fonts/STHeiti Light.ttc', subfontIndex=1))
pdfmetrics.registerFont(TTFont('CNB', '/System/Library/Fonts/STHeiti Medium.ttc', subfontIndex=1))
W, H = landscape(A4)
M = 26
GAP = 18
COL = (W - 2*M - GAP)/2
PANEL = 216
C = canvas.Canvas(str(OUT), pagesize=(W, H))
C.setTitle('2D2V GSIS：迭代、残差、速度场与涡旋（两页图版）')


def panel(path, left, bottom, title):
    C.setFont('CNB', 10)
    C.setFillColor(HexColor('#17324A'))
    C.drawString(left, bottom+PANEL-11, title)
    iw, ih = Image.open(path).size
    width = COL; height = width*ih/iw
    if height > PANEL-19:
        height = PANEL-19; width = height*iw/ih
    C.drawImage(str(path), left+(COL-width)/2, bottom, width=width, height=height,
                preserveAspectRatio=True, mask='auto')


for page, entry in enumerate(assets, 1):
    r = entry['row']; d = int(float(r['delta'])); paths = entry['paths']
    C.setFillColor(HexColor('#1867A8')); C.rect(0, H-6, W, 6, fill=1, stroke=0)
    C.setFillColor(HexColor('#17324A')); C.setFont('CNB', 18)
    C.drawString(M, H-32, f'2D2V GSIS 方腔  |  δ={d}')
    C.setFillColor(HexColor('#1867A8')); C.setFont('CNB', 20)
    C.drawRightString(W-M, H-32, f"{r['iterations']} 次迭代")
    C.setFillColor(HexColor('#50687A')); C.setFont('CN', 8.8)
    C.drawString(M, H-48, '当前配置：空间 32×32，速度 32×32；S=[1,1,1,1]，ω=1，无内部残差修正；宏观更新阈值 1e-5。')
    C.setFont('CN', 8.5)
    C.drawString(M, H-64, f"宏观更新 {float(r['referenceMacroUpdate']):.3e}    原方程残差 {float(r['relativeOriginalResidual']):.3e}    耦合残差 {float(r['relativeCoupledResidual']):.3e}")
    C.setFillColor(HexColor('#A24536'))
    C.drawRightString(W-M, H-64, '宏观更新已达标；原方程残差未达 1e-5')
    top = H-83-PANEL; low = top-9-PANEL
    panel(paths['residual'], M, top, '① 残差与更新历史')
    panel(paths['velocity'], M+COL+GAP, top, '② 速度场：速度大小与方向')
    panel(paths['vortex'], M, low, '③ 整体流线与主涡')
    panel(paths['corners'], M+COL+GAP, low, '④ 底角涡结构放大')
    C.setStrokeColor(HexColor('#D8E2EA')); C.line(M, 39, W-M, 39)
    C.setFillColor(HexColor('#50687A')); C.setFont('CN', 7.6)
    if d == 100:
        old = '历史迭代数（不同配置，S=[1,δ,δ,δ]）：CIS 2803；GSIS 无修正 18；八矩修正、ω=0.1：100。'
    else:
        old = '历史迭代数（S=[1,δ,δ,δ]，ω=1，无内部修正）：32×32 空间 9 次；64×64 空间 18 次。'
    C.drawString(M, 27, old)
    C.drawString(M, 14, '图取自已保存的实际 kinetic 速度；红圈为插值涡心候选。CW=顺时针，CCW=逆时针。')
    C.drawRightString(W-M, 14, f'2026-11-12  |  {page}/2')
    C.showPage()
C.save()
assert len(PdfReader(OUT).pages) == 2
(TMP/'figures.json').write_text(json.dumps([
    {'delta': a['row']['delta'], 'iterations': a['row']['iterations'],
     'primary_center': a['primary_center'], 'bottom_candidates': a['bottom_candidates'],
     'source': str(DATA/f"delta_{a['row']['delta']}_N32x32_stretch2_GSIS")}
    for a in assets], indent=2)+'\n')
print(OUT)
