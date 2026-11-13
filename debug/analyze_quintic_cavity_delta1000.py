"""Independent bilinear velocity-root audit and raw-velocity corner plots."""
from pathlib import Path
import os
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/quintic_cavity_mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'results/cavity_delta1000_quintic_20261112'


def velocity_roots(x, y, U, V):
    """Eliminate one coordinate analytically; inspect all real roots per patch."""
    found = []
    for j in range(len(y)-1):
        for i in range(len(x)-1):
            u, v = U[j:j+2, i:i+2], V[j:j+2, i:i+2]
            if u.min() > 0 or u.max() < 0 or v.min() > 0 or v.max() < 0:
                continue
            def coeff(q):
                return np.array([q[0, 0], q[0, 1]-q[0, 0], q[1, 0]-q[0, 0],
                                 q[1, 1]-q[1, 0]-q[0, 1]+q[0, 0]])
            a, b = coeff(u), coeff(v)
            scale = max(abs(a).max(), abs(b).max())
            if scale == 0:
                continue
            a, b = a/scale, b/scale
            poly = [b[1]*a[3]-b[3]*a[1],
                    b[0]*a[3]+b[1]*a[2]-b[2]*a[1]-b[3]*a[0],
                    b[0]*a[2]-b[2]*a[0]]
            for z in np.roots(np.trim_zeros(poly, 'f')):
                if abs(z.imag) > 1e-9:
                    continue
                s = float(z.real)
                c = a if abs(a[2]+a[3]*s) > abs(b[2]+b[3]*s) else b
                den = c[2]+c[3]*s
                if abs(den) < 1e-15:
                    continue
                t = -(c[0]+c[1]*s)/den
                if not (0 <= s <= 1 and 0 <= t <= 1):
                    continue
                f = np.array([1, s, t, s*t])
                if max(abs(a@f), abs(b@f)) > 1e-8:
                    continue
                J = scale*np.array([[a[1]+a[3]*t, a[2]+a[3]*s],
                                     [b[1]+b[3]*t, b[2]+b[3]*s]])
                J /= [x[i+1]-x[i], y[j+1]-y[j]]
                discriminant = np.trace(J)**2-4*np.linalg.det(J)
                found.append([x[i]+s*(x[i+1]-x[i]), y[j]+t*(y[j+1]-y[j]),
                              J[1, 0]-J[0, 1], discriminant, float(discriminant < 0)])
    return np.array(found).reshape(-1, 5)


def analyze(case):
    d = np.genfromtxt(case/'profile.csv', delimiter=',', names=True)
    x, y = np.unique(d['x']), np.unique(d['y'])
    U, V = (d[k].reshape(len(y), len(x)) for k in ('ux', 'uy'))
    roots = velocity_roots(x, y, U, V)
    np.savetxt(case/'independent_velocity_roots.csv', roots, delimiter=',',
               header='x,y,vorticity,discriminant,elliptic', comments='')
    elliptic = roots[roots[:, 4] == 1]
    print(case.name, 'elliptic roots:\n', elliptic, flush=True)
    iu, iv = (RegularGridInterpolator((y, x), z) for z in (U, V))
    fig, axes = plt.subplots(3, 2, figsize=(10, 13), constrained_layout=True)
    for row, extent in enumerate([.08, .004, .0004]):
        for side in (0, 1):
            ax = axes[row, side]
            # Use distance to each corner on the horizontal axis.
            xx = np.linspace(x[0], extent, 350)
            yy = np.linspace(y[0], extent, 350)
            X, Y = np.meshgrid(xx, yy)
            pts = np.stack([Y, X if side == 0 else 1-X], axis=-1)
            uu, vv = iu(pts), iv(pts)
            if side == 1:
                uu = -uu
            ax.streamplot(xx, yy, uu, vv, density=1.7, color='#2167a0',
                          linewidth=.7, arrowsize=.7)
            distance = elliptic[:, 0] if side == 0 else 1-elliptic[:, 0]
            selected = (distance < extent) & (elliptic[:, 1] < extent)
            ax.plot(distance[selected], elliptic[selected, 1], 'o', ms=6,
                    mec='#b34638', mfc='white', mew=1.5)
            inside = x[x < extent]
            ax.plot(inside, np.zeros_like(inside), '|', color='gray', ms=5)
            ax.set(xlim=(0, extent), ylim=(0, extent), aspect='equal',
                   xlabel='Distance from side wall', ylabel='y',
                   title=f"{'Left' if side == 0 else 'Right'} bottom corner; extent={extent:g}")
            ax.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
            ax.text(.03, .96, f'{selected.sum()} elliptic velocity zero(s)',
                    transform=ax.transAxes, va='top', color='#973d32', fontsize=9,
                    bbox=dict(facecolor='white', alpha=.9, edgecolor='none'))
    fig.suptitle(f'GSIS delta=1000, {len(x)}x{len(y)} quintic grid; raw kinetic velocity\n'
                 'Bilinear interpolation of cell-centre data; right corner shown in reflected coordinates')
    fig.savefig(case/'corner_multiscale_independent.png', dpi=180)
    plt.close(fig)
    h = np.genfromtxt(case/'history.csv', delimiter=',', names=True)
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    for k, label in [('macroRelativeUpdate', 'Returned macro update'),
                     ('relativeKineticResidual', 'Original equation'),
                     ('relativeCoupledResidual', 'Coupled equation')]:
        ax.semilogy(h['iteration'], h[k], label=label)
    first = np.flatnonzero(h['macroRelativeUpdate'] < 1e-5)
    if len(first):
        ax.axvline(h['iteration'][first[0]], ls=':', color='gray', label='First native stop')
    ax.set(xlabel='Outer iteration', ylabel='Update / relative residual', title=case.name)
    ax.grid(alpha=.2);ax.legend()
    fig.savefig(case/'residual_history.png', dpi=180);plt.close(fig)


if __name__ == '__main__':
    for case in sorted(BASE.glob('delta_1000_N*_quintic_GSIS')):
        if (case/'profile.csv').exists():
            analyze(case)
