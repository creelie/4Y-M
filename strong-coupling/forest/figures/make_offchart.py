"""Figure of Theorem "offchart": (a) a fibre of the cut-off average drawn as a surface, with the
product chart around the background U, the small-field points inside it and the points with a large
plaquette outside; (b) the least beta_W at which the bound 10^4 |W| q_kappa(beta_W) of the theorem
is at most 1/2, as a function of kappa, for the 81 blocks of a cube of side three (one block gives the same curve to plotting accuracy).
Run from this directory: python3 make_offchart.py. Writes fig_offchart.png."""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1'

fig = plt.figure(figsize=(9.0, 3.6))
ax = fig.add_axes([-0.04, -0.02, 0.56, 1.04], projection='3d')
clean3d(ax)
ax.computed_zorder = False

def height(x, y):
    return 0.45 * np.sin(0.8 * x) * np.cos(0.6 * y) + 0.08 * (x ** 2 - y ** 2) / 4

n = 121
t = np.linspace(-3, 3, n)
X, Y = np.meshgrid(t, t)
Z = height(X, Y)
cols = np.empty(X.shape + (4,))
base = np.array(matplotlib.colors.to_rgba('#dfe6ee'))
cols[:] = base
chart = (np.abs(X) <= 1.25) & (np.abs(Y) <= 1.25)
cols[chart] = matplotlib.colors.to_rgba('#b9e0c6')
small = X ** 2 + Y ** 2 <= 0.8 ** 2
cols[small] = matplotlib.colors.to_rgba('#8fb3e6')
bad = np.zeros_like(X, bool)
for cx, cy, r in ((2.2, -1.9, 0.45), (-2.1, 2.0, 0.55), (2.4, 1.6, 0.35), (-2.3, -2.2, 0.3)):
    bad |= (X - cx) ** 2 + (Y - cy) ** 2 <= r ** 2
cols[bad] = matplotlib.colors.to_rgba('#efb1a8')
ax.plot_surface(X, Y, Z, facecolors=cols, rstride=2, cstride=2, lw=0.15, edgecolor='#ffffff',
                shade=False, zorder=1)
s = np.linspace(-1.25, 1.25, 50)
for xs, ys in ((s, 0 * s - 1.25), (s, 0 * s + 1.25), (0 * s - 1.25, s), (0 * s + 1.25, s)):
    ax.plot(xs, ys, height(xs, ys) + 0.01, color=GREEN, lw=1.4, zorder=3)
th = np.linspace(0, 2 * np.pi, 120)
cx, cy = 0.8 * np.cos(th), 0.8 * np.sin(th)
ax.plot(cx, cy, height(cx, cy) + 0.012, color=BLUE, lw=1.3, zorder=3)
ax.scatter([0], [0], [height(0, 0) + 0.02], color=INK, s=14, zorder=4)
ax.text(0.15, -0.1, height(0, 0) + 0.25, r'$U$', fontsize=11, zorder=5)
ax.text(1.7, -0.6, height(1.3, -0.6) + 0.15, r'$\Phi(\Omega_r)$', color=GREEN, fontsize=10, zorder=5)
ax.text(-2.5, 2.0, height(-2.1, 2.0) + 0.45, r'large plaquette', color=RED, fontsize=9, zorder=5)
ax.view_init(32, -58)
ax.set_xlim(-3, 3); ax.set_ylim(-3, 3); ax.set_zlim(-1.2, 1.2)
ax.set_title(r'(a) the fibre $F$ and the chart', fontsize=10, y=0.93)

C0 = 0.5 + 2 / 3 * math.log(24 * math.pi ** 3)


def beta_min(kappa, nw):
    """Least beta >= 1 with 10^4 nw q_kappa(beta) <= 1/2, q = exp(6(log b + C0) - b kappa^2/16)."""
    target = math.log(2e4 * nw)
    lo, hi = 1.0, 1.0
    f = lambda b: b * kappa ** 2 / 16 - 6 * (math.log(b) + C0) - target
    while f(hi) < 0:
        hi *= 2
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if f(mid) >= 0:
            hi = mid
        else:
            lo = mid
    return hi


ax2 = fig.add_axes([0.60, 0.17, 0.37, 0.70])
ks = np.logspace(-19, 0, 200)
ax2.plot(ks, [beta_min(k, 81) for k in ks], color=RED, lw=1.6, label=r'$|W|=81$')
ax2.axvline(4.84e-18, color=GREEN, ls='--', lw=1.0)
ax2.text(7.5e-18, 1e5, r'$\kappa=4.8\cdot10^{-18}$', color=GREEN, fontsize=8, rotation=90, va='bottom')
ax2.set_xscale('log'); ax2.set_yscale('log')
ax2.set_xlabel(r'$\kappa$')
ax2.set_ylabel(r'least $\beta_W$ with bound $\leq 1/2$')
ax2.legend(frameon=False, fontsize=8, loc='upper right')
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
ax2.set_title(r'(b) the bound $10^4|W|\,q_\kappa(\beta_W)$', fontsize=10)
print('least beta at kappa = 4.8e-18, |W| = 81:', f'{beta_min(4.8e-18, 81):.3e}')
save(fig, 'fig_offchart.png')
