"""Figure of the large fields under the SU(2) Wilson measure: (a) a schematic
plane of the lattice with clusters of bad unit cubes and the coarse links on
which the cut-off average can be active, all within distance two of a bad
cube; (b) the frustrated exterior of a link e, one staple at I and five at -I;
(c) the defect weight q_rho(beta_W) of Theorem "defects" for three values of
rho, with the level at which the cluster bound 6400 q <= 1/2 starts. Run from
this directory: python3 make_largefields.py. Writes fig_largefields.png."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
GREY, BLUE, RED, ORANGE, GREEN, INK = '#c9c9c9', '#3b6fb6', '#c0392b', '#e67e22', '#2e8b57', '#222222'
C0 = 0.5 + 2 / 3 * np.log(24 * np.pi ** 3)


def plane(ax):
    n = 24
    for i in range(n + 1):
        ax.plot([i, i], [0, n], color=GREY, lw=0.35)
        ax.plot([0, n], [i, i], color=GREY, lw=0.35)
    bad = [(3, 4), (4, 4), (4, 5), (15, 6), (11, 15), (12, 15), (12, 16), (19, 19), (6, 18)]
    active = set()
    for (x, y) in bad:
        for mu in range(2):
            for yy in range(0, n, 2):
                for xx in range(0, n, 2):
                    tx, ty = x - xx, y - yy
                    t, w = (tx, ty) if mu == 0 else (ty, tx)
                    if 0 <= t <= 2 and 0 <= w <= 1:
                        active.add((xx, yy, mu))
    for (xx, yy, mu) in active:
        dx, dy = (2, 0) if mu == 0 else (0, 2)
        ax.plot([xx, xx + dx], [yy, yy + dy], color=ORANGE, lw=2.2, solid_capstyle='round', zorder=2)
    for (x, y) in bad:
        ax.add_patch(Rectangle((x, y), 1, 1, fc=RED, ec=RED, alpha=0.85, zorder=3))
    for xx in range(0, n + 1, 2):
        for yy in range(0, n + 1, 2):
            ax.plot(xx, yy, '.', color=INK, ms=2.2, zorder=4)
    ax.set_xlim(-0.5, n + 0.5); ax.set_ylim(-5.6, n + 2.6); ax.set_aspect('equal'); ax.axis('off')
    ax.add_patch(Rectangle((0.5, -2.6), 1, 1, fc=RED, ec=RED))
    ax.text(2.8, -2.45, r'bad cube $\mathcal{C}_x$', fontsize=8)
    ax.plot([0.0, 2.0], [-4.6, -4.6], color=ORANGE, lw=2.2)
    ax.text(2.8, -4.95, r'coarse link, bad cube in $\mathcal{Z}_{y,\mu}$', fontsize=8)


def frustrated(ax):
    e0, e1 = np.array([0.0, 0.0]), np.array([2.0, 0.0])
    dirs = [np.array([0, 1.4]), np.array([0, -1.4]), np.array([0.85, 0.75]), np.array([-0.85, -0.75]),
            np.array([-0.8, 0.95]), np.array([0.8, -0.95])]
    for k, d in enumerate(dirs):
        c = GREEN if k == 0 else RED
        p = [e0, e0 + d, e1 + d, e1]
        ax.plot([q[0] for q in p], [q[1] for q in p], color=c, lw=1.6 if k else 2.2, alpha=0.95)
        u = d / np.linalg.norm(d)
        far = (e1 if d[0] >= 0 else e0) + d + 0.32 * u + (np.array([0.18, 0]) if d[0] >= 0 else np.array([-0.18, 0]))
        ax.text(far[0], far[1], r'$I$' if k == 0 else r'$-I$', color=c, ha='center', va='center', fontsize=9)
    ax.plot([0, 2], [0, 0], color=INK, lw=3.2)
    ax.text(1.0, -0.28, r'$e$', ha='center', fontsize=11)
    ax.text(1.0, 1.15, r'$p_0$', color=GREEN, ha='center', fontsize=10)
    ax.text(1.0, -2.35, r'staple sum $-4I$:  $U_e\to-I$,  $p_0$ bad', ha='center', fontsize=9)
    ax.set_xlim(-1.6, 3.6); ax.set_ylim(-2.6, 2.1); ax.set_aspect('equal'); ax.axis('off')


def weights(ax):
    b = np.logspace(1, 6.3, 600)
    for rho, c, lab in [(1.0, INK, r'$\rho=1$'), (0.5, BLUE, r'$\rho=1/2$'), (1 / 12, ORANGE, r'$\rho=1/12$')]:
        lq = np.minimum(0.0, 6 * (np.log(b) + C0) - b * rho ** 2 / 16) / np.log(10)
        ax.plot(b, lq, color=c, lw=1.8, label=lab)
    ax.axhline(np.log10(1 / 12800), color=RED, lw=1, ls=(0, (3, 2)))
    ax.text(12, np.log10(1 / 12800) - 1.0, r'$6400\,q=1/2$', color=RED, fontsize=8.5)
    ax.set_xscale('log'); ax.set_ylim(-14, 1)
    ax.set_xlabel(r'$\beta_W$'); ax.set_ylabel(r'$\log_{10}q_\rho(\beta_W)$')
    ax.legend(frameon=False, fontsize=8.5, loc='lower left', handlelength=1.4, borderaxespad=0.1)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)


fig = plt.figure(figsize=(7.2, 2.75))
ax1 = fig.add_axes([0.0, 0.02, 0.31, 0.95]); plane(ax1)
ax2 = fig.add_axes([0.33, 0.06, 0.27, 0.88]); frustrated(ax2)
ax3 = fig.add_axes([0.71, 0.18, 0.28, 0.76]); weights(ax3)
fig.text(0.01, 0.93, '(a)'); fig.text(0.34, 0.93, '(b)'); fig.text(0.63, 0.93, '(c)')
save(fig, 'fig_largefields.png')
