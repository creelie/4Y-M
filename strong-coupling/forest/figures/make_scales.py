"""Figure of the cut-off average across scales on torus-valued fields: a coarse
plaquette as the mean of the 2x2 fluxes at its offsets (drawn in the plane of
the plaquette), and, for the thin-path fields, the coarse energy left by the
k-fold path map and the bound for the k-fold cut-off average. Run from this
directory: python3 make_scales.py. Writes fig_scales.png (flattened, 256
colours)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
GREY, BLUE, RED, INK = '#c9c9c9', '#3b6fb6', '#c0392b', '#222222'


def squares(ax):
    for i in range(4):
        ax.plot([i, i], [0, 3], color=GREY, lw=0.6)
        ax.plot([0, 3], [i, i], color=GREY, lw=0.6)
    ax.add_patch(Rectangle((0, 0), 2, 2, fill=False, ec=INK, lw=2.4))
    cols = ['#3b6fb6', '#5b8fd6', '#2e8b57', '#c0392b']
    for k, (a, b) in enumerate([(0, 0), (1, 0), (0, 1), (1, 1)]):
        ax.add_patch(Rectangle((a + 0.06 * k, b + 0.06 * k), 2, 2, fill=True, fc=cols[k], alpha=0.16,
                               ec=cols[k], lw=1.4))
        ax.plot(a, b, 'o', color=cols[k], ms=4)
    ax.text(3.15, 2.6, 'fine 2x2 squares\nat the offsets o', color=BLUE, fontsize=9)
    ax.text(-0.15, -0.35, r'$2y$', color=INK, fontsize=10)
    ax.text(0.15, 1.55, 'coarse\nplaquette', color=INK, fontsize=8.5)
    ax.text(0.0, -0.75, r'$\Theta^\sharp_P=\frac{1}{16}\sum_o F_{2y+o}$', color=INK, fontsize=10)
    ax.set_xlim(-0.4, 5.2); ax.set_ylim(-1.0, 3.4); ax.set_aspect('equal'); ax.axis('off')


def energies(ax):
    k = np.arange(14, 25)
    s = 2.0 ** k
    ax.semilogy(k, 24 + 0 * k, color=INK, lw=2, label=r'path map, $Q_k=24$')
    ax.semilogy(k, 18 * np.pi ** 2 / s ** 2, color=BLUE, lw=2, label=r'cut-off average, $\leq18\pi^2/s^2$')
    ax.semilogy(k, 24 * s * np.sin(np.pi / (2 * s)) ** 2, color=GREY, lw=1.5, ls=(0, (3, 2)),
                label=r'fine energy $Q$')
    ax.set_xlabel(r'$k$, with $s=2^k$')
    ax.set_ylabel('energy')
    ax.legend(frameon=False, fontsize=8.5, loc='upper right', bbox_to_anchor=(1.0, 0.9))
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)


fig = plt.figure(figsize=(7.2, 3.0))
ax1 = fig.add_axes([0.0, 0.05, 0.42, 0.9]); squares(ax1)
ax2 = fig.add_axes([0.53, 0.16, 0.45, 0.78]); energies(ax2)
fig.text(0.02, 0.92, '(a)'); fig.text(0.44, 0.92, '(b)')
save(fig, 'fig_scales.png')
