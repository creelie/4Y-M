"""Figure of several islands of bad blocks in a window: the islands drawn as
columns over a grid of blocks, the good layers between them, and the pairwise
cross-ratio bounds Delta_uv of Theorem thm:manyislands; and the size of Delta
for two islands against the width w of the good layer between them, with the
prefactor of Lemma lem:lambdatwo. Run from this directory:
python3 make_manyislands.py. Writes fig_manyislands.png (flattened, 256
colours)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
GREY, BLUE, RED, GREEN, INK = '#c9c9c9', '#3b6fb6', '#c0392b', '#2e8b57', '#222222'

ISLANDS = {
    'a': ([(1, 1), (2, 1), (1, 2)], RED),
    'b': ([(12, 11), (13, 11), (13, 12), (12, 12)], RED),
    'q_1': ([(7, 3), (8, 3)], BLUE),
    'q_2': ([(5, 10), (5, 11), (6, 11)], BLUE),
    'q_3': ([(12, 5)], BLUE),
}


def window(ax):
    n = 15
    for i in range(n):
        for j in range(n):
            ax.bar3d(i, j, 0, 0.92, 0.92, 0.04, color='#eef3ea', edgecolor='#d4dccf', linewidth=0.2, shade=False)
    centres = {}
    for name, (cells, col) in ISLANDS.items():
        for (i, j) in cells:
            ax.bar3d(i, j, 0.04, 0.92, 0.92, 0.3, color=col, alpha=0.9, edgecolor=INK, linewidth=0.3)
        c = np.mean(cells, axis=0) + 0.46
        centres[name] = c
        ax.text(c[0] + 1.6, c[1] - 2.4, 0.05, '$I_{%s}$' % name, color=col, fontsize=12, ha='center', zorder=12)
    pairs = [('a', 'b', RED, r'$\Delta_{ab}$', 2.0), ('a', 'q_1', BLUE, '', 1.0), ('q_1', 'b', BLUE, '', 1.0),
             ('a', 'q_2', BLUE, '', 1.0), ('q_2', 'b', BLUE, '', 1.0), ('q_3', 'b', BLUE, '', 0.8),
             ('a', 'q_3', BLUE, '', 1.0)]
    t = np.linspace(0, 1, 60)
    for u, v, col, lab, h in pairs:
        p, q = centres[u], centres[v]
        x = p[0] + (q[0] - p[0]) * t
        y = p[1] + (q[1] - p[1]) * t
        z = 0.34 + 4 * h * t * (1 - t)
        ax.plot(x, y, z, zorder=10, color=col, lw=1.8 if lab else 0.9, ls='-' if lab else (0, (3, 2)))
        if lab:
            ax.text(x[30], y[30], z[30] + 0.25, lab, color=col, fontsize=12, ha='center')
    ax.text2D(0.02, 0.0, 'good layer: every other block is good', transform=ax.transAxes, fontsize=8.5,
              color=GREEN)
    ax.text2D(0.02, 0.07, r'dashed: the couplings $\Delta_{uq}$ through the other islands', transform=ax.transAxes,
              fontsize=8.5, color=BLUE)
    ax.view_init(elev=34, azim=-62)
    clean3d(ax)
    ax.set_zlim(0, 2.6)
    ax.set_title('(a) islands and their pairwise couplings', fontsize=10)


def decay(ax):
    aG = 2.5e-4
    w = np.linspace(10, 1.3e5, 600)
    D = 3.25e21 * np.exp(-2 * aG * (w - 10))
    ax.semilogy(w / 1e3, D, color=RED, lw=2, label=r'$\Delta(w)$ with the evaluated $\Lambda_2$')
    ax.axhline(1, color=INK, lw=0.8, ls=(0, (3, 2)))
    ax.axvline(99.072, color=GREY, lw=1.2)
    ax.text(100.5, 1e12, r'$w=99072$', color=INK, fontsize=9)
    ax.text(5, 3, r'$\Delta=1$', color=INK, fontsize=9)
    ax.set_xlabel(r'width $w$ of the good layer (thousands of blocks)')
    ax.set_ylabel(r'cross-ratio bound')
    ax.set_ylim(1e-8, 1e23)
    ax.set_title('(b) the bound for two islands', fontsize=10)
    ax.legend(frameon=False, fontsize=8.5, loc='upper right')


fig = plt.figure(figsize=(10, 4.2))
ax1 = fig.add_subplot(1, 2, 1, projection='3d', computed_zorder=False)
window(ax1)
ax2 = fig.add_subplot(1, 2, 2)
decay(ax2)
fig.tight_layout()
save(fig, 'fig_manyislands.png')
