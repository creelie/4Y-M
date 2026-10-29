"""Figure of the cut-off average: the sixteen paths of length two from a block
of side two, drawn in three of the four dimensions, and, for SU(2), the coarse
link produced by one link exp(theta X) at the straight path or at another
transverse position. Run from this directory: python3 make_cutoff.py. Writes
fig_cutoff.png (flattened, 256 colours)."""
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
GREY, BLUE, RED, INK, GREEN = '#c9c9c9', '#3b6fb6', '#c0392b', '#222222', '#2e8b57'


def paths(ax):
    segs = []
    for x, y, z in itertools.product(range(4), range(2), range(2)):
        for d in range(3):
            q = np.array([x, y, z]) + np.eye(3, dtype=int)[d]
            if q[0] < 4 and q[1] < 2 and q[2] < 2:
                segs.append([(x, y, z), tuple(q)])
    ax.add_collection3d(Line3DCollection(segs, colors=GREY, linewidths=0.6))
    tree = [[(0, 0, 0), (0, 1, 0)], [(0, 0, 0), (0, 0, 1)], [(0, 1, 0), (0, 1, 1)],
            [(0, 0, 0), (1, 0, 0)], [(1, 0, 0), (1, 1, 0)], [(1, 0, 0), (1, 0, 1)], [(1, 1, 0), (1, 1, 1)],
            [(0, 0, 1), (0, 1, 1)]]
    tr = tree[:3] + tree[4:7]
    tr = tr + [[(p[0] + 2,) + p[1:], (q[0] + 2,) + q[1:]] for p, q in tr]
    ax.add_collection3d(Line3DCollection(tr, colors=GREEN, linewidths=1.6, linestyles=(0, (2, 1.5))))
    for o, a1, a2 in itertools.product(range(2), range(2), range(2)):
        h = 0.07 * (1 if o else -1)
        p = [(o, a1, a2 + h), (o + 2, a1, a2 + h)]
        if (o, a1, a2) == (0, 0, 0):
            ax.plot(*zip(*p), color=INK, lw=2.8)
        else:
            ax.plot(*zip(*p), color=BLUE, lw=1.5)
    for c in [(0, 0, 0), (2, 0, 0)]:
        ax.scatter(*c, color=INK, s=26, depthshade=False)
    ax.text(-0.35, -0.05, -0.32, r'$2y$', color=INK)
    ax.text(1.65, -0.1, -0.42, r'$2y+2e_\mu$', color=INK)
    ax.text(0.0, 0.0, 1.25, r'$B(y)$', color=INK)
    ax.text(2.6, 0.0, 1.5, r'$B(y+e_\mu)$', color=INK)
    ax.text(3.25, 1.0, 0.95, r'$P_x$', color=BLUE)
    ax.text(0.75, -0.05, -0.32, r'$V^{\rm p}$', color=INK)
    ax.view_init(elev=24, azim=-62)
    ax.set_xlim(0, 3); ax.set_ylim(0, 1.6); ax.set_zlim(0, 1.6)
    clean3d(ax); ax.set_box_aspect((2.1, 1, 1))


def chi(t):
    return np.clip(2 - 2 * t, 0, 1)


def response(ax):
    th = np.linspace(0, 1.2, 1200)
    dl = 2 * np.sqrt(2) * np.sin(th / 2)          # delta(e^{theta X}, I), rho = 1
    c = chi(dl)
    ax.axvspan(0, 2 * np.arcsin(0.5 / (2 * np.sqrt(2))), color='#eef3fa', lw=0)
    ax.axvspan(2 * np.arcsin(0.5 / (2 * np.sqrt(2))), 2 * np.arcsin(1 / (2 * np.sqrt(2))), color='#f6f6f6', lw=0)
    ax.plot(th, th, color=GREY, lw=1.0, ls=(0, (3, 2)))
    ax.plot(th, th / 8, color=GREY, lw=1.0, ls=(0, (3, 2)))
    ax.plot(th, th - 7 * c * th / 8, color=INK, lw=2.0)
    ax.plot(th, c * th / 8, color=BLUE, lw=2.0)
    ax.text(0.58, 1.0, 'link on the path', color=INK, fontsize=9)
    ax.text(0.76, 0.06, 'link elsewhere', color=BLUE, fontsize=9)
    ax.text(0.02, 1.08, r'$\chi=1$', color=INK, fontsize=9)
    ax.text(0.39, 1.08, r'$0<\chi<1$', color=INK, fontsize=9)
    ax.text(1.03, 0.2, r'$\theta/8$', color='#888888', fontsize=9)
    ax.set_xlabel(r'angle $\theta$ of the link')
    ax.set_ylabel(r'angle of $V^\sharp$')
    ax.set_xlim(0, 1.2); ax.set_ylim(0, 1.2)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)


fig = plt.figure(figsize=(7.2, 3.0))
ax1 = fig.add_axes([0.0, 0.0, 0.56, 1.0], projection='3d'); ax1.computed_zorder = False
paths(ax1)
ax2 = fig.add_axes([0.65, 0.16, 0.33, 0.78])
response(ax2)
fig.text(0.02, 0.92, '(a)'); fig.text(0.56, 0.92, '(b)')
save(fig, 'fig_cutoff.png')
