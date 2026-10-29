"""Figure of the path blocks: two adjacent blocks with the straight path between
their roots, a thin path and the slice of links between them, and the torus of
the slice variables with its diagonal loop. Run from this directory:
python3 make_pathblocks.py. Writes fig_pathblocks.png (flattened, 256 colours)."""
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
GREY, BLUE, RED, INK = '#c9c9c9', '#3b6fb6', '#c0392b', '#222222'


def blocks(ax):
    s = 4
    segs = []
    for x, y, z in itertools.product(range(2 * s), range(s), range(s)):
        for d in range(3):
            e = np.eye(3, dtype=int)[d]
            p, q = np.array([x, y, z]), np.array([x, y, z]) + e
            if q[0] < 2 * s and q[1] < s and q[2] < s:
                segs.append([p, q])
    ax.add_collection3d(Line3DCollection(segs, colors=GREY, linewidths=0.5))
    slice_ = [[(s - 1, y, z), (s, y, z)] for y in range(s) for z in range(s)]
    ax.add_collection3d(Line3DCollection(slice_, colors=BLUE, linewidths=2.0))
    ax.plot([0, s], [0, 0], [0, 0], color=INK, lw=2.6)
    thin = [[(j, 2, 3), (j + 1, 2, 3)] for j in range(1, s, 2)]
    ax.add_collection3d(Line3DCollection(thin, colors=RED, linewidths=3.0))
    ax.plot([0, s - 1], [2, 2], [3, 3], color=RED, lw=0.9, ls=(0, (3, 2)))
    for c in [(0, 0, 0), (s, 0, 0)]:
        ax.scatter(*c, color=INK, s=28, depthshade=False)
    ax.text(-0.6, -0.8, -0.9, r'$sy$', color=INK)
    ax.text(s - 0.2, -0.9, -0.9, r'$s(y+e_\mu)$', color=INK)
    ax.text(1.4, -0.6, s + 0.4, r'$B(y)$', color=INK)
    ax.text(s + 2.4, -0.6, s + 0.4, r'$B(y+e_\mu)$', color=INK)
    ax.text(s - 1.2, s + 0.6, -1.4, r'slice $T$', color=BLUE)
    ax.text(-3.4, 2.0, 3.0, 'thin path', color=RED)
    ax.view_init(elev=22, azim=-58)
    ax.set_xlim(0, 2 * s - 1); ax.set_ylim(0, s); ax.set_zlim(0, s)
    clean3d(ax); ax.set_box_aspect((2, 1, 1))


def torus(ax):
    t = np.linspace(0, 2 * np.pi, 200)
    ax.plot([0, 2 * np.pi, 2 * np.pi, 0, 0], [0, 0, 2 * np.pi, 2 * np.pi, 0], color=GREY, lw=1)
    ax.plot(t, t, color=INK, lw=2.2)
    ax.plot(t, 0 * t + 0.02, color=BLUE, lw=2.2)
    ax.plot(0 * t + 0.02, t, color=RED, lw=2.2)
    ax.text(3.5, 1.9, 'diagonal,\ndegree 1', color=INK, fontsize=9)
    ax.text(3.6, 0.25, r'$W_1$', color=BLUE)
    ax.text(0.2, 4.1, r'$W_2$', color=RED)
    ax.text(0.25, 6.75, r'$W_1+W_2+\cdots=1$', color=INK)
    ax.set_xlabel(r'$\theta_1$'); ax.set_ylabel(r'$\theta_2$')
    ax.set_xticks([0, np.pi, 2 * np.pi]); ax.set_xticklabels(['0', r'$\pi$', r'$2\pi$'])
    ax.set_yticks([0, np.pi, 2 * np.pi]); ax.set_yticklabels(['0', r'$\pi$', r'$2\pi$'])
    ax.set_aspect('equal'); ax.set_xlim(-0.3, 2 * np.pi + 0.3); ax.set_ylim(-0.3, 2 * np.pi + 1.0)
    for sp in ax.spines.values(): sp.set_visible(False)


fig = plt.figure(figsize=(7.2, 3.0))
ax1 = fig.add_axes([0.0, 0.0, 0.62, 1.0], projection='3d'); ax1.computed_zorder = False
blocks(ax1)
ax2 = fig.add_axes([0.66, 0.14, 0.32, 0.8])
torus(ax2)
fig.text(0.02, 0.92, '(a)'); fig.text(0.63, 0.92, '(b)')
save(fig, 'fig_pathblocks.png')
