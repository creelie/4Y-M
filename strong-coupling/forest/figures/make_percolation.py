"""Figure for Theorem thm:clusterpercolation and Lemma lem:windowclusters:
(a) a schematic of the hierarchical construction in three of the four
dimensions: bad blocks (small cubes) grouped into clusters of the first level
inside boxes of side a_1, and clusters of the first level joined into one of
the second level inside a box of side a_2, the joining distance growing like
log of the cluster size; (b) a torus with a window Lambda: a cluster of the
torus that crosses the boundary of the window is cut into clusters of the
window, which are the units of Lemma lem:windowclusters. Schematic, not to scale.
Run from this directory: python3 make_percolation.py. Writes fig_percolation.png."""
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Rectangle
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1'
rng = np.random.default_rng(7)


def cube(ax, o, s, color):
    x, y, z = o
    v = np.array(list(itertools.product((x, x + s), (y, y + s), (z, z + s))))
    faces = [[v[i] for i in f] for f in ((0, 1, 3, 2), (4, 5, 7, 6), (0, 1, 5, 4),
                                        (2, 3, 7, 6), (0, 2, 6, 4), (1, 3, 7, 5))]
    ax.add_collection3d(Poly3DCollection(faces, facecolor=color, edgecolor='#ffffff', lw=0.3, alpha=0.95))


def frame(ax, o, s, color, lw):
    x, y, z = o
    for a, b in itertools.combinations(list(itertools.product((x, x + s), (y, y + s), (z, z + s))), 2):
        if sum(u != w for u, w in zip(a, b)) == 1:
            ax.plot(*zip(a, b), color=color, lw=lw)


fig = plt.figure(figsize=(9.0, 3.9))
ax = fig.add_axes([-0.03, 0.08, 0.55, 0.94], projection='3d')
clean3d(ax)
A1, A2 = 4.0, 12.0
centres = [(1.5, 1.5, 1.5), (7.5, 2.0, 1.8), (4.0, 7.5, 6.0), (8.0, 8.0, 8.0)]
cols = [BLUE, GREEN, '#8e44ad', '#d68910']
for c, col in zip(centres, cols):
    o = np.array(c) - A1 / 2 + 2
    frame(ax, o, A1, col, 0.8)
    for _ in range(5):
        p = o + rng.integers(0, 3, 3) + 0.2
        cube(ax, p, 0.75, col)
frame(ax, (0, 0, 0), A2, RED, 1.4)
for a, b in zip(centres[:-1], centres[1:]):
    ax.plot(*zip(np.array(a) + 2, np.array(b) + 2), color=RED, lw=1.0, ls='--')
fig.text(0.03, 0.08, r'thin boxes: clusters of level 1, side $a_1$', color=INK, fontsize=9)
fig.text(0.03, 0.03, r'red box and dashes: one cluster of level 2, side $a_2$', color=RED, fontsize=9)
ax.view_init(20, -60)
ax.set_xlim(0, A2); ax.set_ylim(0, A2); ax.set_zlim(0, A2)
ax.set_title(r'(a) clusters joined across levels', fontsize=10, y=0.97)

bx = fig.add_axes([0.58, 0.08, 0.4, 0.8])
bx.set_aspect('equal'); bx.axis('off')
bx.add_patch(Rectangle((0, 0), 10, 10, fill=False, ec=INK, lw=1.0))
bx.add_patch(Rectangle((2.5, 2.5), 5, 5, fc='#eaf1fb', ec=BLUE, lw=1.4))
path = np.array([[0.8, 6.0], [2.0, 6.2], [3.0, 5.6], [4.0, 6.0], [5.0, 5.4], [6.2, 5.8], [7.0, 5.2],
                 [8.2, 5.6], [9.2, 5.0]])
for (x0, y0), (x1, y1) in zip(path[:-1], path[1:]):
    inside = 2.5 <= (x0 + x1) / 2 <= 7.5
    bx.plot([x0, x1], [y0, y1], color=GREEN if inside else GREY, lw=2.0)
for x, y in path:
    inside = 2.5 <= x <= 7.5
    bx.add_patch(Rectangle((x - 0.22, y - 0.22), 0.44, 0.44, fc=GREEN if inside else GREY, ec='none'))
small = [[3.4, 3.9], [3.9, 3.8], [6.4, 4.1]]
for x, y in small:
    bx.add_patch(Rectangle((x - 0.22, y - 0.22), 0.44, 0.44, fc='#8e44ad', ec='none'))
bx.text(2.6, 7.75, r'window $\Lambda$', color=BLUE, fontsize=9)
bx.text(0.2, 9.3, r'torus $(\mathbb{Z}/L\mathbb{Z})^4$', color=INK, fontsize=9)
bx.text(0.15, 4.0, 'a torus\ncluster', color=GREY, fontsize=9)
bx.text(3.0, 6.8, r'its part in $\Lambda$', color=GREEN, fontsize=9)
bx.text(2.75, 2.85, 'other window clusters', color='#8e44ad', fontsize=8)
bx.set_xlim(-0.2, 10.2); bx.set_ylim(-0.2, 10.2)
bx.set_title(r'(b) clusters of a window', fontsize=10)
save(fig, 'fig_percolation.png')
