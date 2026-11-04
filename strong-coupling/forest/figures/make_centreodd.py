"""Figure for Theorem thm:centreodd and Lemma lem:boxenergy: (a) the spatial
torus (Z/LZ)^3 of the box drawn as a cube with identified faces, the winding
loops W_1 and W_2, and the plane x_1 = const across which the centre map
zeta_1 multiplies the links of direction 1 by -1, so that Tr W_1 and
Tr W_1 W_2 change sign; (b) the energy bound eps(b) = 3(log b + 5)/(2b) of
Lemma lem:boxenergy against the threshold 1/(1600 L^4) of
Theorem thm:centreodd(ii) at L = 3.
Run from this directory: python3 make_centreodd.py. Writes fig_centreodd.png."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1'

fig = plt.figure(figsize=(9.0, 3.8))
ax = fig.add_axes([-0.02, 0.06, 0.52, 0.96], projection='3d')
clean3d(ax)
ax.computed_zorder = False
L = 3
for i in range(L + 1):
    for j in range(L + 1):
        for (a, b_, c) in (((0, L), i, j), (i, (0, L), j), (i, j, (0, L))):
            xs = a if isinstance(a, tuple) else (a, a)
            ys = b_ if isinstance(b_, tuple) else (b_, b_)
            zs = c if isinstance(c, tuple) else (c, c)
            ax.plot(xs, ys, zs, color=GREY, lw=0.5, zorder=1)
plane = [[(1.5, 0, 0), (1.5, L, 0), (1.5, L, L), (1.5, 0, L)]]
ax.add_collection3d(Poly3DCollection(plane, facecolor='#f6d5d0', edgecolor=RED, alpha=0.5, lw=1.0))
ax.plot([0, L], [1, 1], [2, 2], color=BLUE, lw=2.4, zorder=5)
ax.plot([1, 1], [0, L], [2, 2], color=GREEN, lw=2.4, zorder=5)
ax.scatter([1.5], [1], [2], color=RED, s=28, zorder=6)
ax.text(L + 0.1, 1, 2.15, r'$W_1$', color=BLUE, fontsize=11, zorder=7)
ax.text(1, L + 0.15, 2.2, r'$W_2$', color=GREEN, fontsize=11, zorder=7)
fig.text(0.04, 0.02, r'shaded plane: $\zeta_1$ multiplies the links of direction 1 crossing it by $-1$', color=RED, fontsize=9)
ax.view_init(22, -62)
ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
ax.set_title(r'(a) the box $(\mathbb{Z}/L\mathbb{Z})^3$ and the winding loops', fontsize=10, y=0.97)

bx = fig.add_axes([0.6, 0.17, 0.37, 0.7])
b = np.logspace(0, 9, 400)
eps = 3 * (np.log(b) + 5) / (2 * b)
thr = 1 / (1600 * L ** 4)
bx.loglog(b, eps, color=BLUE, lw=1.6, label=r'$\varepsilon(b)=3(\log b+5)/(2b)$')
bx.axhline(thr, color=RED, lw=1.2, ls='--')
bstar = b[np.argmax(eps <= thr)]
bx.axvline(bstar, color=GREY, lw=0.8, ls=':')
bx.text(1.5, thr * 1.6, r'$1/(1600L^4)$, $L=3$', color=RED, fontsize=9)
bx.set_xlabel(r'coupling $b$')
bx.spines[['top', 'right']].set_visible(False)
bx.legend(frameon=False, fontsize=9, loc='upper right')
bx.set_title(r'(b) the plaquette energy of the box', fontsize=10)
save(fig, 'fig_centreodd.png')
print('first b on the grid with eps <= threshold:', '%.3g' % bstar)
