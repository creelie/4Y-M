"""Figure for Theorem thm:chessceiling and Corollary cor:chessceiling:
(a) the dyadic scales stacked as plates of side 2^j: the scales at which a
chessboard bound for the islands is below one (green), and the scale of the
cube of D(n, alpha) (red), at least three dyadic scales higher for the labels
of the island theorems and at least twenty for the calm labels;
(b) the condition 16^{j+1} < beta_W eta^2 / (96 (log beta_W + C_0)) of part
(ii): the left side grows by 16 per scale and meets the right side, which does
not depend on j, before the scale of the cube. Schematic, not to scale.
Run from this directory: python3 make_chessceiling.py. Writes fig_chessceiling.png."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1'

fig = plt.figure(figsize=(9.0, 4.0))
ax = fig.add_axes([0.12, 0.0, 0.46, 1.0], projection='3d')
clean3d(ax)
levels = 7
jmax = 3
for j in range(levels):
    side = 0.8 + 0.55 * j
    z = 0.9 * j
    c = -side / 2
    v = [(c, c, z), (c + side, c, z), (c + side, c + side, z), (c, c + side, z)]
    if j <= jmax:
        col, al = GREEN, 0.55
    elif j < levels - 1:
        col, al = GREY, 0.25
    else:
        col, al = RED, 0.55
    ax.add_collection3d(Poly3DCollection([v], facecolor=col, edgecolor=INK, lw=0.5, alpha=al))
    ax.text(c + side + 0.1, c + side / 2, z, r'$j=%s$' % (str(j) if j < levels - 1 else r'j_{\rm cube}'), fontsize=8, color=INK)
fig.text(0.02, 0.25, 'bound below one', color=GREEN, fontsize=9)
fig.text(0.02, 0.56, r'$\geq 3$ scales (islands)' + '\n' + r'$\geq 20$ scales (calm labels)', color=INK, fontsize=9)
fig.text(0.02, 0.80, 'cube of $D(n,\\alpha)$', color=RED, fontsize=9)
ax.set_xlim(-2.5, 2.5); ax.set_ylim(-2.5, 2.5); ax.set_zlim(0, 5.6)
ax.view_init(elev=18, azim=-55)
fig.text(0.02, 0.93, '(a)', fontsize=11)

bx = fig.add_axes([0.6, 0.15, 0.37, 0.75])
j = np.arange(0, 7)
lhs = 16.0 ** (j + 1)
rhs = 16.0 ** 4.4
bx.semilogy(j, lhs, color=BLUE, marker='o', ms=4, label=r'$16^{j+1}$')
bx.axhline(rhs, color=GREEN, lw=1.2, label=r'$\beta_W\tilde\eta^{2}/(96(\log\beta_W+C_0))$')
bx.axvspan(-0.3, 3.4, color=GREEN, alpha=0.12)
bx.axvline(6, color=RED, lw=1.2, ls='--')
bx.text(6.05, 30, r'$j_{\rm cube}$', color=RED)
bx.set_xlabel(r'scale $j$'); bx.set_xlim(-0.3, 6.8)
bx.legend(loc='upper left', fontsize=8, frameon=False)
for sp in ('top', 'right'):
    bx.spines[sp].set_visible(False)
fig.text(0.56, 0.93, '(b)', fontsize=11)
save(fig, 'fig_chessceiling.png')
