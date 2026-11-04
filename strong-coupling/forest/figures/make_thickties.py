"""Figure for Proposition prop:thickties and Lemma lem:ties:
(a) two branch values f_{x_a}, f_{x'_a} over two coordinates of xi when
u_a = w_a = u with d(u, I) = s: over the open set U_0 they differ by at most
2s, so the near-ties form a set of positive measure, not a thin band;
(b) the set U_0 in the coordinates delta(xi_a, I) and max_{b != a}
delta(xi_b, I), and, for comparison, a thin band of width t around a tie
hypersurface, the shape that a bound Haar(T(t)) <= M t would need.
Schematic, not to scale.
Run from this directory: python3 make_thickties.py. Writes fig_thickties.png."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1'

fig = plt.figure(figsize=(9.0, 3.9))
ax = fig.add_axes([-0.02, 0.02, 0.55, 0.96], projection='3d')
clean3d(ax)
x, y = np.meshgrid(np.linspace(0, 1, 70), np.linspace(0, 0.8, 60))
dx = np.maximum(np.maximum(17 / 32 - x, x - 31 / 32), 0)
dy = np.maximum(y - 15 / 32, 0)
dist = np.hypot(dx, dy)
w = 1 - np.exp(-(dist / 0.12) ** 2)
f1 = 0.55 + 0.45 * x - 0.25 * y ** 2
s = 0.02
f2 = f1 + s * np.sin(9 * x) * np.cos(7 * y) + w * 0.28
ax.plot_surface(x, y, f1, color=BLUE, alpha=0.55, lw=0, rstride=2, cstride=2)
ax.plot_surface(x, y, f2, color=GREEN, alpha=0.45, lw=0, rstride=2, cstride=2)
xs = np.array([17 / 32, 31 / 32, 31 / 32, 17 / 32, 17 / 32])
ys = np.array([0, 0, 15 / 32, 15 / 32, 0])
ax.plot(xs, ys, np.full(5, 0.35), color=RED, lw=1.4)
for xx, yy in zip(xs[:4], ys[:4]):
    ax.plot([xx, xx], [yy, yy], [0.35, 0.55 + 0.45 * xx - 0.25 * yy ** 2], color=RED, lw=0.6, ls=':')
ax.text(0.62, 0.05, 0.25, r'$U_0$', color=RED, fontsize=11)
ax.text(0.0, -0.05, 0.42, r'$f_{x_a}$', color=BLUE, fontsize=11)
ax.text(0.0, 0.8, 1.0, r'$f_{x^\prime_a}$', color=GREEN, fontsize=11)
fig.text(0.36, 0.33, r'$|f_{x_a}-f_{x^\prime_a}|\leq 2s$ over $U_0$', color=INK, fontsize=9)
ax.set_zlim(0.35, 1.2)
ax.view_init(elev=16, azim=-70)
fig.text(0.02, 0.93, '(a)', fontsize=11)

bx = fig.add_axes([0.6, 0.15, 0.37, 0.75])
bx.add_patch(Rectangle((17 / 32, 0), 14 / 32, 15 / 32, facecolor=RED, alpha=0.25, edgecolor=RED))
bx.text(0.66, 0.2, r'$U_0\subset T(t)$', color=RED)
tt = np.linspace(0, 1, 100)
curve = 0.75 - 0.35 * (tt - 0.1) ** 2
bx.fill_between(tt[:40], curve[:40] - 0.025, curve[:40] + 0.025, color=BLUE, alpha=0.35)
bx.plot(tt[:40], curve[:40], color=BLUE, lw=1)
bx.text(0.02, 0.83, r'band of width $t$', color=BLUE, fontsize=9)
bx.set_xlim(0, 1); bx.set_ylim(0, 1)
bx.set_xlabel(r'$\delta(\xi_a,I)$'); bx.set_ylabel(r'$\max_{b\neq a}\delta(\xi_b,I)$')
bx.set_xticks([0, 17 / 32, 31 / 32]); bx.set_xticklabels(['0', r'$\frac{17}{32}$', r'$\frac{31}{32}$'])
bx.set_yticks([0, 15 / 32]); bx.set_yticklabels(['0', r'$\frac{15}{32}$'])
for sp in ('top', 'right'):
    bx.spines[sp].set_visible(False)
fig.text(0.56, 0.93, '(b)', fontsize=11)
save(fig, 'fig_thickties.png')
