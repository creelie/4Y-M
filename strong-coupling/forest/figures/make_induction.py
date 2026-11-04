"""Figure for Theorem thm:allscales and Proposition prop:scalecount: (a) the
tori of the scales k = 0, ..., 4 drawn as stacked grids, the statement P_k on
each, the condition (C) on each step, and the steps whose coupling lies below
the first-step large-field threshold lambda = 809300; (b) the lower bound
lower bound log10 log10 Lambda_k >= log10 log10 48^(6^k) of Theorem
thm:allscales(iv), against the uniform bound sigma_k <= 0.4865 of (ii).
Run from this directory: python3 make_induction.py. Writes fig_induction.png."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import proj3d
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1'

fig = plt.figure(figsize=(9.0, 3.9))
ax = fig.add_axes([-0.06, -0.04, 0.6, 1.08], projection='3d')
clean3d(ax)
ax.set_box_aspect((1, 1, 1.25))
ax.computed_zorder = False
K = 5
for k in range(K):
    z = 1.15 * k
    n = 16 // 2 ** k
    s = np.linspace(-2, 2, n + 1)
    col = RED if k >= 3 else BLUE
    face = '#f6d5d0' if k >= 3 else '#d6e3f5'
    X, Y = np.meshgrid([-2, 2], [-2, 2])
    ax.plot_surface(X, Y, z + 0 * X, color=face, alpha=0.55, shade=False, zorder=2 * k)
    for v in s:
        ax.plot([v, v], [-2, 2], [z, z], color=col, lw=0.45, alpha=0.8, zorder=2 * k + 1)
        ax.plot([-2, 2], [v, v], [z, z], color=col, lw=0.45, alpha=0.8, zorder=2 * k + 1)
    if k < K - 1:
        ax.quiver(0, 0, z + 0.12, 0, 0, 0.9, color=GREEN, lw=1.3, arrow_length_ratio=0.25, zorder=20)
ax.view_init(18, -55)
ax.set_xlim(-2.4, 2.4); ax.set_ylim(-2.4, 2.4); ax.set_zlim(0, 1.15 * (K - 1) + 0.5)
ax.set_title(r'(a) the scales and the condition (C)', fontsize=10, y=0.97)
fig.canvas.draw()
M = ax.get_proj()
for k in range(K):
    z = 1.15 * k
    corners = [proj3d.proj_transform(cx, cy, z, M)[:2] for cx in (-2, 2) for cy in (-2, 2)]
    lx, ly = min(corners)
    ax.annotate(r'$P_{%d}$' % k, xy=(lx, ly), xytext=(-16, -4), textcoords='offset points',
                color=INK, fontsize=10, zorder=50)
    if k < K - 1:
        ax_, ay_ = proj3d.proj_transform(0, 0, z + 0.55, M)[:2]
        ax.annotate(r'(C)', xy=(ax_, ay_), xytext=(6, -3), textcoords='offset points',
                    color=GREEN, fontsize=9, zorder=50)
fig.text(0.01, 0.1, r'blue: $\lambda_k\geq809300$', color=BLUE, fontsize=9)
fig.text(0.01, 0.05, r'red: $\lambda_k<809300$', color=RED, fontsize=9)

bx = fig.add_axes([0.6, 0.17, 0.37, 0.7])
k = np.arange(0, 5)
loglogL = np.log10(np.log10(48.0)) + k * np.log10(6.0)
bx.plot(k, loglogL, 'o-', color=RED, lw=1.4, ms=4, label=r'$\log_{10}\log_{10}48^{6^k}\leq\log_{10}\log_{10}\Lambda_k$')
bx.axhline(0.4865, color=GREEN, lw=1.2, ls='--')
bx.text(2.6, 0.6, r'$\sigma_k\leq0.4865$ for all $k$', color=GREEN, fontsize=9, zorder=50)
bx.set_xlabel(r'scale $k$')
bx.set_xticks(k)
bx.spines[['top', 'right']].set_visible(False)
bx.legend(frameon=False, fontsize=9, loc='upper left')
bx.set_title(r'(b) $\sigma_k$ stays bounded, $\Lambda_k$ does not', fontsize=10)
save(fig, 'fig_induction.png')
