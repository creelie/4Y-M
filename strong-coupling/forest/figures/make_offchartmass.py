"""Figure for Theorem thm:offchartmass and Corollary cor:offchartmass: (a) the
lower bound 1 - exp(-15|W| e^{-10^{-24} beta_W}) on the off-chart mass p, as a
surface over beta_W/10^24 and log10|W|; (b) the curve beta_W = 10^24 (4.40 +
ln|W|) above which the bound of Theorem thm:badcross improves on the trivial
one, with the region where the lower bound on p exceeds 1/2.
Run from this directory: python3 make_offchartmass.py. Writes fig_offchartmass.png."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK = '#3b6fb6', '#c0392b', '#2e8b57', '#222222'

fig = plt.figure(figsize=(9.0, 3.7))
ax = fig.add_axes([-0.03, 0.0, 0.55, 1.0], projection='3d')
b = np.linspace(0, 16, 121)
w = np.linspace(0, 4, 81)
B, Wl = np.meshgrid(b, w)
P = 1 - np.exp(-15 * 10 ** Wl * np.exp(-B))
ax.plot_surface(B, Wl, P, cmap='RdBu_r', vmin=-0.2, vmax=1.2, rstride=2, cstride=2,
                lw=0.1, edgecolor='#ffffff', alpha=0.95)
wl = np.linspace(0, 4, 50)
line = 4.40 + wl * np.log(10)
ax.plot(line, wl, 1 - np.exp(-15 * 10 ** wl * np.exp(-line)), color=GREEN, lw=2.0, zorder=10)
ax.set_xlabel(r'$\beta_W/10^{24}$', labelpad=2)
ax.set_ylabel(r'$\log_{10}|W|$', labelpad=2)
ax.set_zticks([0, 0.5, 1])
ax.set_zlim(0, 1)
ax.view_init(28, -58)
ax.set_box_aspect((1.2, 1, 0.75))
ax.tick_params(labelsize=8, pad=0)
ax.set_title(r'(a) lower bound on the off-chart mass $p$', fontsize=10, y=0.98)

bx = fig.add_axes([0.6, 0.17, 0.37, 0.7])
half = np.log(15 * 10 ** wl / np.log(2))
bx.fill_between(wl, 0, half, color='#f3c9c3', alpha=0.8, lw=0)
bx.plot(wl, line, color=GREEN, lw=1.6)
bx.text(0.15, 1.2, r'$p>1/2$', color=RED, fontsize=10)
bx.text(0.15, 14.4, r'$\beta_W=10^{24}(4.40+\ln|W|)$', color=GREEN, fontsize=9)
bx.set_xlabel(r'$\log_{10}|W|$')
bx.set_ylabel(r'$\beta_W/10^{24}$')
bx.set_xlim(0, 4); bx.set_ylim(0, 16)
bx.spines[['top', 'right']].set_visible(False)
bx.set_title(r'(b) where the crossing bound is useful', fontsize=10)
save(fig, 'fig_offchartmass.png')
