"""Figure for Proposition prop:scalelosses and Corollary cor:smearedtrial:
(a) the product of the factors 1 - c_j over the scales between the lattice
spacing and the separation sigma along the one-loop SU(2) tuning, for losses
of first order (c_j = c g_j^2) and of second order (c_j = c g_j^4) in the
running coupling, as a surface over the number of scales J and the
constant c; (b) the section at c = 0.1, with beta_W = 20 at the scale sigma, with the rate (log 1/a)^(-4c/h).
Run from this directory: python3 make_scalelosses.py. Writes
fig_scalelosses.png (flattened, 256 colours)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, INK = '#3b6fb6', '#c0392b', '#222222'

b = 11 / (3 * np.pi ** 2)
h = b * np.log(2)
beta_sigma = 20.0         # beta_W at the scale sigma


def products(J, c):
    j = np.arange(J + 1)
    beta = beta_sigma + h * (J - j)
    g2 = 4 / beta
    first = np.prod(1 - np.minimum(c * g2, 0.5))
    second = np.prod(1 - np.minimum(c * g2 ** 2, 0.5))
    return first, second


Js = np.arange(0, 2001, 50)
cs = np.linspace(0.02, 0.2, 19)
F = np.zeros((len(cs), len(Js)))
S = np.zeros_like(F)
for a, c in enumerate(cs):
    for k, J in enumerate(Js):
        F[a, k], S[a, k] = products(J, c)

fig = plt.figure(figsize=(10, 4.2))
ax = fig.add_subplot(1, 2, 1, projection='3d')
X, Y = np.meshgrid(Js, cs)
ax.plot_surface(X, Y, S, color=BLUE, alpha=0.55, edgecolor=BLUE, linewidth=0.25)
ax.plot_surface(X, Y, F, color=RED, alpha=0.75, edgecolor='#7b241c', linewidth=0.25)
ax.set_xlabel(r'scales $J=\log_2(\sigma/a_k)$', labelpad=4)
ax.set_ylabel(r'$c$', labelpad=2)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r'$\prod_j(1-c_j)$', labelpad=6, rotation=90)
ax.set_zlim(0, 1)
ax.view_init(elev=22, azim=-58)
ax.set_box_aspect((1.3, 1, 0.8))
ax.text2D(0.66, 0.74, r'$c_j=c\,g_j^4$', color=BLUE, fontsize=11, transform=ax.transAxes)
ax.text2D(0.66, 0.42, r'$c_j=c\,g_j^2$', color=RED, fontsize=11, transform=ax.transAxes)
ax.set_title('(a)', loc='left')

ax2 = fig.add_subplot(1, 2, 2)
c = 0.1
Jfine = np.arange(0, 2001)
f = np.array([products(J, c)[0] for J in Jfine])
s = np.array([products(J, c)[1] for J in Jfine])
ax2.semilogy(Jfine, f, color=RED, lw=1.6, label=r'first order, $c_j=0.1\,g_j^2$')
ax2.semilogy(Jfine, s, color=BLUE, lw=1.6, label=r'second order, $c_j=0.1\,g_j^4$')
rate = (beta_sigma / (beta_sigma + h * (Jfine + 1))) ** (4 * c / h)
ax2.semilogy(Jfine, rate, color=INK, lw=0.9, ls='--', label=r'$(\beta_\sigma/(\beta_\sigma+b\log(2\sigma/a_k)))^{4c/h}$')
ax2.set_xlabel(r'scales $J=\log_2(\sigma/a_k)$')
ax2.set_ylabel('product over the scales')
ax2.set_ylim(1e-3, 1.5)
ax2.legend(frameon=False, fontsize=8.5, loc='lower left')
ax2.spines[['top', 'right']].set_visible(False)
ax2.set_title('(b)', loc='left')
fig.tight_layout()
save(fig, 'fig_scalelosses.png')
