"""Figure of the extension of Theorem "effectiveupper": (a) one direction of a
block of side L = 4, the block centres c_z as nodes, the hat functions and
the straight path of length L from a site, which meets three segments
between centres; the averages over the block are the weights 1/8, 3/4, 1/8;
(b) the symbol 3/4 + cos(theta)/4 of this filter, the same at every scale,
against the modulus of the filter obtained with the coarse sites as nodes,
which tends to zero at theta = pi. Run from this directory:
python3 make_whitney.py. Writes fig_whitney.png."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
GREY, BLUE, RED, ORANGE, GREEN, INK = '#c9c9c9', '#3b6fb6', '#c0392b', '#e67e22', '#2e8b57', '#222222'


def blocks(ax):
    L = 4
    xs = np.arange(-6, 10)
    for b in range(-2, 3):
        ax.axvspan(L * b - 0.5, L * b + L - 0.5, color=GREY if b % 2 else '#ececec', alpha=0.6, lw=0)
    ax.plot(xs, np.zeros_like(xs), 'o', color=INK, ms=3)
    t = np.linspace(-8, 12, 801)
    cols = {-1: GREEN, 0: BLUE, 1: ORANGE}
    for z in (-1, 0, 1):
        c = L * z + (L - 1) / 2
        ax.plot(t, 1.1 * np.maximum(0, 1 - np.abs(t - c) / L), color=cols[z], lw=1.6)
        ax.plot([c], [0], marker='D', color=cols[z], ms=6)
        ax.text(c, -0.2, r'$c_{%d}$' % z if z >= 0 else r'$c_{-1}$', ha='center', va='top', fontsize=10, color=cols[z])
    s = 2
    ax.annotate('', xy=(s + L, 1.35), xytext=(s, 1.35),
                arrowprops=dict(arrowstyle='->', color=RED, lw=1.6))
    ax.text(s + L / 2, 1.42, r'path of length $L$ from $x$', ha='center', va='bottom', color=RED, fontsize=9.5)
    ax.plot([s, s], [0, 1.35], ls=':', color=RED, lw=1)
    for z, wgt in zip((-1, 0, 1), (r'$\frac{1}{8}$', r'$\frac{3}{4}$', r'$\frac{1}{8}$')):
        ax.text(L * z + (L - 1) / 2, 0.55 if z == 0 else 0.35, wgt, ha='center', fontsize=12, color=cols[z],
                bbox=dict(fc='white', ec='none', alpha=0.8, pad=1))
    ax.text(1.5, -0.55, r'block $B_k(0)$, $L=4$', ha='center', fontsize=9.5)
    ax.set_xlim(-6.5, 9.5)
    ax.set_ylim(-0.75, 1.75)
    ax.axis('off')
    ax.set_title('(a) nodes at the block centres', fontsize=10.5)


def symbol(ax):
    th = np.linspace(0, np.pi, 400)
    ax.plot(th, 0.75 + 0.25 * np.cos(th), color=BLUE, lw=2, label=r'centres: $\frac{3}{4}+\frac{1}{4}\cos\theta$, every $L$')
    for L, ls in ((4, '--'), (16, ':')):
        a = (L + 1) / (2 * L)
        ax.plot(th, np.abs(a + (1 - a) * np.exp(1j * th)), color=RED, ls=ls, lw=1.6,
                label=r'coarse sites, $L=%d$' % L)
    ax.axhline(0.5, color=GREY, lw=0.8)
    ax.set_xlim(0, np.pi)
    ax.set_ylim(0, 1.05)
    ax.set_xticks([0, np.pi / 2, np.pi])
    ax.set_xticklabels(['0', r'$\pi/2$', r'$\pi$'])
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel('modulus of the filter symbol')
    ax.legend(loc='lower left', fontsize=8.5, frameon=False)
    ax.set_title('(b) the filter in one direction', fontsize=10.5)


fig, axs = plt.subplots(1, 2, figsize=(7.6, 2.9), gridspec_kw={'width_ratios': [1.45, 1]})
blocks(axs[0])
symbol(axs[1])
fig.tight_layout()
save(fig, 'fig_whitney.png')
