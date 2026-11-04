"""Figure for the single-link coefficients across scales (ssec:bpbalaban).

fig_bpbalaban: (a) one step of the straight-path map with blocks of side
b = 4, drawn in a plane. Each coarse link is the ordered product of the b
fine links of its path, marked at the link with base point b y + floor(b/2)
e_nu. In the axial gauge the law of the block variables is the marginal on
the marks; two marks are at distance at least floor(b/2), and the influence
of one mark on another runs through walks of fine links, each step weighted
by a Dobrushin coefficient (Lemma lem:bpbalmarginal, Theorem thm:bpbalstep).
(b) the lower bound 2 - 2 sqrt(3/B) - (3 pi^3/4) e^(1/2) B^(3/2) e^(-B),
B = 6 beta_W, on the Dobrushin sum of the SU(2) Wilson measure that gauge
covariance forces (Corollary cor:bpbalwilson), against beta_W, with the
threshold 1 of the Dobrushin condition. (c) the lower bound
2 - 4 (1 - eps_b)/pi, eps_b = (1 + 1/(4 beta_W))^(-b), at the coarsest scale,
b = 2^k, along the illustrative tuning beta_W,k = 2 + (11/(3 pi^2)) k log 2
(Proposition prop:bpbalrandom): above 1 the coarse law cannot satisfy the
Dobrushin condition; below 1 gauge covariance no longer excludes it, and
the bound tends to 2 - 4/pi.

Run from this directory: python3 make_balaban.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def wilson_bound(bw):
    B = 6 * bw
    return 2 - 2 * np.sqrt(3 / B) - 0.75 * np.pi ** 3 * np.exp(0.5) * B ** 1.5 * np.exp(-B)


def coarse_bound(k, beta0=2.0):
    beta = beta0 + 11 / (3 * np.pi ** 2) * k * np.log(2)
    eps = (1 + 1 / (4 * beta)) ** (-(2.0 ** k))
    return 2 - 4 * (1 - eps) / np.pi


def fig_bpbalaban():
    fig = plt.figure(figsize=(10.0, 3.5))
    ax = fig.add_axes([0.01, 0.06, 0.34, 0.86])
    b, nb = 4, 3
    n = b * nb
    for i in range(n + 1):
        ax.plot([0, n], [i, i], color=GREY, lw=0.4, zorder=1)
        ax.plot([i, i], [0, n], color=GREY, lw=0.4, zorder=1)
    marks = []
    for yx in range(nb):
        for yy in range(nb):
            x0, y0 = b * yx, b * yy
            ax.plot([x0, x0 + b], [y0, y0], color=INK, lw=1.6, zorder=2)
            ax.plot([x0, x0], [y0, y0 + b], color=INK, lw=1.6, zorder=2)
            mx = (x0 + b // 2, y0)
            my = (x0, y0 + b // 2)
            ax.plot([mx[0], mx[0] + 1], [mx[1], mx[1]], color=RED, lw=3.2, zorder=3, solid_capstyle='butt')
            ax.plot([my[0], my[0]], [my[1], my[1] + 1], color=RED, lw=3.2, zorder=3, solid_capstyle='butt')
            marks += [mx, my]
    ax.plot([0, n], [n, n], color=INK, lw=1.6); ax.plot([n, n], [0, n], color=INK, lw=1.6)
    # a walk from mark e to mark j through unmarked fine links
    walk = [(6.5, 4.0), (7.0, 4.5), (7.0, 5.5), (7.5, 6.0), (8.0, 6.5)]
    for (p, q) in zip(walk[:-1], walk[1:]):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=8, color=BLUE, lw=1.1, zorder=4))
    ax.text(6.3, 3.35, r'$l(e)$', color=RED, fontsize=9)
    ax.text(8.25, 6.55, r'$l(j)$', color=RED, fontsize=9)
    ax.text(8.6, 4.55, r'$C_{ek_1}C_{k_1k_2}\cdots$', color=BLUE, fontsize=8.5)
    ax.annotate('', xy=(4.0, -0.55), xytext=(8.0, -0.55), arrowprops=dict(arrowstyle='<->', color=INK, lw=0.8))
    ax.text(5.3, -1.25, r'side $b$', fontsize=8.5)
    ax.set_xlim(-0.6, n + 0.8); ax.set_ylim(-1.6, n + 0.9); ax.set_aspect('equal'); ax.axis('off')
    ax.text(-0.5, n + 0.4, '(a)', fontsize=11)
    ax.text(0.0, -1.25, 'marks', color=RED, fontsize=8.5)

    bx = fig.add_axes([0.42, 0.17, 0.25, 0.70])
    bw = np.linspace(0.6, 12, 400)
    bx.plot(bw, wilson_bound(bw), color=INK, lw=1.5, label='lower bound on $c_\\mu$')
    bx.axhline(1, color=RED, lw=0.9, ls='--', label='Dobrushin threshold')
    bx.axvline(8 / 3, color=GREY, lw=0.8, ls=':')
    bx.text(8 / 3 + 0.25, 1.85, r'$\beta_W=\frac{8}{3}$', fontsize=8.5, color='#444444')
    bx.set_xlabel(r'$\beta_W$'); bx.set_ylim(0, 2.05)
    bx.set_title(r'(b)  finest level, $\mathrm{SU}(2)$ Wilson', fontsize=10, loc='left')
    bx.legend(fontsize=8, frameon=False, loc='lower right', bbox_to_anchor=(1.0, 0.12))

    cx = fig.add_axes([0.73, 0.17, 0.25, 0.70])
    ks = np.arange(1, 13)
    cx.plot(ks, coarse_bound(ks), 'o-', color=BLUE, lw=1.3, ms=3.5, label='lower bound on $c$')
    cx.axhline(1, color=RED, lw=0.9, ls='--')
    cx.axhline(2 - 4 / np.pi, color=GREY, lw=0.8, ls=':')
    cx.text(7.2, 2 - 4 / np.pi - 0.17, r'$2-4/\pi$', fontsize=8.5, color='#444444')
    cx.set_xticks([2, 4, 6, 8, 10, 12]); cx.set_xlabel(r'$k$ (blocks of side $2^k$)'); cx.set_ylim(0, 2.05)
    cx.set_title(r'(c)  coarsest scale along the tuning', fontsize=10, loc='left')
    cx.legend(fontsize=8, frameon=False, loc='upper right')
    save(fig, 'fig_bpbalaban.png')


if __name__ == '__main__':
    fig_bpbalaban()
