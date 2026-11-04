"""Figure for the smoothing along the paths (ssec:bppaths).

fig_bppaths: (a) one step of the straight-path map with blocks of side
b = 6, drawn in a plane. The strict interiors of the paths (red) lie in no
common plaquette, so given the other fine links (the end links of the paths
in black and the links off the paths in grey) they are independent one-link
Wilson laws, and each coarse link is a product of b - 2 of them between two
fixed end links (Lemma lem:bppathint). (b) the multipliers
a_d(B) = I_d(B)/I_1(B) of the one-link law on the irreducible representations
of dimension d = 2, ..., 6, all below the bound 2B/(2B+3)
(Lemma lem:bppathmodes). (c) the bounds eps, K^2 and h of
Theorem thm:bppathmixed for blocks of side 2^k along the illustrative tuning
beta_W,k = 2 + (11/(3 pi^2)) k log 2, on a logarithmic scale; the curves are
computed from their closed forms through log10, and h, K^2
are drawn where eps < 1/2.

Run from this directory: python3 make_paths.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import ive
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def tuning(k, beta0=2.0):
    return beta0 + 11 / (3 * np.pi ** 2) * k * np.log(2)


def log10_bounds(k):
    beta = tuning(k)
    b = 2.0 ** k
    B = 6 * beta
    q = 4 * beta / (4 * beta + 1)
    le = np.log10(48.0) + 2 * np.log10(1 + B) + (b - 4) * np.log10(q)
    eps = 10.0 ** np.maximum(le, -300)
    # valid where eps < 1/2: K = B eps/(1-eps), h = (B^2+B) eps/(1-eps) + K^2
    r = 1 - np.minimum(eps, 0.5)
    lK2 = 2 * (np.log10(B) + le - np.log10(r))
    lh = np.log10((B ** 2 + B) * 10.0 ** np.maximum(le, -300) / r + 10.0 ** np.maximum(lK2, -300))
    lh = np.where(le > -150, lh, np.log10(B ** 2 + B) + le - np.log10(r))
    return le, lK2, lh


def fig_bppaths():
    fig = plt.figure(figsize=(10.0, 3.6))
    ax = fig.add_axes([0.0, 0.05, 0.34, 0.88])
    b, nb = 6, 2
    n = b * nb
    for i in range(n + 1):
        ax.plot([0, n], [i, i], color=GREY, lw=0.5, zorder=1)
        ax.plot([i, i], [0, n], color=GREY, lw=0.5, zorder=1)
    for yx in range(nb + 1):
        for yy in range(nb + 1):
            x0, y0 = b * yx, b * yy
            if yx < nb:
                ax.plot([x0, x0 + 1], [y0, y0], color=INK, lw=2.0, zorder=2, solid_capstyle='butt')
                ax.plot([x0 + b - 1, x0 + b], [y0, y0], color=INK, lw=2.0, zorder=2, solid_capstyle='butt')
                ax.plot([x0 + 1, x0 + b - 1], [y0, y0], color=RED, lw=3.0, zorder=3, solid_capstyle='butt')
            if yy < nb:
                ax.plot([x0, x0], [y0, y0 + 1], color=INK, lw=2.0, zorder=2, solid_capstyle='butt')
                ax.plot([x0, x0], [y0 + b - 1, y0 + b], color=INK, lw=2.0, zorder=2, solid_capstyle='butt')
                ax.plot([x0, x0], [y0 + 1, y0 + b - 1], color=RED, lw=3.0, zorder=3, solid_capstyle='butt')
    for yx in range(nb + 1):
        for yy in range(nb + 1):
            ax.plot([b * yx], [b * yy], 'o', color=GOLD, ms=4.5, zorder=4, mec=INK, mew=0.5)
    ax.text(-0.3, -0.9, r'red: interiors $Y_1,\ldots,Y_{b-2}$;  black: end links;  grey: off the paths', fontsize=7.4)
    ax.set_xlim(-0.6, n + 0.6); ax.set_ylim(-1.3, n + 0.9); ax.set_aspect('equal'); ax.axis('off')
    ax.text(-0.5, n + 0.4, '(a)', fontsize=11)

    bx = fig.add_axes([0.41, 0.17, 0.25, 0.70])
    B = np.linspace(0.01, 40, 400)
    for d, ls in zip(range(2, 7), ['-', '--', '-.', ':', (0, (5, 1))]):
        bx.plot(B, ive(d, B) / ive(1, B), color=BLUE, lw=1.2, ls=ls, label=r'$a_{%d}$' % d)
    bx.plot(B, 2 * B / (2 * B + 3), color=RED, lw=1.0, label=r'$2B/(2B+3)$')
    bx.set_xlabel(r'$B$'); bx.set_ylim(0, 1.02)
    bx.set_title(r'(b)  multipliers of $f_B$', fontsize=10, loc='left')
    bx.legend(fontsize=7.5, frameon=False, loc='lower right', ncol=2)

    cx = fig.add_axes([0.735, 0.17, 0.245, 0.70])
    ks = np.arange(6, 15)
    le, lK2, lh = log10_bounds(ks)
    cx.plot(ks, le, 'o-', color=BLUE, lw=1.2, ms=3.2, label=r'$\log_{10}\epsilon$')
    ok = le < np.log10(0.5)
    cx.plot(ks[ok], lh[ok], 's--', color=GREEN, lw=1.0, ms=3.0, label=r'$\log_{10}h$')
    cx.plot(ks[ok], lK2[ok], '^:', color=RED, lw=1.0, ms=3.0, label=r'$\log_{10}K^2$')
    cx.axhline(0, color=GREY, lw=0.7)
    cx.set_yscale('symlog', linthresh=10)
    cx.set_xlabel(r'$k$ (blocks of side $2^k$)')
    cx.set_title(r'(c)  along the tuning', fontsize=10, loc='left')
    cx.legend(fontsize=7.5, frameon=False, loc='lower left')
    save(fig, 'fig_bppaths.png')


if __name__ == '__main__':
    fig_bppaths()
