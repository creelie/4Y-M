"""Figure for the fibre route to the coarsest condition (ssec:bpfibre).

fig_bpfibre: (a) the fibre law of one step of the straight-path map with
blocks of side b = 4, drawn in a plane. The coarse data fix the marks (red);
the other links of the paths are fixed to I (black) and the links off the
paths (grey) are integrated. In the Gaussian model the minimiser carries a
current K Theta that vanishes off the paths, is constant along each path
with the value (Gamma v)_e and is conserved at each block corner
(Proposition prop:bpfibgauss (ii)). (b) the proved bounds
1/b^2 <= Gamma_ee <= 6/b and the lower bound 2/b^2 on the row sums gamma_b
(Propositions prop:bpfibgauss, prop:bpfibbloch), against the decay
gamma_b = O(1/b) that condition (G) asks for, drawn with an arbitrary
constant. (c) the lower bound 4 beta delta^2/(pi^2 b^2) on the action of
the background that coarse data of SU(2) force on the square of a coarse
plaquette, for delta = sqrt(2) pi, b = 2^k, along the illustrative tuning
beta_k = 2 + (11/(3 pi^2)) k log 2 (Proposition prop:bpfibrough): it tends
to zero, so it obstructs nothing.

Run from this directory: python3 make_fibre.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def tuning(k, beta0=2.0):
    return beta0 + 11 / (3 * np.pi ** 2) * k * np.log(2)


def fig_bpfibre():
    fig = plt.figure(figsize=(10.0, 3.6))
    ax = fig.add_axes([0.0, 0.05, 0.36, 0.88])
    b, nb = 4, 2
    n = b * nb
    for i in range(n + 1):
        ax.plot([0, n], [i, i], color=GREY, lw=0.5, zorder=1)
        ax.plot([i, i], [0, n], color=GREY, lw=0.5, zorder=1)
    for yx in range(nb + 1):
        for yy in range(nb + 1):
            x0, y0 = b * yx, b * yy
            if yx < nb:
                ax.plot([x0, x0 + b], [y0, y0], color=INK, lw=1.8, zorder=2)
                ax.plot([x0 + b // 2, x0 + b // 2 + 1], [y0, y0], color=RED, lw=3.4, zorder=3, solid_capstyle='butt')
            if yy < nb:
                ax.plot([x0, x0], [y0, y0 + b], color=INK, lw=1.8, zorder=2)
                ax.plot([x0, x0], [y0 + b // 2, y0 + b // 2 + 1], color=RED, lw=3.4, zorder=3, solid_capstyle='butt')
    # current through the corner (4,4)
    cx, cy = 4, 4
    arrows = [((0.6, 4), (3.4, 4)), ((4.6, 4), (7.4, 4)), ((4, 0.6), (4, 3.4)), ((4, 4.6), (4, 7.4))]
    for p, q in arrows:
        off = 0.32
        if p[1] == q[1]:
            p2, q2 = (p[0], p[1] + off), (q[0], q[1] + off)
        else:
            p2, q2 = (p[0] - off, p[1]), (q[0] - off, q[1])
        ax.add_patch(FancyArrowPatch(p2, q2, arrowstyle='-|>', mutation_scale=9, color=BLUE, lw=1.3, zorder=4))
    ax.plot([cx], [cy], 'o', color=GOLD, ms=7, zorder=5, mec=INK, mew=0.6)
    ax.text(1.0, 4.55, r'$(\Gamma v)_{e}$', color=BLUE, fontsize=9)
    ax.text(5.1, 4.55, r'$(\Gamma v)_{e^\prime}$', color=BLUE, fontsize=9)
    ax.text(4.25, 2.2, r'$\sum_\nu=0$', color=GOLD, fontsize=9.5)
    ax.text(-0.3, -0.9, 'marks (red): coarse data;  black: fixed to $I$;  grey: integrated', fontsize=7.6)
    ax.set_xlim(-0.6, n + 0.6); ax.set_ylim(-1.3, n + 0.9); ax.set_aspect('equal'); ax.axis('off')
    ax.text(-0.5, n + 0.4, '(a)', fontsize=11)

    bx = fig.add_axes([0.42, 0.17, 0.25, 0.70])
    bb = np.logspace(np.log10(3), np.log10(1024), 200)
    bx.loglog(bb, 6 / bb, color=INK, lw=1.4, label=r'$6/b\geq\Gamma_{ee}$')
    bx.loglog(bb, 1 / bb ** 2, color=INK, lw=1.4, ls='--', label=r'$1/b^2\leq\Gamma_{ee}$')
    bx.loglog(bb, 2 / bb ** 2, color=BLUE, lw=1.3, label=r'$2/b^2\leq\gamma_b$')
    bx.loglog(bb, 20 / bb, color=RED, lw=1.0, ls=':', label=r'(G): $\gamma_b=O(1/b)$')
    bx.set_xlabel(r'$b$')
    bx.set_title(r'(b)  the effective form, Gaussian model', fontsize=10, loc='left')
    bx.legend(fontsize=7.5, frameon=False, loc='lower left')

    cx_ = fig.add_axes([0.745, 0.17, 0.235, 0.70])
    ks = np.arange(1, 13)
    D2 = 2 * np.pi ** 2
    bound = 4 * tuning(ks) * D2 / (np.pi ** 2 * 4.0 ** ks)
    cx_.semilogy(ks, bound, 'o-', color=BLUE, lw=1.3, ms=3.5)
    cx_.set_xticks([2, 4, 6, 8, 10, 12]); cx_.set_xlabel(r'$k$ (blocks of side $2^k$)')
    cx_.set_title(r'(c)  forced action, $\mathrm{SU}(2)$', fontsize=10, loc='left')
    save(fig, 'fig_bpfibre.png')


if __name__ == '__main__':
    fig_bpfibre()
