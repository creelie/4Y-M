"""Figures for the split inequality at strong coupling and for the
coarsest-scale condition.

fig_bpsplitstrong: (a) the exterior information of one fine cube. Changing a
coarse link e moves the fibre law through the marked link l(e) inside the
cube, by at most q per partner sum, and the cube, resampled given its
exterior, feels the change through its boundary links, again by at most q;
two cubes at distance s see each other through fibre correlations that
decay like xi^(s-n'-1)/(1-xi) (Theorems thm:bpexterior, thm:bpsplitstrong).
(b) the constants C_A and C_B (bn/n')^4 of Theorem thm:bpsplitstrong for
SU(2) against |beta_W|, at n' = 8 and in the limit n' -> infinity, with the
refinement value 1 (Corollary cor:bpsplitranges).

fig_bphidden: (a) the law of Proposition prop:bphidden(ii): one sign S,
drawn once for the torus, chooses between two product tilts of the selected
plaquettes. (b) the law of the spatial average F under that measure, drawn
with normal approximations of its two components (an illustration only),
for L = 8, 16, 32, epsilon = 1/2, SU(2): the two modes at +-epsilon/4 do
not merge as L grows, while resampling one cube moves F by at most
2 r_x / M.

Run from this directory: python3 make_split.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def constants(beta, n):
    xi = 4.5 * beta
    c = beta / 4
    q = xi / (1 - xi)
    if n is None:
        alpha, rn = 0.0, 16.0
    else:
        alpha = c * (18 * n - 3) / ((1 - 18 * c) * n * n)
        rn = ((2 * n + 3) ** 4 + 10 * (n + 2) ** 3 * xi * (8 - 5 * xi + 4 * xi ** 2 - xi ** 3) / (1 - xi) ** 5) / n ** 4
    eta = q * np.sqrt(rn / (1 - alpha))
    return 1 / ((1 - eta) * (1 - alpha)), 1 / (1 - q * q - eta)


def fig_bpsplitstrong():
    fig = plt.figure(figsize=(9.6, 4.1))
    ax = fig.add_axes([0.0, 0.02, 0.56, 0.93])
    ax.set_xlim(-0.6, 26.6); ax.set_ylim(-1.6, 13.2); ax.set_aspect('equal'); ax.axis('off')
    W, H, b = 26, 12, 4
    for x in range(W + 1):
        ax.plot([x, x], [0, H], color=GREY, lw=0.3, alpha=0.6)
    for y in range(H + 1):
        ax.plot([0, W], [y, y], color=GREY, lw=0.3, alpha=0.6)
    for x in range(0, W + 1, b):
        ax.plot([x, x], [0, H], color=INK, lw=1.1, alpha=0.55)
    for y in range(0, H + 1, b):
        ax.plot([0, W], [y, y], color=INK, lw=1.1, alpha=0.55)
    cubes = [(2, 3, 6), (17, 3, 6)]
    for (x0, y0, s) in cubes:
        ax.add_patch(Rectangle((x0, y0), s, s, fc=RED, ec=RED, alpha=0.10, lw=0))
        ax.add_patch(Rectangle((x0, y0), s, s, fc='none', ec=RED, lw=1.6))
        ax.add_patch(Rectangle((x0 - 1, y0 - 1), s + 2, s + 2, fc='none', ec=RED, lw=0.9, ls=(0, (3, 2))))
    # coarse link e: path along y = 4 from x = 4 to 8, marked link inside the cube
    ax.plot([4, 8], [4, 4], color=BLUE, lw=3.2, solid_capstyle='butt')
    ax.plot([6, 7], [4, 4], color=GOLD, lw=4.2, solid_capstyle='butt')
    ax.text(4.2, 3.25, r'$\pi(e)$', color=BLUE, fontsize=9)
    ax.text(6.0, 4.35, r'$l(e)$', color='#9a7400', fontsize=9)
    ax.text(2.1, 9.35, r"$V'_{x'}$", color=RED, fontsize=10)
    ax.text(17.1, 9.35, r"$V'_{x''}$", color=RED, fontsize=10)
    ax.text(0.9, 1.1, r"$\partial V'_{x'}$", color=RED, fontsize=8.5)
    a1 = FancyArrowPatch((6.5, 4.2), (8.6, 7.6), arrowstyle='-|>', mutation_scale=10, color=GREEN, lw=1.2,
                         connectionstyle='arc3,rad=-0.3')
    ax.add_patch(a1)
    ax.text(8.9, 7.7, r'fibre: $q$', color=GREEN, fontsize=9)
    a2 = FancyArrowPatch((8.6, 7.1), (5.2, 6.4), arrowstyle='-|>', mutation_scale=10, color=GREEN, lw=1.2,
                         connectionstyle='arc3,rad=-0.3')
    ax.add_patch(a2)
    ax.text(3.6, 10.6, r'exterior information $\leq q^2$', color=GREEN, fontsize=9)
    ax.annotate('', xy=(16.0, 6.0), xytext=(9.0, 6.0),
                arrowprops=dict(arrowstyle='<|-|>', color=INK, lw=1.0, ls=(0, (4, 2))))
    ax.text(12.5, 6.35, r"$\xi^{s-n'-1}/(1-\xi)$", ha='center', fontsize=9, color=INK)
    ax.text(12.5, 5.2, 'fibre decay', ha='center', fontsize=8.5, color='#555555')
    ax.text(0.0, -1.2, r'thin: fine links;  thick: paths of the coarse links ($b=4$)', fontsize=8.5, color='#444444')
    ax.text(-0.5, 12.6, '(a)', fontsize=11)

    bx = fig.add_axes([0.64, 0.15, 0.34, 0.75])
    beta = np.linspace(0, 0.012, 200)
    for n, ls, lab in ((8, '-', r"$n'=8$"), (None, '--', r"$n'\to\infty$")):
        ca, cb = np.array([constants(t, n) for t in beta]).T
        bx.plot(beta, ca, ls, color=RED, lw=1.5, label=r'$C_A$, ' + lab)
        bx.plot(beta, cb, ls, color=BLUE, lw=1.5, label=r"$C_B(bn/n')^4$, " + lab)
    bx.axhline(1, color=GREEN, lw=1.0)
    bx.text(0.0062, 1.015, 'Haar refinements', color=GREEN, fontsize=8.5)
    bx.axvline(0.01, color=GREY, lw=0.8, ls=':')
    bx.set_xlabel(r'$|\beta_W|$'); bx.set_xlim(0, 0.012); bx.set_ylim(0.95, 1.75)
    bx.legend(fontsize=7.6, frameon=False, loc='upper left')
    bx.set_title('(b)  SU(2)', fontsize=10, loc='left')
    save(fig, 'fig_bpsplitstrong.png')


def fig_bphidden():
    fig = plt.figure(figsize=(9.6, 4.2))
    ax = fig.add_axes([0.0, 0.0, 0.52, 1.0], projection='3d')
    clean3d(ax)
    ax.set_box_aspect((1, 1, 0.8))
    L = 8
    rng = np.random.default_rng(3)
    for z, col, lab in ((2.6, RED, r'$S=+1$:  $\prod_i(1+\epsilon X_i)$'), (0.0, BLUE, r'$S=-1$:  $\prod_i(1-\epsilon X_i)$')):
        for x in range(L + 1):
            ax.plot([x, x], [0, L], [z, z], color=GREY, lw=0.35)
            ax.plot([0, L], [x, x], [z, z], color=GREY, lw=0.35)
        X, Y = np.meshgrid([0, L], [0, L])
        ax.plot_surface(X, Y, np.full_like(X, z, dtype=float), color=col, alpha=0.05, linewidth=0, shade=False)
        for i in range(0, L, 2):
            for j in range(0, L, 2):
                a = 0.25 + 0.55 * rng.random()
                ax.add_collection3d(matplotlib.collections.PolyCollection(
                    [[(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)]], facecolors=col, alpha=a, edgecolors=col),
                    zs=z, zdir='z')
        ax.text(L + 0.4, 0, z + 0.15, lab, fontsize=9, color=col)
    ax.text(-1.4, -1.0, 1.25, r'one sign $S$', fontsize=9.5, color=INK)
    ax.quiver(-0.5, -0.5, 1.5, 0, 0, 0.9, color=INK, lw=1.2, arrow_length_ratio=0.25)
    ax.quiver(-0.5, -0.5, 1.1, 0, 0, -0.9, color=INK, lw=1.2, arrow_length_ratio=0.25)
    ax.view_init(elev=22, azim=-62)
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, 3.0)
    ax.text2D(0.02, 0.93, '(a)', transform=ax.transAxes, fontsize=11)

    bx = fig.add_axes([0.60, 0.20, 0.37, 0.68])
    eps, s2 = 0.5, 0.25
    m, v = eps * s2, s2 - (eps * s2) ** 2
    x = np.linspace(-0.3, 0.3, 800)
    for Lt, col in ((8, GREY), (16, BLUE), (32, INK)):
        M = (Lt / 2) ** 4
        sd = np.sqrt(v / M)
        dens = 0.5 * (np.exp(-(x - m) ** 2 / (2 * sd * sd)) + np.exp(-(x + m) ** 2 / (2 * sd * sd)))
        bx.plot(x, dens / dens.max(), color=col, lw=1.3, label=r'$L=%d$' % Lt)
    bx.axvline(m, color=RED, lw=0.8, ls=':'); bx.axvline(-m, color=BLUE, lw=0.8, ls=':')
    bx.set_xlabel(r'$F=M^{-1}\sum_i X_i$'); bx.set_yticks([]); bx.set_ylabel('density, unit peak', fontsize=9); bx.set_ylim(0, 1.35)
    bx.legend(fontsize=8.5, frameon=False, loc='upper center')
    bx.set_title(r'(b)  law of $F$, $\mathrm{SU}(2)$, $\epsilon=1/2$', fontsize=10, loc='left')
    fig.text(0.60, 0.015, r'two modes at $\pm\epsilon s_2^2$; one cube moves $F$ by at most $2r_x/M$',
             fontsize=8.3, color='#444444')
    save(fig, 'fig_bphidden.png')


if __name__ == '__main__':
    fig_bpsplitstrong()
    fig_bphidden()
