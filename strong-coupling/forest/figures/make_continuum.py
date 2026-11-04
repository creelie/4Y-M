"""Figures for the averaged condition across scales and for the continuum
limit at separated points.

fig_bptower: (a) three lattices of spacings a_0, a_0/b, a_0/b^2 drawn as
stacked planes, each with a cube of the same physical side ell; the
averaged condition descends to block variables for every gauge-invariant
measure (Theorem thm:bpdescent, green arrows) and climbs one scale under the
split inequality S (Theorem thm:bpup, gold arrows). (b) the measures on the
fine lattice: Haar refinements satisfy S (Theorem thm:bprefine), the Wilson
measures at beta != 0 are not Haar refinements (Proposition
prop:bpnotrefine), so for them S is an input.

fig_limitlogic: the conditions on the lattice fields and what they give in
the limit, in three layers: inputs on reflected norms of packets (bottom),
lattice bounds (middle), properties of the Schwinger functions (top).
Solid arrows are implications proved in the text, crossed arrows are
non-implications proved by examples.

Run from this directory: python3 make_continuum.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon, FancyArrowPatch
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def _plane(ax, z, n, side, colour):
    h = side / n
    xs = np.linspace(0, side, n + 1)
    for x in xs:
        ax.plot([x, x], [0, side], [z, z], color=GREY, lw=0.35, alpha=0.8)
        ax.plot([0, side], [x, x], [z, z], color=GREY, lw=0.35, alpha=0.8)
    X, Y = np.meshgrid([0, side], [0, side])
    ax.plot_surface(X, Y, np.full_like(X, z, dtype=float), color=colour, alpha=0.10, linewidth=0, shade=False)
    return h


def _cube(ax, x0, y0, z0, s, hgt, colour):
    v = np.array([[x0, y0], [x0 + s, y0], [x0 + s, y0 + s], [x0, y0 + s], [x0, y0]])
    ax.plot(v[:, 0], v[:, 1], [z0] * 5, color=colour, lw=1.6)
    ax.plot(v[:, 0], v[:, 1], [z0 + hgt] * 5, color=colour, lw=1.0, alpha=0.8)
    for p in v[:4]:
        ax.plot([p[0], p[0]], [p[1], p[1]], [z0, z0 + hgt], color=colour, lw=1.0, alpha=0.8)


def fig_bptower():
    fig = plt.figure(figsize=(9.6, 4.4))
    ax = fig.add_axes([0.0, 0.0, 0.56, 1.0], projection='3d')
    clean3d(ax)
    ax.set_box_aspect((1, 1, 1.05))
    side, ell = 8.0, 3.0
    levels = [(0.0, 8, r'$a_0$'), (2.2, 16, r'$a_0/b$'), (4.4, 32, r'$a_0/b^2$')]
    for z, n, lab in levels:
        h = _plane(ax, z, n, side, BLUE)
        _cube(ax, 2.5, 2.5, z, ell, 0.35, RED)
        ax.text(side + 0.6, 0, z, lab, fontsize=10, color=INK)
        ax.text(2.5 + ell + 0.3, 2.5, z + 0.05, r'$n=%d$' % int(round(ell / h)), fontsize=9, color=RED)
    for (z0, _, _), (z1, _, _) in zip(levels[:-1], levels[1:]):
        ax.quiver(0.4, 0.4, z1 - 0.15, 0, 0, -(z1 - z0) + 0.45, color=GREEN, lw=1.8, arrow_length_ratio=0.18)
        ax.quiver(7.2, 7.2, z0 + 0.25, 0, 0, (z1 - z0) - 0.45, color=GOLD, lw=1.8, arrow_length_ratio=0.18)
    fig.text(0.01, 0.47, 'descent\n(every\nmeasure)', color=GREEN, fontsize=8.5)
    fig.text(0.445, 0.50, 'upward step\n(under S)', color=GOLD, fontsize=8.5)
    ax.text2D(0.03, 0.93, '(a)', transform=ax.transAxes, fontsize=11)
    ax.view_init(elev=17, azim=-58)
    ax.set_xlim(0, side); ax.set_ylim(0, side); ax.set_zlim(0, 5.2)

    bx = fig.add_axes([0.58, 0.08, 0.40, 0.84])
    bx.set_xlim(0, 10); bx.set_ylim(0, 10); bx.axis('off')
    bx.text(-0.2, 9.75, '(b)', fontsize=11)
    outer = FancyBboxPatch((0.4, 0.6), 9.2, 8.4, boxstyle='round,pad=0.1', fc='#f4f6f9', ec=GREY, lw=1.0)
    bx.add_patch(outer)
    bx.text(5.0, 8.5, 'gauge-invariant measures on the fine lattice', ha='center', fontsize=8.5, color=INK)
    ref = FancyBboxPatch((1.0, 4.3), 4.2, 3.4, boxstyle='round,pad=0.1', fc='#e6f2ea', ec=GREEN, lw=1.4)
    bx.add_patch(ref)
    bx.text(3.1, 6.7, 'Haar refinements', ha='center', fontsize=9, color=GREEN)
    bx.text(3.1, 5.6, r'S holds', ha='center', fontsize=9, color=GREEN)
    bx.text(3.1, 4.8, 'at every scale', ha='center', fontsize=8, color=GREEN)
    wil = FancyBboxPatch((5.6, 1.4), 3.6, 3.0, boxstyle='round,pad=0.1', fc='#fbeaea', ec=RED, lw=1.4)
    bx.add_patch(wil)
    bx.text(7.4, 3.5, r'Wilson, $\beta\neq0$', ha='center', fontsize=9, color=RED)
    bx.text(7.4, 2.5, 'not a refinement;', ha='center', fontsize=8, color=RED)
    bx.text(7.4, 1.9, 'S open', ha='center', fontsize=8, color=RED)
    bx.plot([1.6], [2.2], 'o', color=BLUE, ms=6)
    bx.text(2.0, 2.05, r'Haar, $\beta=0$', fontsize=8.5, color=BLUE)
    bx.annotate('', xy=(2.0, 4.25), xytext=(1.7, 2.5), arrowprops=dict(arrowstyle='->', color=BLUE, lw=1.0))
    save(fig, 'fig_bptower.png')


def _node(ax, x, y, text, colour, w=1.7, h=0.62, fs=9):
    p = FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.06', fc='white', ec=colour, lw=1.3, zorder=3)
    ax.add_patch(p)
    ax.text(x, y, text, ha='center', va='center', fontsize=fs, color=colour, zorder=4)


def _arrow(ax, a, b, colour=INK, cross=False, rad=0.0):
    ar = FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=10, color=colour, lw=1.0,
                         connectionstyle='arc3,rad=%g' % rad, shrinkA=14, shrinkB=14, zorder=2,
                         linestyle='--' if cross else '-')
    ax.add_patch(ar)
    if cross:
        m = (np.array(a) + np.array(b)) / 2
        ax.plot([m[0] - 0.13, m[0] + 0.13], [m[1] - 0.13, m[1] + 0.13], color=RED, lw=1.6, zorder=5)
        ax.plot([m[0] - 0.13, m[0] + 0.13], [m[1] + 0.13, m[1] - 0.13], color=RED, lw=1.6, zorder=5)


def fig_limitlogic():
    fig, ax = plt.subplots(figsize=(9.6, 5.6))
    ax.set_xlim(0, 13.2); ax.set_ylim(-0.45, 8.4); ax.axis('off')
    layers = [(0.35, 'packets under the Wilson measure', '#fbeaea'),
              (3.05, 'lattice bounds', '#eef3fa'),
              (5.75, 'Schwinger functions', '#e6f2ea')]
    for y0, lab, col in layers:
        ax.add_patch(Polygon([[0.5, y0], [12.2, y0], [12.9, y0 + 1.75], [1.2, y0 + 1.75]], closed=True,
                             fc=col, ec=GREY, lw=0.8, zorder=0))
        if y0 > 5:
            ax.text(12.1, y0 + 0.1, lab, fontsize=8.5, color='#555555', style='italic', ha='right')
        else:
            ax.text(0.75, y0 + 0.1, lab, fontsize=8.5, color='#555555', style='italic')
    yb, ym, yt = 1.15, 3.85, 6.78
    # inputs
    _node(ax, 1.9, yb, '(SD), (NE)', RED)
    _node(ax, 6.0, yb, r'(SD$_{\rm p}$)', RED, w=1.2)
    _node(ax, 8.4, yb, r'(SB$_n$), (RI$_n$)', RED, w=1.9)
    _node(ax, 10.7, yb, r'(MC), (O$_{\rm s}$)', RED, w=1.8)
    # lattice bounds
    _node(ax, 2.2, ym, r'(M$_{\rm s}$)', BLUE, w=1.1)
    _node(ax, 4.0, ym, r'(RN$_m$)', BLUE, w=1.2)
    _node(ax, 6.0, ym, r'(G$_{\rm s}$)', BLUE, w=1.1)
    _node(ax, 8.4, ym, r'(W$_{\rm s}$)', BLUE, w=1.1)
    _node(ax, 10.7, ym, r'(U$_{\rm s}$)', BLUE, w=1.1)
    # limit
    _node(ax, 2.6, yt, 'separated limits:\nRP, gap, clustering', GREEN, w=2.6, h=0.95, fs=8.5)
    _node(ax, 6.0, yt, 'tempered,\nOS axioms', GREEN, w=1.9, h=0.95, fs=8.5)
    _node(ax, 9.6, yt, 'convergent,\nO(4)-invariant', GREEN, w=2.2, h=0.95, fs=8.5)
    # proved implications
    _arrow(ax, (1.9, yb), (2.2, ym), GREEN)
    ax.text(1.95, 2.45, '$n=2$;\nevery $n$ with (GD)', fontsize=8, color=GREEN, ha='right')
    _arrow(ax, (6.0, yb), (6.0, ym), GREEN)
    ax.text(6.12, 2.45, '$n=2$;\nevery $n$ with (GD)', fontsize=8, color=GREEN)
    _arrow(ax, (8.4, yb), (8.4, ym), GREEN)
    _arrow(ax, (10.7, yb), (10.7, ym), GREEN)
    ax.annotate('', xy=(3.45, ym), xytext=(2.75, ym), arrowprops=dict(arrowstyle='<|-|>', color=GREEN, lw=1.0))
    ax.annotate('', xy=(5.45, ym), xytext=(4.6, ym), arrowprops=dict(arrowstyle='<|-|>', color=GREEN, lw=1.0))
    _arrow(ax, (2.2, ym), (2.6, yt - 0.2), GREEN)
    _arrow(ax, (6.0, ym), (6.0, yt - 0.2), GREEN)
    _arrow(ax, (8.4, ym), (9.3, yt - 0.2), GREEN)
    _arrow(ax, (10.7, ym), (9.75, yt - 0.2), GREEN)
    # obstructions
    _arrow(ax, (1.9, yb), (4.0, ym), INK, cross=True, rad=0.05)
    ax.text(3.25, 2.05, r'$m\geq2$', fontsize=8, color=RED)
    _arrow(ax, (6.0, ym), (8.4, ym), INK, cross=True)
    ax.text(12.85, 8.05, 'the averaged condition along the tuning is assumed throughout', fontsize=8, color='#555555', ha='right')
    ax.plot([11.5, 13.0], [0.0, 0.0], color='none')
    ax.text(0.6, -0.35, 'red: open inputs     blue: lattice bounds     green arrows: proved     crossed: fails in examples',
            fontsize=8, color='#444444')
    save(fig, 'fig_limitlogic.png')


if __name__ == '__main__':
    fig_bptower()
    fig_limitlogic()
