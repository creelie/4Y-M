"""Schematic figures of the continuum part of the forest-gauge paper.
Run from this directory: python3 make_schematics.py. Writes flattened 256-colour PNG files."""
import os, itertools, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from scipy.special import ive
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
GREEN, BLUE, RED, GREY = '#1e8449', '#2e6fb7', '#c0392b', '0.55'


def cuboid(ax, o, s, color, alpha=0.35, lw=0.6, edge=None):
    x, y, z = o; dx, dy, dz = s
    v = np.array([[x, y, z], [x + dx, y, z], [x + dx, y + dy, z], [x, y + dy, z],
                  [x, y, z + dz], [x + dx, y, z + dz], [x + dx, y + dy, z + dz], [x, y + dy, z + dz]])
    faces = [[v[i] for i in f] for f in ([0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4],
                                          [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4])]
    ax.add_collection3d(Poly3DCollection(faces, facecolors=color, edgecolors=edge or color,
                                         linewidths=lw, alpha=alpha))


def clean(ax, aspect=(1, 1, 0.8)):
    ax.set_axis_off()
    ax.set_box_aspect(aspect)


def fig_logic():
    """Three layers: proved at fixed spacing, proved implications, open inputs."""
    fig = plt.figure(figsize=(7.4, 5.0))
    ax = fig.add_axes([0.0, 0.0, 0.62, 1.0], projection='3d'); clean(ax, (1.6, 1, 1.25))
    ax.computed_zorder = False
    layers = [(0.0, GREEN), (1.6, BLUE), (3.2, RED)]
    for z, c in layers:
        P = [[(-0.3, -0.3, z), (6.3, -0.3, z), (6.3, 2.3, z), (-0.3, 2.3, z)]]
        ax.add_collection3d(Poly3DCollection(P, facecolors=c, alpha=0.06, edgecolors=c, linewidths=0.5))
        for i in range(3):
            cuboid(ax, (2.1 * i, 0.4, z), (1.6, 1.2, 0.45), c, alpha=0.45)
    for i in range(3):
        x = 2.1 * i + 1.25
        ax.plot([x, x], [1.0, 1.0], [3.2, 2.12], color='0.25', lw=1.1, zorder=8)
        ax.plot([x], [1.0], [2.12], marker='v', color='0.25', ms=5, zorder=8)
    for i in range(3):
        for z, c, s in ((0.0, GREEN, 'P'), (1.6, BLUE, 'I'), (3.2, RED, 'O')):
            ax.text(2.1 * i + 0.55, 1.0, z + 0.68, '%s%d' % (s, i + 1), ha='center', fontsize=9.5,
                    color=c, weight='bold', zorder=10)
    ax.set_xlim(-0.3, 6.3); ax.set_ylim(-0.3, 2.3); ax.set_zlim(0, 4.2)
    ax.view_init(elev=14, azim=-72)
    lines = [
        (RED, 'open inputs'),
        (RED, 'O1  cube estimate (D) along the tuning'),
        (RED, 'O2  slice ratio in four dimensions'),
        (RED, 'O3  tree bounds (T1)-(T3)'),
        (BLUE, 'proved implications'),
        (BLUE, 'I1  (D) gives a gap for every limit'),
        (BLUE, 'I2  ratio iff nontrivial flowed limit'),
        (BLUE, 'I3  tree bounds give the ratio'),
        (GREEN, 'proved at fixed spacing'),
        (GREEN, 'P1  gap for every simple group'),
        (GREEN, 'P2  SU(2), SU(3) with blocks'),
        (GREEN, 'P3  (D) at strong coupling'),
        ('0.25', 'arrows: an open input and the'),
        ('0.25', 'implication that consumes it'),
    ]
    y = 0.93
    for k, (c, s) in enumerate(lines):
        bold = k % 4 == 0 and k < 12
        fig.text(0.6, y, s, color=c, fontsize=9.5 if bold else 9, weight='bold' if bold else 'normal')
        y -= 0.058 if not (k % 4 == 3 and k < 12) else 0.1
    save(fig, 'fig_logic.png')


def fig_packet():
    """The packet f_0, its reflection and its translate in a three-dimensional slice of spacetime."""
    fig = plt.figure(figsize=(7.0, 4.4))
    ax = fig.add_axes([0.0, 0.08, 1.0, 0.92], projection='3d'); ax.set_axis_off()
    ax.set_box_aspect((1, 1, 1.3), zoom=1.2)
    ax.computed_zorder = False
    L = 2.2
    P = [[(-L, -L, 0), (L, -L, 0), (L, L, 0), (-L, L, 0)]]
    ax.add_collection3d(Poly3DCollection(P, facecolors='0.7', alpha=0.25, edgecolors='0.4', linewidths=0.6))
    u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 20))
    def blob(tc, color, alpha):
        x = 0.9 * np.cos(u) * np.sin(v); y = 0.9 * np.sin(u) * np.sin(v); t = tc + 0.35 * np.cos(v)
        ax.plot_surface(x, y, t, color=color, alpha=alpha, linewidth=0, shade=True)
    blob(1.05, BLUE, 0.55)
    blob(-1.05, GREY, 0.45)
    blob(2.55, GREEN, 0.5)
    ax.plot([1.6, 1.6], [1.6, 1.6], [-0.7, 1.4], color='0.3', lw=0.8, ls=':')
    for t, lab in ((0.7, r'$t_0$'), (1.4, r'$t_1$')):
        ax.plot([1.45, 1.75], [1.6, 1.6], [t, t], color='0.3', lw=0.9)
        ax.text(1.85, 1.6, t, lab, fontsize=10, va='center')
    ax.plot([-1.3, -1.3], [-1.3, -1.3], [1.05, 2.55], color='k', lw=1.0)
    ax.text(-1.55, -1.3, 1.8, r'$\sigma$', fontsize=11, ha='right')
    for y, lab, c in ((0.817, r'$\tau_\sigma f_0$', GREEN), (0.582, r'$f_0$', BLUE),
                      (0.26, r'$\vartheta f_0$', '0.25')):
        fig.text(0.6, y, lab, color=c, fontsize=11, va='center')
    fig.text(0.215, 0.425, r'$t=0$', fontsize=10, color='0.25')
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-1.6, 3.3)
    ax.view_init(elev=12, azim=-55)
    fig.text(0.02, 0.04, r'$C_k(\sigma)=\langle\Psi_{A_k},\mathbb{T}_k^{\sigma/a_k}\Psi_{A_k}\rangle$:'
             r' the reflected packet $\vartheta f_0$ paired with the packet translated by $\sigma$',
             fontsize=9.5)
    save(fig, 'fig_packet.png')


def su2_energy(a, j, g2=1.0, ell=1.0):
    beta = 4.0 / (g2 * a * a)
    return -(ell / a ** 2) * math.log(ive(int(2 * j + 1), beta) / ive(1, beta))


def fig_twod():
    """Two-dimensional SU(N): the cylinder lattice and the convergence of the energies."""
    fig = plt.figure(figsize=(7.4, 3.5))
    ax = fig.add_axes([0.0, 0.02, 0.47, 0.96], projection='3d'); clean(ax, (1, 1, 1.5))
    ax.computed_zorder = False
    L, T = 12, 6
    th = 2 * np.pi * np.arange(L) / L
    segs = []
    for t in range(T):
        for i in range(L):
            p = (math.cos(th[i]), math.sin(th[i]), t); q = (math.cos(th[(i + 1) % L]), math.sin(th[(i + 1) % L]), t)
            segs.append([p, q])
            if t + 1 < T:
                segs.append([p, (p[0], p[1], t + 1)])
    ax.add_collection3d(Line3DCollection(segs, colors='0.72', linewidths=0.7))
    t = 2
    ring = [(math.cos(s), math.sin(s), t) for s in np.linspace(0, 2 * np.pi, 120)]
    ax.plot(*zip(*ring), color=BLUE, lw=2.4)
    ring = [(math.cos(s), math.sin(s), t + 1) for s in np.linspace(0, 2 * np.pi, 120)]
    ax.plot(*zip(*ring), color=GREEN, lw=2.4)
    for i in (1, 2, 3):
        a0, a1 = th[i], th[i + 1]
        quad = [[(math.cos(a0), math.sin(a0), t), (math.cos(a1), math.sin(a1), t),
                 (math.cos(a1), math.sin(a1), t + 1), (math.cos(a0), math.sin(a0), t + 1)]]
        ax.add_collection3d(Poly3DCollection(quad, facecolors=RED, alpha=0.45, edgecolors=RED))
    ax.text2D(0.02, 0.62, r'$V^{\prime}=QV$', color=GREEN, fontsize=11, transform=ax.transAxes)
    ax.text2D(0.12, 0.40, r'$V$', color=BLUE, fontsize=11, transform=ax.transAxes)
    ax.text2D(0.78, 0.56, r'$W_x$', color=RED, fontsize=11, transform=ax.transAxes)
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_zlim(0, T - 1)
    ax.view_init(elev=18, azim=-60)
    ax.text2D(0.0, 0.95, '(a) holonomy on the circle', transform=ax.transAxes)
    bx = fig.add_axes([0.57, 0.17, 0.4, 0.7])
    a = np.geomspace(0.002, 0.25, 60)
    for j, c in ((0.5, BLUE), (1.0, GREEN), (1.5, RED)):
        E = [su2_energy(x, j) for x in a]
        bx.plot(a, E, color=c, lw=1.6, label=r'$j=%s$' % ('1/2' if j == 0.5 else ('1' if j == 1 else '3/2')))
        bx.axhline(j * (j + 1) / 2, color=c, lw=0.6, ls='--')
    bx.set_xscale('log'); bx.set_xlabel(r'lattice spacing $a$'); bx.set_ylabel(r'$E_{j,k}$')
    bx.set_title(r'(b) $\mathrm{SU}(2)$, $g^2=\ell=1$; limits $j(j+1)/2$', fontsize=9.5, loc='left')
    bx.legend(frameon=False, fontsize=9, loc='upper left')
    bx.set_ylim(0, 2.3)
    save(fig, 'fig_twod.png')


def fig_slice():
    """D_k(r)/D_k(r_0) over (r, a): two-dimensional SU(2) against a field whose mass is fixed in lattice units."""
    fig = plt.figure(figsize=(7.4, 3.6))
    r = np.linspace(0.5, 2.0, 30); la = np.linspace(1, 8, 30)
    R, LA = np.meshgrid(r, la); A = 2.0 ** (-LA)
    E2 = np.vectorize(lambda x: su2_energy(x, 0.5))(A)
    Z1 = np.exp(-(R - 0.5) * E2)
    mk = np.arccosh(1 + 0.5 ** 2 / 2) / A
    Z2 = np.exp(-(R - 0.5) * mk)
    for k, (Z, title, c) in enumerate(((Z1, r'(a) $\mathrm{SU}(2)$ in two dimensions', BLUE),
                                        (Z2, r'(b) mass $1/2$ in lattice units', RED))):
        ax = fig.add_subplot(1, 2, k + 1, projection='3d')
        ax.plot_surface(R, LA, Z, color=c, alpha=0.75, linewidth=0.2, edgecolor='w', shade=True)
        ax.set_xlabel(r'$r$', labelpad=-6); ax.set_ylabel(r'$\log_2(1/a)$', labelpad=-4)
        ax.tick_params(labelsize=7, pad=-2)
        ax.set_zlim(0, 1.05)
        ax.view_init(elev=24, azim=-128)
        ax.text2D(0.0, 0.97, title, transform=ax.transAxes, fontsize=9.5)
    fig.text(0.5, 0.02, r'height: $D_k(r)/D_k(1/2)$ for $g^2=\ell=1$', ha='center', fontsize=9.5)
    fig.subplots_adjust(left=0.0, right=0.97, wspace=0.05, top=0.97, bottom=0.08)
    save(fig, 'fig_slice.png')


def fig_blocking():
    """Block averages between the reflection plane and the packet, at three block sizes."""
    fig = plt.figure(figsize=(7.4, 3.4))
    D = 4.0
    for k, (n, lab, c) in enumerate(((4, r'$b=d/4$', GREEN), (2, r'$b=d/2$', BLUE), (1, r'$b=d$', RED))):
        ax = fig.add_axes([k / 3, 0.1, 1 / 3, 0.88], projection='3d'); clean(ax, (1, 1, 1.25))
        ax.computed_zorder = False
        P = [[(0, 0, 0), (D, 0, 0), (D, D, 0), (0, D, 0)]]
        ax.add_collection3d(Poly3DCollection(P, facecolors='0.7', alpha=0.3, edgecolors='0.4', linewidths=0.6))
        s_ = D / n
        for i, j, l in itertools.product(range(n), range(n), range(n)):
            if i == n - 1 or j == 0 or l == n - 1:
                cuboid(ax, (i * s_, j * s_, l * s_), (s_, s_, s_ * 0.98), c, alpha=0.07, lw=0.5, edge=c)
        u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 20))
        ax.plot_surface(D / 2 + 1.1 * np.cos(u) * np.sin(v), D / 2 + 1.1 * np.sin(u) * np.sin(v),
                        D + 0.75 + 0.35 * np.cos(v), color='#8e44ad', alpha=0.6, linewidth=0)
        ax.plot([D, D], [D, D], [0, D + 0.75], color='k', lw=1.0)
        ax.set_xlim(0, D); ax.set_ylim(0, D); ax.set_zlim(0, D + 1.2)
        ax.view_init(elev=14, azim=-58)
        ax.text2D(0.42, 0.02, lab, color=c, fontsize=11, transform=ax.transAxes)
        if k == 0:
            ax.text2D(0.72, 0.76, r'$F$', color='#8e44ad', fontsize=11, transform=ax.transAxes)
            ax.text2D(0.91, 0.5, r'$d$', fontsize=11, transform=ax.transAxes)
            ax.text2D(0.02, 0.2, r'$t=0$', fontsize=10, color='0.25', transform=ax.transAxes)
    fig.text(0.5, 0.02, r'block averages of side $b$ between the plane $t=0$ and the packet $F$ at height $d$',
             ha='center', fontsize=9.5)
    save(fig, 'fig_blocking.png')


if __name__ == '__main__':
    fig_logic(); fig_packet(); fig_twod(); fig_slice(); fig_blocking()
