"""Figures of the forest-gauge paper. Run from this directory: python3 make_figures.py.
Reads the certified boxes in ../../gauge-fixed/checks/su3/deg2/r3_2/cover_boxes.txt and the face
values in ../../gauge-fixed/checks/su3/deg2/face2.txt. Writes flattened 256-colour PNG files."""
import os, itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
DEG2 = os.path.join(HERE, '..', '..', 'gauge-fixed', 'checks', 'su3', 'deg2')
plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
TABLE = [0, 0, 2, 3, 3, 1, 2, 1, 1, 2, 1, 3, 3, 2, 0, 0]
sigma = lambda x: TABLE[(x[0] % 2) + 2 * (x[1] % 2) + 4 * (x[2] % 2) + 8 * (x[3] % 2)]

def save(fig, name, width_px=1600):
    path = os.path.join(HERE, name)
    w_in = fig.get_size_inches()[0]
    fig.savefig(path, dpi=width_px / w_in, facecolor='white')
    plt.close(fig)
    im = Image.open(path).convert('RGBA')
    bg = Image.new('RGBA', im.size, (255, 255, 255, 255))
    bg.alpha_composite(im)
    bg.convert('RGB').quantize(colors=256, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).save(path, optimize=True)
    print('wrote', name, Image.open(path).size)

def clean3d(ax):
    ax.set_axis_off()
    ax.set_box_aspect((1, 1, 0.9))

def fig_linefield():
    fig = plt.figure(figsize=(7.2, 3.4))
    ax = fig.add_subplot(1, 2, 1, projection='3d'); clean3d(ax)
    ax.computed_zorder = False
    n = 3
    segs = []
    for x0, x1, x2 in itertools.product(range(n), repeat=3):
        for m in range(3):
            y = [x0, x1, x2]; y[m] += 1
            if y[m] < n:
                segs.append([(x0, x1, x2), tuple(y)])
    ax.add_collection3d(Line3DCollection(segs, colors='0.78', linewidths=0.6, zorder=1))
    for x0, x1, x2 in itertools.product(range(n), repeat=3):
        s = sigma((x0, x1, x2, 0))
        if s < 3:
            e = np.zeros(3); e[s] = 1.0; p0 = np.array([x0, x1, x2], float)
            ax.plot(*zip(p0, p0 + 1.0 * e), color='#c0392b', linewidth=3.0, solid_capstyle='butt', zorder=5)
            ax.scatter(*(p0 + 0.62 * e), s=16, marker='o', c='#c0392b', depthshade=False, zorder=6)
        else:
            ax.scatter([x0], [x1], [x2], s=34, marker='o', facecolors='#e67e22', edgecolors='#7e3f0a', linewidths=0.6, depthshade=False, zorder=6)
    ax.scatter(*np.array(list(itertools.product(range(n), repeat=3))).T, s=5, c='k', depthshade=False, zorder=7)
    ax.set_xlim(0, n - 0.4); ax.set_ylim(0, n - 0.4); ax.set_zlim(0, n - 0.4)
    ax.view_init(elev=22, azim=-58)
    ax.text2D(0.02, 0.95, '(a) the line field in the slice $x_3=0$', transform=ax.transAxes)
    ax = fig.add_subplot(1, 2, 2, projection='3d'); clean3d(ax)
    P = np.array([(0, 0, 0), (1.6, 0, 0), (1.6, 0, 1.6), (0, 0, 1.6)])
    edges = [(0, 1, '#c0392b', r'$e_1$ (fixed)'), (1, 2, '#2e6fb7', r'$e_2$:  $\frac{21}{10}$'),
             (3, 2, '#1e8449', r'$e_3$:  $\frac{5}{4}$'), (0, 3, '#2e6fb7', r'$e_4$:  $\frac{21}{10}$')]
    for i, j, c, lab in edges:
        a, b = P[i], P[j]
        ax.plot(*zip(a, b), color=c, linewidth=3.2)
        mid = (a + b) / 2; off = mid - P.mean(0)
        ax.text(*(mid + 0.38 * off / np.linalg.norm(off)), lab, color=c, ha='center', va='center', fontsize=10)
    for k in (2, 3):
        ax.scatter(*P[k], s=60, c='#8e44ad', depthshade=False)
    ax.plot([0.25, 1.35, 1.35, 0.25, 0.25], [0, 0, 0, 0, 0], [0.25, 0.25, 1.35, 1.35, 0.25], color='0.6', linewidth=0.6, linestyle=':')
    ax.text(0.8, 0, 0.8, r'$p$', ha='center', va='center', fontsize=13)
    ax.scatter(*P.T, s=12, c='k', depthshade=False)
    ax.set_xlim(-0.4, 2.0); ax.set_ylim(-1, 1); ax.set_zlim(-0.4, 2.0)
    ax.view_init(elev=18, azim=-70)
    ax.text2D(0.02, 0.95, '(b) one plaquette in the gauge', transform=ax.transAxes)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0)
    save(fig, 'fig_linefield.png')

def fig_ranges():
    from matplotlib.patches import Patch
    fig = plt.figure(figsize=(7.2, 4.4))
    ax = fig.add_axes([-0.06, 0.2, 1.1, 0.86], projection='3d')
    meth = {'SZZ': ('ungauged curvature (Shen, Zhu, Zhu)', '#d5d8dc'),
            'curv': ('ungauged curvature, sharp', '#95a5a6'),
            'quarter': ('ungauged Dobrushin, sharp quarter', '#566573'),
            'fcurv': ('forest gauge, curvature', '#2e6fb7'),
            'fdob': ('forest gauge, single-link Dobrushin', '#c0392b'),
            'fblock': ('forest gauge, plaquette blocks', '#e67e22')}
    col = {k: i for i, k in enumerate(meth)}
    rows = [('SU(2)', [('curv', 4 / 32, '1/8'), ('quarter', 2 / 9, '2/9'), ('fdob', 1 / 3, '1/3'), ('fblock', 0.345, '0.345')]),
            ('SU(3)', [('SZZ', 9 / 48, '3/16'), ('curv', 9 / 32, '9/32'), ('fcurv', 9 / 29, '9/29'), ('fdob', 3 / 4, '3/4'), ('fblock', 0.76, '0.76')])]
    ax.computed_zorder = False
    for yi, (name, items) in enumerate(rows):
        y0 = 1.5 * yi
        for key, v, lab in items:
            x0 = 1.15 * col[key]
            z = 10 + 10 * (1 - yi) + col[key]
            ax.bar3d(x0, y0, 0, 0.6, 0.6, v, color=meth[key][1], alpha=1.0, shade=True, edgecolor='0.3', linewidth=0.3, zorder=z).set_clip_on(False)
            ax.text(x0 + 0.3, y0 + 0.3, v + 0.02, f'${lab}$', ha='center', va='bottom', fontsize=9, zorder=100,
                    bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.85))
        ax.text(-0.35, y0 + 0.3, 0, name, fontsize=11, ha='right', va='center', zorder=100)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zlim(0, 0.8); ax.set_zticks([0, 0.25, 0.5, 0.75])
    ax.tick_params(labelsize=8)
    ax.set_xlim(-0.2, 6.6); ax.set_ylim(-0.2, 2.3)
    ax.set_box_aspect((2.8, 1.1, 1.0), zoom=1.25)
    ax.view_init(elev=22, azim=-62); ax.set_proj_type('ortho')
    ax.xaxis.pane.fill = ax.yaxis.pane.fill = False
    ax.set_zlabel(r'$\beta_W$', labelpad=2)
    fig.legend(handles=[Patch(color=c, label=l) for l, c in meth.values()], loc='lower center',
               bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=8.5, frameon=False)
    save(fig, 'fig_ranges.png')

def fig_cover():
    A = np.array([[float(t) for t in l.split()] for l in open(os.path.join(DEG2, 'r3_2', 'cover_boxes.txt')) if not l.startswith('#')])
    lo, hi, val = A[:, :4], A[:, 4:8], A[:, 9]
    c = (lo + hi) / 2
    r = 3 / 2
    fig = plt.figure(figsize=(7.2, 3.8))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    sc = ax.scatter(c[:, 0], c[:, 1], np.hypot(c[:, 2], c[:, 3]), c=val, cmap='viridis', vmin=0.84, vmax=1.0, s=10 * (hi[:, 0] - lo[:, 0]) / (r / 8) + 3, depthshade=False)
    ax.set_xlabel(r'$a_1$', labelpad=-6); ax.set_ylabel(r'$a_2$', labelpad=-6)
    ax.tick_params(labelsize=7, pad=-2)
    ax.view_init(elev=24, azim=-128)
    ax.set_title(f'(a) all {len(val)} boxes; height $|z_3|$', fontsize=10)
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    sel = val > 0.99
    ax.scatter(c[sel, 0], c[sel, 1], c[sel, 2], c=val[sel], cmap='viridis', vmin=0.84, vmax=1.0, s=12, depthshade=False)
    ax.set_xlabel(r'$a_1$', labelpad=-6); ax.set_ylabel(r'$a_2$', labelpad=-6)
    ax.tick_params(labelsize=7, pad=-2)
    ax.view_init(elev=24, azim=-128)
    ax.set_title('(b) bound above 0.99; height $x_3$', fontsize=10)
    cax = fig.add_axes([0.905, 0.22, 0.016, 0.56])
    cb = fig.colorbar(sc, cax=cax, ticks=[0.84, 0.88, 0.92, 0.96, 1.0])
    cb.set_label('certified bound on the box', fontsize=9); cb.ax.tick_params(labelsize=7)
    fig.subplots_adjust(left=0.05, right=0.87, bottom=0.02, top=0.92, wspace=0.08)
    save(fig, 'fig_cover.png')

def fig_face():
    data = {}
    for l in open(os.path.join(DEG2, 'face2.txt')):
        if not l.startswith('r='): continue
        a = l.replace(':', ' ').split()
        r = float(a[0][2:]); q = float(a[1].split('=')[1])
        v1 = float(a[2].split('=')[1]); v2 = float(a[3].split('=')[1])
        data.setdefault(r, []).append((q, v1, v2))
    fig = plt.figure(figsize=(7.2, 3.6))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    cols = {1.44: '#2e6fb7', 1.5: '#1e8449', 1.54: '#c0392b'}
    for r, pts in sorted(data.items()):
        pts = sorted(pts); q = np.array([p[0] for p in pts])
        v1 = np.array([p[1] for p in pts]); v2 = np.array([p[2] for p in pts])
        c = cols.get(r, 'k')
        ax.plot(q, r + 0 * q, v1, color=c, linewidth=1.4, linestyle='--')
        ax.plot(q, r + 0 * q, v2, color=c, linewidth=2.2, marker='o', markersize=3)
        for qq, vv in zip(q, v2):
            ax.plot([qq, qq], [r, r], [0.90, vv], color=c, linewidth=0.4, alpha=0.5)
    Q, R = np.meshgrid(np.linspace(0.5, 1.0, 2), np.linspace(1.42, 1.56, 2))
    ax.plot_surface(Q, R, 1 + 0 * Q, color='0.6', alpha=0.25)
    ax.set_xlabel(r'$\rho/r$', labelpad=-4); ax.set_ylabel(r'$r$', labelpad=-2); ax.set_zlabel('bound', labelpad=2)
    ax.tick_params(labelsize=7, pad=-2); ax.set_zlim(0.90, 1.10)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0]); ax.set_yticks([1.44, 1.50, 1.54]); ax.set_zticks([0.92, 0.96, 1.0, 1.04, 1.08])
    ax.view_init(elev=16, azim=-62)
    ax.set_box_aspect((2.2, 1.1, 1.0), zoom=1.12)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color='k', linewidth=2.2, marker='o', markersize=3, label='degree-two corrector'),
                        Line2D([], [], color='k', linewidth=1.4, linestyle='--', label='blended corrector')]
               + [Line2D([], [], color=c, linewidth=2.2, label=f'$r={r:.2f}$') for r, c in sorted(cols.items())],
               loc='upper left', bbox_to_anchor=(0.06, 0.97), fontsize=8.5, frameon=False)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    save(fig, 'fig_face.png')

GROUPS = [('SU(2)', 2, 1.0), ('SU(3)', 3, 1.5), ('Sp(2)', 4, 1.5), ('SO(5)', 5, 0.75), ('G$_2$', 7, 1.0),
          ('F$_4$', 26, 0.75), ('E$_6$', 27, 1.0), ('E$_7$', 56, 0.75), ('E$_8$', 248, 0.25)]

def fig_groups():
    """Gap m_G as a function of beta_W = N^2 beta for nine groups (Theorem thm:groups)."""
    fig = plt.figure(figsize=(7.2, 4.0))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    cmap = plt.get_cmap('viridis')
    for k, (name, N, kap) in enumerate(GROUPS):
        top = 2 * kap * N / 29  # range of beta_W
        bw = np.geomspace(top * 1e-3, top * (1 - 1e-6), 400)
        u = bw / (N * kap)  # N|beta| / kappa_G
        m = np.arccosh(1 + (1 / u - 29 / 2) / 12)
        x = np.log10(bw); y = np.full_like(x, k); c = cmap(k / (len(GROUPS) - 1))
        ax.plot(x, y, m, color=c, linewidth=1.8)
        verts = [list(zip(x, y, m)) + [(x[-1], k, 0), (x[0], k, 0)]]
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        ax.add_collection3d(Poly3DCollection(verts, facecolors=[c], alpha=0.28, edgecolors='none'))
        ax.scatter([np.log10(top)], [k], [0], color=c, s=14, depthshade=False)
    ax.set_yticks(range(len(GROUPS))); ax.set_yticklabels([g[0] for g in GROUPS], fontsize=7.5)
    ax.set_xticks([-3, -2, -1, 0]); ax.set_xticklabels(['$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '$1$'], fontsize=7.5)
    ax.set_xlabel(r'$\beta_W$', labelpad=0); ax.set_zlabel(r'gap $m_G$', labelpad=0)
    ax.tick_params(axis='z', labelsize=7.5)
    ax.set_zlim(0, 6.5)
    ax.view_init(elev=24, azim=-58)
    ax.set_box_aspect((1.7, 1.5, 0.9), zoom=1.1)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    save(fig, 'fig_groups.png')

def fig_cubes():
    """(a) averaged cubes: V (n=4) and translates by 2Z^d, drawn in three dimensions;
    (b) one cube of physical side l at three lattice spacings."""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(7.2, 3.5))
    ax = fig.add_subplot(1, 2, 1, projection='3d'); clean3d(ax)
    n, L = 4, 7
    def box(o, s, c, a, lw):
        o = np.array(o, float)
        V = [o + s * np.array(v) for v in itertools.product((0, 1), repeat=3)]
        faces = [[V[0], V[1], V[3], V[2]], [V[4], V[5], V[7], V[6]], [V[0], V[1], V[5], V[4]],
                 [V[2], V[3], V[7], V[6]], [V[0], V[2], V[6], V[4]], [V[1], V[3], V[7], V[5]]]
        ax.add_collection3d(Poly3DCollection(faces, facecolors=c, alpha=a, edgecolors=c, linewidths=lw))
    pts = np.array(list(itertools.product(range(L), repeat=3)))
    ax.scatter(*pts.T, s=2, c='0.55', depthshade=False)
    box((0, 0, 0), n - 1, '#2e6fb7', 0.10, 1.4)
    box((2, 0, 0), n - 1, '#e67e22', 0.06, 0.9)
    box((2, 2, 2), n - 1, '#1e8449', 0.06, 0.9)
    j0 = np.array([3.0, 3.0, 4.0])
    ax.plot(*zip(j0, j0 + np.array([1, 0, 0])), color='#c0392b', linewidth=3.2)
    ax.text(*(j0 + np.array([0.6, -0.2, 0.55])), r'$j$', color='#c0392b', fontsize=11)
    ax.text(0.2, -0.9, -0.6, r'$V$', color='#2e6fb7', fontsize=11)
    ax.text(5.2, -0.9, -0.6, r'$V+2e_1$', color='#e67e22', fontsize=9)
    ax.text2D(0.62, 0.86, r'$V+2(e_1+e_2+e_3)$', color='#1e8449', fontsize=9, transform=ax.transAxes)
    ax.set_xlim(0, L - 1); ax.set_ylim(0, L - 1); ax.set_zlim(0, L - 1)
    ax.view_init(elev=20, azim=-60)
    ax.text2D(0.0, 0.97, '(a) averaged cubes, $n=4$', transform=ax.transAxes)
    ax = fig.add_subplot(1, 2, 2, projection='3d'); clean3d(ax)
    ell = 1.0
    for k, (a, c) in enumerate([(0.25, '#2e6fb7'), (0.125, '#1e8449'), (0.0625, '#c0392b')]):
        z = 1.2 * k
        g = np.arange(0, 1.5 + 1e-9, a)
        segs = [[(x, 0, z), (x, 1.5, z)] for x in g] + [[(0, y, z), (1.5, y, z)] for y in g]
        ax.add_collection3d(Line3DCollection(segs, colors=c, linewidths=0.35, alpha=0.55))
        sq = [(0.25, 0.25, z), (0.25 + ell, 0.25, z), (0.25 + ell, 0.25 + ell, z), (0.25, 0.25 + ell, z)]
        ax.add_collection3d(Poly3DCollection([sq], facecolors=c, alpha=0.18, edgecolors=c, linewidths=1.8))
        lab = (f'$a/{2**k}$' if k else '$a$') + f', $n={int(round(ell / a))}$'
        ax.text(1.62, 0.0, z, lab, color=c, fontsize=9.5)
    ax.plot([0.25, 0.25 + ell], [0.12, 0.12], [-0.05, -0.05], color='k', linewidth=1.0)
    ax.text(0.55, -0.25, -0.35, r'$\ell=na$', fontsize=10)
    ax.set_xlim(0, 1.6); ax.set_ylim(0, 1.6); ax.set_zlim(-0.3, 2.6)
    ax.view_init(elev=18, azim=-62)
    ax.text2D(0.0, 0.97, '(b) one physical scale', transform=ax.transAxes)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0)
    save(fig, 'fig_cubes.png')

if __name__ == '__main__':
    fig_linefield(); fig_ranges(); fig_cover(); fig_face(); fig_groups(); fig_cubes()
