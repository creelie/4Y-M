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
            ('SU(3)', [('SZZ', 9 / 48, '3/16'), ('curv', 9 / 32, '9/32'), ('fcurv', 9 / 29, '9/29'), ('fdob', 3 / 4, '3/4'), ('fblock', 3 / 4, '3/4')])]
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

if __name__ == '__main__':
    fig_linefield(); fig_ranges(); fig_cover(); fig_face()
