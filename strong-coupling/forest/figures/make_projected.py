"""Figure for Definition def:projavg and Theorem thm:projinverse:
(a) the projected average on a sphere standing for SU(2): the paths P_x
around the base path P_{2y}, their projections to the tangent plane (the Lie
algebra), the mean A of the projections, and the point P_{2y} exp A;
(b) the map g -> c exp(-A(g)) contracts with constant 7/8, so its iterates
converge to the unique solution G of the equation hat V = c, and hat V moves
distances by a factor between 1/8 and 15/8. Schematic, not to scale.
Run from this directory: python3 make_projected.py. Writes fig_projected.png."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1'
rng = np.random.default_rng(3)

fig = plt.figure(figsize=(9.0, 4.0))
ax = fig.add_axes([-0.04, 0.02, 0.56, 0.98], projection='3d')
clean3d(ax)
ax.set_box_aspect((1, 1, 1))
u, v = np.mgrid[0:2 * np.pi:60j, 0:np.pi:30j]
ax.plot_surface(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v),
                color='#dfe6ee', alpha=0.25, lw=0, shade=True)
base = np.array([0.0, 0.0, 1.0])
# tangent plane at the base point
s = np.linspace(-1.0, 1.0, 2)
X, Y = np.meshgrid(s, s)
ax.plot_surface(X, Y, np.ones_like(X), color=BLUE, alpha=0.12, lw=0)
ax.scatter(*base, color=INK, s=28, depthshade=False)
ax.text(-0.45, -0.25, 1.02, r'$P_{2y}$', color=INK, fontsize=11)
pts = []
for k in range(7):
    t = rng.uniform(0.55, 0.95); ph = -1.1 + 2.9 * k / 6 + rng.uniform(-0.15, 0.15)
    p = np.array([np.sin(t) * np.cos(ph), np.sin(t) * np.sin(ph), np.cos(t)])
    q = np.array([p[0], p[1], 1.0])
    pts.append(q)
    ax.scatter(*p, color=GREEN, s=14, depthshade=False)
    ax.plot(*zip(p, q), color=GREY, lw=0.8, ls=':')
    ax.scatter(*q, color=BLUE, s=10, depthshade=False)
    ax.plot(*zip(base, q), color=BLUE, lw=0.6, alpha=0.6)
A = np.mean(pts, axis=0)
ax.quiver(*base, *(A - base), color=RED, lw=1.8, arrow_length_ratio=0.25)
r = np.linalg.norm(A[:2])
tip = np.array([A[0] / r * np.sin(r), A[1] / r * np.sin(r), np.cos(r)])
ax.scatter(*tip, color=RED, s=30, depthshade=False)
ax.text(*(A + np.array([0.04, 0.04, 0.08])), r'$A$', color=RED, fontsize=11)
ax.text(*(tip + np.array([0.55, 0.35, -0.05])), r'$\widehat V=P_{2y}e^{A}$', color=RED, fontsize=11)
ax.text(0.25, -1.05, 0.55, r'paths $P_x$', color=GREEN, fontsize=10)
ax.text(-1.0, 0.9, 1.08, r'tangent space', color=BLUE, fontsize=10)
ax.view_init(elev=38, azim=-58)
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-0.6, 1.2)
fig.text(0.02, 0.93, '(a)', fontsize=11)

bx = fig.add_axes([0.57, 0.14, 0.40, 0.76])
g = np.array([1.6, 1.2]); G = np.array([0.0, 0.0])
traj = [g]
rot = np.array([[np.cos(0.5), -np.sin(0.5)], [np.sin(0.5), np.cos(0.5)]])
for _ in range(12):
    g = 0.875 * rot @ g
    traj.append(g)
traj = np.array(traj)
for rad in (0.4, 0.8, 1.2, 1.6, 2.0):
    th = np.linspace(0, 2 * np.pi, 200)
    bx.plot(rad * np.cos(th), rad * np.sin(th), color=GREY, lw=0.4, alpha=0.6)
bx.plot(traj[:, 0], traj[:, 1], color=BLUE, lw=1.0, marker='o', ms=3)
bx.scatter(*G, color=RED, s=40, zorder=5)
bx.annotate(r'$G$', G, xytext=(-0.45, -0.35), color=RED)
bx.annotate(r'$g$', traj[0], xytext=(traj[0][0] + 0.08, traj[0][1] + 0.08), color=BLUE)
bx.text(-2.0, 1.75, r'$g\mapsto c_k\,e^{-A(g)}$, constant $\frac{7}{8}$', fontsize=9)
bx.text(-2.0, -2.05, r'$\frac{1}{8}\,d(g,g^\prime)\leq d(\widehat V_k g,\widehat V_k g^\prime)\leq\frac{15}{8}\,d(g,g^\prime)$', fontsize=9)
bx.set_aspect('equal'); bx.set_xlim(-2.2, 2.2); bx.set_ylim(-2.3, 2.1); bx.axis('off')
fig.text(0.56, 0.93, '(b)', fontsize=11)
save(fig, 'fig_projected.png')
