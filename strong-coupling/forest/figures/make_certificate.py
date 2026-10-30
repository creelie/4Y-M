"""Figure of Appendix "app:certificate": (a) the least eigenvalue of the plaquette matrix M_p(phi)
over the torus of phases, for the plaquette whose minimum is smallest; (b) the bracket of the definition
of h for the 48 free links of a period cell, against 9 + sqrt(17), 14.384 and 29/2; (c) the constants of
the block gap, Theorem "blockgap" (52), Proposition "blockgapsharp" (5000/537) and the bulk Bloch value
1/c* = 1/0.1074877694868. Run from this directory: python3 make_certificate.py (about one minute; it
imports strong-coupling/gauge-fixed/checks/romanov/romanov_certificate.py and its data). The
eigenvalues are in floating point. Writes fig_certificate.png."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save

HERE = os.path.dirname(os.path.abspath(__file__))
ROM = os.path.join(HERE, '..', '..', 'gauge-fixed', 'checks', 'romanov')
sys.path.insert(0, ROM)
sys.argv = [sys.argv[0], os.path.join(ROM, 'cert4d_h14.42.json')]
import romanov_certificate as rc  # noqa: E402

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, ORANGE, GREEN, INK = '#3b6fb6', '#c0392b', '#e67e22', '#2e8b57', '#222222'


def plaquette(j):
    p, c, f = rc.PL[j]
    k = [2 * float(rc.G[s][r][q]) for s, r, q in c]
    idx = [4 * i + q for i in f for q in range(4)]
    R = lambda M: M[np.ix_(idx, idx)]
    D = np.zeros((16, 16))
    for i in f:
        D[4 * i:4 * i + 4, 4 * i:4 * i + 4] = float(rc.A[j][i]) * np.eye(4)
    K0 = R(D - k[0] * rc.C14 - k[1] * rc.C12 - k[3] * rc.C34 - k[2] * rc.CB[0])
    F = [R(x) for x in rc.FB]

    def lmin(p1, p2):
        dl = p1 - p2
        M = (K0[None] - np.cos(p1)[:, None, None] * F[0] - np.sin(p1)[:, None, None] * F[1]
             - np.cos(p2)[:, None, None] * F[2] - np.sin(p2)[:, None, None] * F[3]
             - k[2] * (np.cos(dl)[:, None, None] * R(rc.CB[1]) + np.sin(dl)[:, None, None] * R(rc.CB[2])))
        return np.linalg.eigvalsh(M)[:, 0]
    return lmin


n = 96
t = np.linspace(0, 2 * np.pi, n + 1)
P1, P2 = np.meshgrid(t, t, indexing='ij')
best = None
for j in range(len(rc.PL)):
    lam = plaquette(j)(P1.ravel(), P2.ravel()).reshape(P1.shape)
    if best is None or lam.min() < best[1].min():
        best = (j, lam)
j, lam = best
print('plaquette', j, 'least eigenvalue on the grid', lam.min())

fig = plt.figure(figsize=(9.6, 3.2))
ax = fig.add_axes([-0.03, 0.0, 0.38, 1.0], projection='3d')
ax.set_box_aspect((1, 1, 0.75))
ax.plot_surface(P1, P2, lam, cmap='viridis', rstride=2, cstride=2, lw=0.1, edgecolor='#444444', alpha=0.95)
i0 = np.unravel_index(lam.argmin(), lam.shape)
ax.scatter([P1[i0]], [P2[i0]], [lam[i0]], color=RED, s=18, depthshade=False)
ax.set_xlabel(r'$\varphi_1$', labelpad=-6)
ax.set_ylabel(r'$\varphi_2$', labelpad=-6)
ax.set_xticks([0, np.pi, 2 * np.pi]); ax.set_xticklabels(['0', r'$\pi$', r'$2\pi$'])
ax.set_yticks([0, np.pi, 2 * np.pi]); ax.set_yticklabels(['0', r'$\pi$', r'$2\pi$'])
ax.tick_params(pad=-2, labelsize=8)
ax.view_init(26, -55)
ax.set_title(r'(a) least eigenvalue of $M_p(\varphi)$', fontsize=10, y=0.97)

tot = sorted(float(v) for v in rc.tot.values())
ax2 = fig.add_axes([0.40, 0.18, 0.27, 0.68])
ax2.bar(range(len(tot)), tot, color=BLUE, width=0.8)
for y, col, lab in ((14.5, RED, r'$\frac{29}{2}$'), (14.420001, ORANGE, r'$14.420001$'),
                    (14.384, GREEN, r'$14.384$'), (9 + np.sqrt(17), INK, r'$9+\sqrt{17}$')):
    ax2.axhline(y, color=col, lw=1.0, ls='--' if col != ORANGE else '-')
    ax2.text(len(tot) + 0.5, y, lab, color=col, fontsize=8, va='center')
ax2.set_ylim(12.6, 14.7)
ax2.set_xlim(-1, len(tot))
ax2.set_xlabel('free links of a period cell')
ax2.set_title('(b) totals of the certificate', fontsize=10)
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)

ax3 = fig.add_axes([0.80, 0.18, 0.18, 0.68])
vals = [(1 / 0.1074877694868, GREEN, 'Bloch'), (5000 / 537, ORANGE, r'$\frac{5000}{537}$'), (52, RED, '52')]
for x, (v, col, lab) in enumerate(vals):
    ax3.bar(x, v, color=col, width=0.6)
    ax3.text(x, v * 1.08, lab, ha='center', fontsize=9, color=col)
ax3.set_yscale('log')
ax3.set_ylim(5, 100)
ax3.set_xticks([])
ax3.set_ylabel(r'$C$ in $\sum\phi^2\leq C\,Q$')
ax3.set_title('(c) block gap', fontsize=10)
ax3.spines['top'].set_visible(False); ax3.spines['right'].set_visible(False)
save(fig, 'fig_certificate.png')
