"""Figure for Proposition prop:curvedsharp: the ratio of the covariant
Hessian of the plaquette action to half the covariant plaquette form, at the
field u of the proposition and the SU(2) background of constant curvature in
the plane {0,1}, from the closed forms of its proof. (a) The ratio as a
surface over the window side n and kappa (n+1)^2, for the two orientations;
(b) the section n = 63 with the lower bound of part (i) and the threshold 31
of part (ii).
Run from this directory: python3 make_curvedsharp.py. Writes
fig_curvedsharp.png (flattened, 256 colours)."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, INK = '#3b6fb6', '#c0392b', '#222222'


def ratio(n, K, sig):
    """Hess P(0)[u,u] / (Q^{V0}(u)/2) at kappa (n+1)^2 = K."""
    x = np.pi / (2 * n + 2)
    psi = ((n + 1) / 2) ** 2
    s2 = K / 2 / (n + 1) ** 2
    c2 = 1 - s2
    cot2 = 1 / np.tan(x) ** 2
    sigma = 2 / np.tan(x) - np.sin(2 * x)
    u2 = n * (n + 1) * psi
    q1 = 12 * np.sin(x) ** 2 * u2
    qv = q1 + psi * s2 * (4 * n * (n + 1) * np.cos(2 * x) + 8 * sig * c2 * cot2)
    hess = (qv / 2 - 4 * s2 ** 2 * n * (n + 1) * np.cos(x) ** 2 * psi
            + 2 * sig * s2 * c2 * (sigma ** 2 - 2 * cot2) * psi)
    return hess / (qv / 2)


ns = np.arange(64, 129, 4)
Ks = np.linspace(0, 40, 41)
N, KK = np.meshgrid(ns, Ks)
Rp = ratio(N, KK, 1)
Rm = ratio(N, KK, -1)

fig = plt.figure(figsize=(10, 4.2))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_surface(N, KK, Rp, color=RED, alpha=0.7, edgecolor='#7b241c', linewidth=0.25)
ax.plot_surface(N, KK, Rm, color=BLUE, alpha=0.6, edgecolor=BLUE, linewidth=0.25)
ax.plot_surface(N, KK, 0 * KK, color='#999999', alpha=0.15, linewidth=0)
ax.set_xlabel(r'window side $n$', labelpad=4)
ax.set_ylabel(r'$\kappa(n+1)^2$', labelpad=4)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel(r'$2\,\mathcal{H}P(u,u)/Q^{V_0}(u)$', labelpad=8, rotation=90)
ax.view_init(elev=20, azim=-50)
ax.set_box_aspect((1.1, 1.2, 0.85))
ax.text2D(0.22, 0.84, r'$\varsigma=+1$', color=RED, fontsize=11, transform=ax.transAxes)
ax.text2D(0.58, 0.33, r'$\varsigma=-1$', color=BLUE, fontsize=11, transform=ax.transAxes)
ax.set_title('(a)', loc='left')

ax2 = fig.add_subplot(1, 2, 2)
n = 63
K = np.linspace(0, 40, 801)
ax2.plot(K, ratio(n, K, 1), color=RED, lw=1.6, label=r'$\varsigma=+1$')
ax2.plot(K, ratio(n, K, -1), color=BLUE, lw=1.6, label=r'$\varsigma=-1$')
Kb = K[K <= 2]
ax2.plot(Kb, 1 + 0.037 * Kb, color=INK, lw=0.9, ls='--', label=r'$1+0.037\,\kappa(n+1)^2$, part (i)')
ax2.axhline(0, color='#888888', lw=0.6)
ax2.axvline(31, color=INK, lw=0.6, ls=':')
ax2.text(31.6, 0.62, r'$\kappa(n+1)^2=31$', fontsize=9)
ax2.set_xlabel(r'$\kappa(n+1)^2$, at $n=63$')
ax2.set_ylabel(r'$2\,\mathcal{H}P(u,u)/Q^{V_0}(u)$')
ax2.legend(frameon=False, fontsize=8.5, loc='lower left')
ax2.spines[['top', 'right']].set_visible(False)
ax2.set_title('(b)', loc='left')
fig.tight_layout()
save(fig, 'fig_curvedsharp.png')
print('ratio at n=63, K=31, sigma=-1:', ratio(63, 31, -1), ' zero crossing:',
      K[np.argmax(ratio(n, K, -1) < 0)])
