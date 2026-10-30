"""Figure of Theorem "effectivelocal": (a) the multiplier K_P of the effective form S_3
for the plane P = {1,2}, as a function of the two transverse momenta (theta_3, theta_4) at
theta_1 = theta_2 = 0, equal to 1 at the origin; (b) the largest value of the kernel of K_P
at l1 distance d on the coarse lattice, k = 1, 2, 3, in floating point on a grid of 12^4
momenta. Run from this directory: python3 make_locality.py (about four minutes; it imports
the functions of strong-coupling/gauge-fixed/checks/blocks/effective_locality.py). Writes
fig_locality.png."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'gauge-fixed', 'checks', 'blocks'))
sys.path.insert(0, HERE)
import effective_locality as el

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, ORANGE, GREEN, INK = '#3b6fb6', '#c0392b', '#e67e22', '#2e8b57', '#222222'

fig = plt.figure(figsize=(8.0, 3.3))
ax = fig.add_axes([-0.04, 0.0, 0.56, 1.0], projection='3d')
ax.set_box_aspect((1, 1, 0.8))
t = np.linspace(-np.pi, np.pi, 41)
T3, T4 = np.meshgrid(t, t, indexing='ij')
th = np.stack([0 * T3, 0 * T3, T3, T4], -1).reshape(-1, 4)
at, Y = el.blocks(th, 3)
Kp = (at[:, 2] * at[:, 3] / Y).real.reshape(T3.shape)
ax.plot_surface(T3, T4, Kp, cmap='viridis', rstride=1, cstride=1, lw=0.1, edgecolor='#444444', alpha=0.95)
ax.scatter([0], [0], [1], color=RED, s=18, depthshade=False)
ax.set_xlabel(r'$\theta_3$', labelpad=-6)
ax.set_ylabel(r'$\theta_4$', labelpad=-6)
ax.set_zlabel(r'$K_P$', labelpad=-4)
ax.set_xticks([-3, 0, 3]); ax.set_yticks([-3, 0, 3])
ax.tick_params(pad=-2, labelsize=8)
ax.view_init(24, -58)
ax.set_title(r'(a) $K_P$, $P=\{1,2\}$, $k=3$', fontsize=10, y=0.97)

ax2 = fig.add_axes([0.62, 0.17, 0.36, 0.7])
N = 12
js = np.arange(N) * 2 * np.pi / N
X = np.stack(np.meshgrid(js, js, js, js, indexing='ij'), -1).reshape(-1, 4)
rr = np.indices((N,) * 4)
rr = np.minimum(rr, N - rr).sum(0)
for k, col in zip((1, 2, 3), (GREEN, BLUE, ORANGE)):
    at, Y = el.blocks(X, k)
    KP = (np.prod(at[:, 2:], -1) / Y).real
    Kr = np.fft.ifftn(KP.reshape((N,) * 4)).real
    mag = [np.abs(Kr)[rr == d].max() for d in range(7)]
    ax2.semilogy(range(7), mag, 'o-', color=col, ms=4, lw=1.4, label=r'$k=%d$' % k)
ax2.set_xlabel(r'$|y|_1$')
ax2.set_ylabel(r'$\max|\mathcal{K}_P(y)|$')
ax2.legend(frameon=False, fontsize=9)
ax2.set_title('(b) kernel of $K_P$ (floating point)', fontsize=10)
ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
save(fig, 'fig_locality.png')
