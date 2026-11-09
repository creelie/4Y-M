"""Figure for the cube condition of free massive fields (ssec:fffree).

fig_freefield: (a) the response of one cube to a constant shift of its
boundary data, the function mu^2 (-Delta_D + mu^2)^{-1} 1 on a square cube of
side n = 14 with s = n mu = 2, drawn as a surface; its mean over the cube is
the exact constant of Theorem ffs:thm:exact (in d = 2). (b) the exact
constant varpi_n(s) = mu^2 n^{-d} <1_V, (-Delta_D + mu^2)^{-1} 1_V> in d = 4
for n = 2, 4, 8 and its limit s^2 <1, (-Delta + s^2)^{-1} 1> on the unit
cube (Theorem ffs:thm:limit), against s = m l, with the lower bound
(5/6)^4 s^2/(s^2 + 40) of Theorem thm:ffscalar.

Run from this directory: python3 make_freefield.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def path_modes(n):
    u = np.arange(1, n + 1)
    j = np.arange(1, n + 1)
    phi = np.sqrt(2 / (n + 1)) * np.sin(np.pi * np.outer(j, u) / (n + 1))
    lam = 4 * np.sin(np.pi * j / (2 * (n + 1))) ** 2
    c = phi.sum(axis=1) ** 2
    return lam, c


def varpi_lattice(n, s, d=4):
    mu2 = (s / n) ** 2
    lam, c = path_modes(n)
    L = mu2
    C = 1.0
    for _ in range(d):
        L = np.add.outer(L, lam) if np.ndim(L) else L + lam
        C = np.multiply.outer(C, c) if np.ndim(C) else C * c
    return mu2 * n ** (-d) * np.sum(C / L)


def heat_content(t, kmax=4001):
    k = np.arange(1, kmax + 1, 2)[:, None]
    return (8 / np.pi ** 2) * np.sum(np.exp(-np.pi ** 2 * k ** 2 * t) / k ** 2, axis=0)


def varpi_limit(s, d=4):
    x = np.linspace(np.log(1e-9), np.log(60.0 / max(s * s, 1e-3) + 20), 6000)
    t = np.exp(x)
    f = s * s * np.exp(-s * s * t) * heat_content(t) ** d * t
    return np.trapezoid(f, x)


# sanity: n = 2 gives s^2/(s^2 + 4d)
assert abs(varpi_lattice(2, 1.0) - 1 / 17) < 1e-12
assert abs(varpi_limit(1.0) - 0.012954) < 2e-5, varpi_limit(1.0)

fig = plt.figure(figsize=(9.6, 3.9))
ax = fig.add_subplot(1, 2, 1, projection='3d')
n, s = 14, 2.0
mu2 = (s / n) ** 2
N = n * n
idx = lambda i, j: i * n + j
A = np.zeros((N, N))
for i in range(n):
    for j in range(n):
        A[idx(i, j), idx(i, j)] = 4 + mu2
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            a, b = i + di, j + dj
            if 0 <= a < n and 0 <= b < n:
                A[idx(i, j), idx(a, b)] = -1
w = mu2 * np.linalg.solve(A, np.ones(N)).reshape(n, n)
X, Y = np.meshgrid(np.arange(n), np.arange(n))
ax.plot_surface(X, Y, w, cmap='Blues', edgecolor=INK, linewidth=0.25, alpha=0.95)
ax.set_xlabel('$y_1$', labelpad=-6)
ax.set_ylabel('$y_2$', labelpad=-6)
ax.set_zlabel('response', labelpad=-2)
ax.tick_params(labelsize=7, pad=-2)
ax.view_init(elev=28, azim=-58)
ax.set_title('(a) one cube, $n=14$, $s=2$: mean %.3f' % w.mean(), fontsize=10)

ax2 = fig.add_subplot(1, 2, 2)
ss = np.linspace(0.05, 4, 120)
for nn, col in ((2, GOLD), (4, GREEN), (8, RED)):
    ax2.plot(ss, [varpi_lattice(nn, x) for x in ss], color=col, lw=1.4, label='exact, $n=%d$' % nn)
ax2.plot(ss, [varpi_limit(x) for x in ss], color=BLUE, lw=2.0, label=r'limit $n\to\infty$')
ax2.plot(ss, (5 / 6) ** 4 * ss ** 2 / (ss ** 2 + 40), color=INK, ls='--', lw=1.2, label=r'$(5/6)^4s^2/(s^2+40)$')
ax2.set_xlabel(r'$s=m\ell$')
ax2.set_ylabel(r'$\varpi$ in $d=4$')
ax2.set_xlim(0, 4)
ax2.set_ylim(0, None)
ax2.legend(frameon=False, fontsize=8, loc='upper left')
ax2.set_title('(b) the exact constant and the lower bound', fontsize=10)
fig.tight_layout()
save(fig, 'fig_freefield.png')
