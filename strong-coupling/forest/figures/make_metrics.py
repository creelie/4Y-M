"""Figure for the metrics across scales (ssec:bpmetrics).

fig_bpmetrics: (a) one column of the conditional-mean reconstruction
R = A^{-1} P^T (P A^{-1} P^T)^{-1} of the block averages P over blocks of
side 2 on the periodic 18 x 18 lattice, for A = m^2 - Delta with m^2 = 1/10,
drawn as bars over the fine sites: positive on its own block, negative on
the ring of neighbouring blocks (Proposition prop:bpnegative). (b) The
products prod_{i<k} theta_i, at least (sqrt 2)^k, against 1 + 2 zeta_k,
which stays below 4 + 2 sqrt 3, for the closed forms of Remark
rem:bpchain, on a logarithmic scale.

Run from this directory: python3 make_metrics.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def column(n=18, b=2, m2=0.1):
    sites = [(i, j) for i in range(n) for j in range(n)]
    idx = {s: a for a, s in enumerate(sites)}
    N = len(sites)
    A = np.eye(N) * (m2 + 4.0)
    for (i, j), a in idx.items():
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            A[a, idx[((i + di) % n, (j + dj) % n)]] -= 1.0
    M = n // b
    P = np.zeros((M * M, N))
    for (i, j), a in idx.items():
        P[(i // b) * M + (j // b), a] = 1.0 / (b * b)
    Ainv = np.linalg.inv(A)
    R = Ainv @ P.T @ np.linalg.inv(P @ Ainv @ P.T)
    w = (M // 2) * M + M // 2
    return n, R[:, w].reshape(n, n), Ainv.min()


def chain(K=20):
    rho = lambda k: np.sqrt((4 + 2 * 0.25 ** k) / 3)
    gam = lambda k: 1 + 2 / (rho(k) * (2 + rho(k)))
    tau = lambda k: 2 / (rho(k) * (2 + rho(k)))
    beta = [0.0, (np.sqrt(2) - 1) / 2, (tau(1) + gam(1) * (np.sqrt(2) - 1) / 2) / 2]
    for k in range(2, K):
        beta.append((1 / (2 * rho(k)) + (1 + 1 / rho(k)) * beta[k]) / 2)
    prod = np.cumprod([1.0] + [gam(i) for i in range(K)])
    return np.arange(K + 1), prod, 1 + 2 * np.array(beta[:K + 1])


def fig_bpmetrics():
    fig = plt.figure(figsize=(11, 4.4))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    n, C, amin = column()
    assert amin > 0 and C.min() < 0
    lo, hi = 3, n - 3
    X, Y = np.meshgrid(np.arange(lo, hi), np.arange(lo, hi), indexing='ij')
    Z = C[lo:hi, lo:hi]
    cols = np.where(Z.ravel() >= 0, BLUE, RED)
    ax.bar3d(X.ravel(), Y.ravel(), np.minimum(Z.ravel(), 0), 0.8, 0.8, np.abs(Z.ravel()),
             color=cols, shade=True, alpha=0.95, edgecolor='white', linewidth=0.2)
    ax.set_box_aspect((1, 1, 0.6))
    ax.view_init(elev=28, azim=-55)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_zlabel('entry of $R$', labelpad=6)
    ax.set_title('(a) a column of the conditional mean', fontsize=10)

    ax2 = fig.add_subplot(1, 2, 2)
    k, prod, theta = chain()
    ax2.semilogy(k, prod, 'o-', color=RED, ms=3.5, label=r'$\prod_{i<k}\theta_i$ (site-additive)')
    ax2.semilogy(k, np.sqrt(2.0) ** k, '--', color=GREY, label=r'$(\sqrt{2})^k$')
    ax2.semilogy(k, theta, 's-', color=BLUE, ms=3.5, label=r'$1+2\zeta_k$ (value and gradient)')
    ax2.axhline(4 + 2 * np.sqrt(3), color=GREEN, lw=1, ls=':', label=r'$4+2\sqrt{3}$')
    ax2.set_xlabel('number of scales $k$')
    ax2.set_xlim(0, k[-1])
    ax2.legend(frameon=False, fontsize=9, loc='upper left')
    ax2.set_title('(b) loss across scales', fontsize=10)
    fig.tight_layout()
    save(fig, 'fig_bpmetrics.png')


if __name__ == '__main__':
    fig_bpmetrics()
