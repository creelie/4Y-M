"""Figures for the stability, weak-coupling and mixing subsections.

fig_stability (ssec:rSstab): the bound s^2/(s^2+(1+t)16d/9), d = 4, on the
cube constant of the free field perturbed by the convex term (t/2)|grad phi|^2
(Theorem thm:rSnonmono), as a surface over (s, t) and along t for fixed s.
fig_maxwell (rD:sec): the bound 8 pi^2 (n+1)^2 / L^2 of Theorem rD:maxwell on
any constant of BP(n, varpi) for the Maxwell field.
fig_wilson (rDtwo:sec): the lower bound 1/(1200 beta^2) on the averaged
conditional variance of a Wilson loop at one link and the upper bound
min(4, 2241 r^8 / beta^2) on its variance (Theorem rDtwo:loops, 4 | L,
l -> 0).
fig_mixing (rW9:ssec): the constant varpi = Rb(1-x)^4 - (Rb-1)(1+x)^4,
x = s/n, of Theorem rW9:main for large n, with the level 1/3.

Run from this directory: python3 make_newsections.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def surface(ax, X, Y, Z, xl, yl, zl, cmap=cm.viridis):
    ax.plot_surface(X, Y, Z, cmap=cmap, linewidth=0, antialiased=True, alpha=0.95, rstride=1, cstride=1)
    ax.set_xlabel(xl, labelpad=4)
    ax.set_ylabel(yl, labelpad=4)
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel(zl, labelpad=6, rotation=90)
    ax.view_init(elev=24, azim=-128)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor((1, 1, 1, 0))
        a.pane.set_edgecolor(GREY)
    ax.tick_params(labelsize=8, pad=1)


def stability():
    d = 4
    fig = plt.figure(figsize=(9.2, 3.8))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    s = np.linspace(0.2, 8, 60)
    t = np.linspace(0, 20, 60)
    S, T = np.meshgrid(s, t)
    Z = S**2 / (S**2 + (1 + T) * 16 * d / 9)
    surface(ax, S, T, Z, r'$s=m\ell$', r'$t$', r'bound on $\varpi(\mu_t)$')
    ax.set_title('(a)', loc='left')
    bx = fig.add_subplot(1, 2, 2)
    tt = np.linspace(0, 40, 400)
    for sv, c in [(1, BLUE), (2, GREEN), (4, GOLD), (8, RED)]:
        bx.plot(tt, sv**2 / (sv**2 + (1 + tt) * 16 * d / 9), color=c, lw=1.6, label=rf'$s={sv}$')
    bx.set_xlabel(r'$t$')
    bx.set_ylabel(r'$s^2/(s^2+(1+t)\,16d/9)$')
    bx.set_ylim(0, 1.0)
    bx.legend(frameon=False)
    bx.set_title('(b)', loc='left')
    bx.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    save(fig, 'fig_stability.png')


def maxwell():
    fig = plt.figure(figsize=(9.2, 3.8))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    n = np.arange(2, 17)
    L = np.linspace(40, 400, 50)
    N, LL = np.meshgrid(n, L)
    Z = np.log10(8 * np.pi**2 * (N + 1)**2 / LL**2)
    surface(ax, N, LL, Z, r'$n$', r'$L$', r'$\log_{10}$ bound on $\varpi$', cmap=cm.plasma)
    ax.set_title('(a)', loc='left')
    bx = fig.add_subplot(1, 2, 2)
    Lv = np.logspace(1.3, 3.5, 300)
    for nv, c in [(2, BLUE), (4, GREEN), (8, GOLD), (16, RED)]:
        Lm = Lv[Lv >= nv + 2]
        bx.loglog(Lm, np.minimum(1, 8 * np.pi**2 * (nv + 1)**2 / Lm**2), color=c, lw=1.6, label=rf'$n={nv}$')
    bx.set_xlabel(r'$L$')
    bx.set_ylabel(r'$\min(1,8\pi^2(n+1)^2/L^2)$')
    bx.legend(frameon=False)
    bx.set_title('(b)', loc='left')
    bx.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    save(fig, 'fig_maxwell.png')


def wilson():
    fig = plt.figure(figsize=(4.6, 3.4))
    ax = fig.add_subplot(1, 1, 1)
    beta = np.logspace(np.log10(11), 3, 200)
    ax.loglog(beta, 1 / (1200 * beta**2), color=INK, lw=1.8, label=r'lower: $1/(1200\beta^2)$')
    for r, c in [(1, BLUE), (2, GREEN), (4, RED)]:
        ax.loglog(beta, np.minimum(4, 2241 * r**8 / beta**2), color=c, lw=1.6, ls='--',
                  label=rf'upper, $r={r}$')
    ax.set_xlabel(r'$\beta$')
    ax.set_ylabel('variance bounds')
    ax.set_ylim(1e-10, 30)
    ax.legend(frameon=False, fontsize=8, loc='lower left')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    save(fig, 'fig_wilson.png', width_px=1000)


def mixing():
    fig = plt.figure(figsize=(9.2, 3.8))
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    x = np.linspace(0, 0.2, 60)
    R = np.linspace(1, 1.6, 60)
    X, RR = np.meshgrid(x, R)
    Z = np.clip(RR * (1 - X)**4 - (RR - 1) * (1 + X)**4, 0, None)
    surface(ax, X, RR, Z, r'$s\,/\,n$', r'$\bar R$', r'$\varpi$', cmap=cm.cividis)
    ax.set_title('(a)', loc='left')
    bx = fig.add_subplot(1, 2, 2)
    xx = np.linspace(0, 0.2, 400)
    for Rv, c in [(1.0, BLUE), (1.1, GREEN), (1.25, GOLD), (1.5, RED)]:
        bx.plot(xx, np.clip(Rv * (1 - xx)**4 - (Rv - 1) * (1 + xx)**4, 0, None), color=c, lw=1.6, label=rf'$\bar R={Rv:g}$')
    bx.axhline(1 / 3, color=GREY, lw=1, ls=':')
    bx.text(0.198, 1 / 3 + 0.02, r'$\varpi=1/3$', ha='right', fontsize=9, color=INK)
    bx.axvline(1 / 8, color=GREY, lw=1, ls='--')
    bx.text(1 / 8 + 0.003, 0.9, r'$8s=n$', fontsize=9, color=INK)
    bx.set_xlabel(r'$s/n$')
    bx.set_ylabel(r'$\varpi$')
    bx.set_ylim(0, 1.02)
    bx.legend(frameon=False, loc='lower left')
    bx.set_title('(b)', loc='left')
    bx.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    save(fig, 'fig_mixing.png')


if __name__ == '__main__':
    stability()
    maxwell()
    wilson()
    mixing()
