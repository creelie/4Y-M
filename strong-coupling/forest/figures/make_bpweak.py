"""Figures for the averaged condition at strong coupling, at fixed cube size
and at weak coupling, and for the derivative in the ends.

fig_bpranges: (a) the bound of Theorem thm:bpfixed(ii),
log10(38400 n^4 (log beta + 6.11)/beta), as a surface over the cube side n
and log10 beta, with the plane varpi = 1; (b) the plane of the Wilson
coupling beta_W = 4 beta and the cube side n for SU(2): the strong-coupling
range of Corollary cor:bpranges, the region where Theorem thm:bpfixed
excludes BP(n, 1/2) on the torus of side 2n+4, and the continuum tuning
n ~ exp(3 pi^2 beta_W / 11), schematic (the power prefactor is omitted).

fig_bpsectors: (a) a three-dimensional section of the torus with the slices
{y_nu = t}; the centre twist multiplies the links that cross one slice by z,
a Polyakov loop crosses every slice once (charged), a local box misses a
slice (neutral) (Proposition prop:bpcentre); (b) the spectral measures of the
cube chain, schematic: a local function has no mass in (1 - epsilon, 1]
under BP_loc, while charged functions may carry mass near 1
(Theorem thm:bplocal, Corollary cor:bpcharged).

fig_endslinear: log10 of the two derivative constants of Theorem
thm:endlinear over beta_W and the tie length b_tie: L_* with the factor
e^{72 beta_W} (using A_tie < 1.5756 b_tie) and L_lin, linear in beta_W.

fig_bpneutral: (a) the great circle of SU(2) through I, u, -I, -u with the
holonomy flow sin(2 theta) d/d theta of Lemma lem:bpflow, which fixes the
centre and the traceless elements; (b) the Haar law of T = Tr P_0, density
sqrt(4 - T^2)/(2 pi), and the weight T^2 (4 - T^2) of Theorem thm:bpneutral,
which vanishes at T in {-2, 0, 2}.

Run from this directory: python3 make_bpweak.py."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from make_figures import save, clean3d

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
BLUE, RED, GREEN, INK, GREY, GOLD = '#3b6fb6', '#c0392b', '#2e8b57', '#222222', '#9aa5b1', '#d4a017'


def fig_bpranges():
    fig = plt.figure(figsize=(9.6, 4.2))
    ax = fig.add_axes([0.0, 0.02, 0.5, 0.96], projection='3d')
    n = np.linspace(2, 40, 60)
    lb = np.linspace(4, 16, 60)
    N, LB = np.meshgrid(n, lb)
    beta = 10.0 ** LB
    Z = np.log10(38400 * N ** 4 * (np.log(beta) + 6.11) / beta)
    cols = np.where(Z[..., None] < 0, np.array([0.75, 0.22, 0.17, 0.85]), np.array([0.23, 0.44, 0.71, 0.55]))
    ax.plot_surface(N, LB, Z, facecolors=cols, rstride=1, cstride=1, linewidth=0, antialiased=True, shade=False)
    P = np.zeros_like(Z)
    ax.plot_surface(N, LB, P, color=GREY, alpha=0.18, linewidth=0)
    nn = np.linspace(2, 40, 200)
    lbc = []
    for v in nn:
        f = lambda t: 38400 * v ** 4 * (np.log(10 ** t) + 6.11) / 10 ** t - 1
        lo, hi = 4.0, 16.0
        if f(lo) < 0 or f(hi) > 0:
            lbc.append(np.nan); continue
        for _ in range(60):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
        lbc.append(lo)
    ax.plot(nn, lbc, 0 * nn, color=INK, lw=1.4)
    ax.set_xlabel('cube side $n$', labelpad=2)
    ax.set_ylabel(r'$\log_{10}\beta$', labelpad=2)
    ax.set_zlabel(r'$\log_{10}$ bound on $\varpi$', labelpad=2)
    ax.set_box_aspect((1, 1, 0.7))
    ax.view_init(elev=24, azim=-128)
    ax.tick_params(labelsize=7, pad=0)
    fig.text(0.03, 0.92, '(a)', fontsize=11)
    fig.text(0.08, 0.86, r'red: the bound is below $\varpi=1$', color=RED, fontsize=8)

    bx = fig.add_axes([0.6, 0.15, 0.37, 0.75])
    bw = np.logspace(-1.5, 14.0, 800)
    b = bw / 4.0
    # Theorem thm:bpfixed excludes BP(n, 1/2) when 38400 n^4 (log b + 6.11)/b < 1/2, b >= 1
    with np.errstate(invalid='ignore'):
        nobs = np.where(b >= 1, (b / (2 * 38400 * (np.log(b) + 6.11))) ** 0.25, np.nan)
    bx.fill_between(bw, 2, np.where(nobs > 2, nobs, 2), where=nobs > 2, color=RED, alpha=0.25, lw=0)
    bx.plot(bw, nobs, color=RED, lw=1.2)
    bx.axvspan(bw[0], 2 / 9, color=GREEN, alpha=0.18, lw=0)
    tun = 2.0 * np.exp(np.minimum(3 * np.pi ** 2 * bw / 11, 50.0))
    bx.plot(bw, tun, color=BLUE, lw=1.3, ls='--')
    bx.set_xscale('log'); bx.set_yscale('log')
    bx.set_xlim(bw[0], bw[-1]); bx.set_ylim(2, 1e4)
    bx.set_xlabel(r'$\beta_W$ for $\mathrm{SU}(2)$'); bx.set_ylabel('cube side $n$')
    bx.text(0.075, 60, 'BP proved', color=GREEN, fontsize=8, rotation=90, ha='center')
    bx.text(1e10, 3.0, 'BP$(n,1/2)$ fails\non $\\Lambda_{2n+4}$', color=RED, fontsize=8)
    bx.text(12.0, 1.5e3, 'tuning', color=BLUE, fontsize=8)
    bx.text(1e5, 400, 'open', color=INK, fontsize=9)
    for sp in ('top', 'right'):
        bx.spines[sp].set_visible(False)
    fig.text(0.55, 0.92, '(b)', fontsize=11)
    save(fig, 'fig_bpranges.png')


def fig_bpsectors():
    fig = plt.figure(figsize=(9.6, 4.2))
    ax = fig.add_axes([0.0, 0.0, 0.52, 1.0], projection='3d')
    clean3d(ax)
    ax.computed_zorder = False
    L = 6
    # slices {y_0 = t}, drawn as planes x = t
    for t in range(L + 1):
        v = [(t, 0, 0), (t, L, 0), (t, L, L), (t, 0, L)]
        col = GOLD if t == 2 else GREY
        al = 0.32 if t == 2 else 0.07
        ax.add_collection3d(Poly3DCollection([v], facecolor=col, edgecolor=GREY, lw=0.4, alpha=al))
    # links crossing slice t=2 in direction 0: from x=2 to x=3 are multiplied by z
    segs = [[(2, y, z), (3, y, z)] for y in range(L + 1) for z in range(L + 1)]
    ax.add_collection3d(Line3DCollection(segs, colors=GOLD, lw=1.0))
    # Polyakov loop in direction 0
    ax.plot([0, L], [4.5, 4.5], [4.5, 4.5], color=RED, lw=2.6)
    ax.text(L + 0.2, 4.5, 4.9, r'Polyakov loop: $\mathcal{A}_{e_0}$', color=RED, fontsize=9)
    # local box avoiding slice t=2
    b0, s = 3.5, 1.6
    cub = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]) * s + np.array([b0, 0.8, 0.8])
    faces = [[cub[i] for i in f] for f in ([0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4])]
    ax.add_collection3d(Poly3DCollection(faces, facecolor=GREEN, edgecolor=INK, lw=0.5, alpha=0.45))
    ax.text(b0 + s + 0.3, 0.4, 0.5, r'local box: $\mathcal{A}_0$', color=GREEN, fontsize=9)
    ax.text(-1.0, 0.0, -1.6, r'twist $Z^0_t$: links across $\{y_0=t\}$ times $z$', color=GOLD, fontsize=9)
    ax.set_xlim(0, L); ax.set_ylim(0, L); ax.set_zlim(0, L)
    ax.view_init(elev=20, azim=-62)
    fig.text(0.03, 0.92, '(a)', fontsize=11)

    bx = fig.add_axes([0.6, 0.18, 0.37, 0.68])
    lam = np.linspace(0, 1, 500)
    eps = 0.18
    loc = np.exp(-((lam - 0.42) / 0.16) ** 2) * (lam <= 1 - eps)
    chg = 0.55 * np.exp(-((lam - 0.97) / 0.02) ** 2) + 0.4 * np.exp(-((lam - 0.5) / 0.2) ** 2)
    bx.fill_between(lam, 0, loc, color=GREEN, alpha=0.45, lw=0, label='local, neutral')
    bx.plot(lam, chg, color=RED, lw=1.3, label='charged')
    bx.axvspan(1 - eps, 1, color=GREY, alpha=0.22, lw=0)
    bx.text(1 - eps + 0.01, 1.05, r'$(1-\epsilon,1]$', fontsize=8)
    bx.set_xlim(0, 1.02); bx.set_ylim(0, 1.2)
    bx.set_xlabel(r'spectrum of $\mathcal{K}$'); bx.set_yticks([])
    bx.legend(loc='upper left', fontsize=8, frameon=False)
    for sp in ('top', 'right', 'left'):
        bx.spines[sp].set_visible(False)
    fig.text(0.55, 0.92, '(b)', fontsize=11)
    save(fig, 'fig_bpsectors.png')


def fig_endslinear():
    fig = plt.figure(figsize=(7.4, 4.6))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], projection='3d')
    bw = np.linspace(0.0, 0.6, 50)
    bt = np.linspace(0.0, 20.0, 50)
    B, T = np.meshgrid(bw, bt)
    Ls = np.log10(128.07 * B + (180.4 + 50.6 * 1.5756 * T) * np.exp(72 * B))
    Ll = np.log10(7192 + 443620 * B + T * (79.69 + 27216 * B))
    ax.plot_surface(B, T, Ls, color=RED, alpha=0.8, linewidth=0)
    ax.plot_surface(B, T, Ll, color=BLUE, alpha=0.8, linewidth=0)
    ax.set_xlabel(r'$\beta_W$', labelpad=2)
    ax.set_ylabel(r'$\mathfrak{b}_{\mathrm{tie}}$', labelpad=2)
    ax.set_zlabel(r'$\log_{10}L$', labelpad=2)
    ax.text(0.6, 0, 20, r'$L_*$ (factor $e^{72\beta_W}$)', color=RED, fontsize=9)
    ax.text(0.6, 20, 4.2, r'$L_{\mathrm{lin}}$', color=BLUE, fontsize=9)
    ax.set_box_aspect((1, 1, 0.75))
    ax.view_init(elev=18, azim=-50)
    ax.tick_params(labelsize=7, pad=0)
    save(fig, 'fig_endslinear.png')


def fig_bpneutral():
    fig = plt.figure(figsize=(9.6, 4.0))
    ax = fig.add_axes([0.02, 0.06, 0.42, 0.88])
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(th), np.sin(th), color=GREY, lw=1.2)
    for t in np.linspace(0, 2 * np.pi, 24, endpoint=False) + np.pi / 24:
        v = np.sin(2 * t) * 0.32
        x, y = np.cos(t), np.sin(t)
        ax.annotate('', xy=(x - v * np.sin(t), y + v * np.cos(t)), xytext=(x, y),
                    arrowprops=dict(arrowstyle='->', color=BLUE, lw=1.1))
    for (x, y, lab, c) in ((1, 0, r'$I$', RED), (-1, 0, r'$-I$', RED), (0, 1, r'$u$', GOLD), (0, -1, r'$-u$', GOLD)):
        ax.plot(x, y, 'o', color=c, ms=7, zorder=5)
        ax.text(1.22 * x - (0.06 if x < 0 else 0), 1.18 * y - 0.04, lab, color=c, fontsize=11, ha='center')
    ax.text(0, 0.05, r'$\cos\theta\,I+\sin\theta\,u$', ha='center', fontsize=10)
    ax.text(0, -0.2, r'velocity $\sin 2\theta\,\partial_\theta$', ha='center', fontsize=9, color=BLUE)
    ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.4, 1.4); ax.set_aspect('equal'); ax.axis('off')
    fig.text(0.03, 0.92, '(a)', fontsize=11)

    bx = fig.add_axes([0.56, 0.16, 0.41, 0.72])
    T = np.linspace(-2, 2, 600)
    dens = np.sqrt(np.clip(4 - T ** 2, 0, None)) / (2 * np.pi)
    w = T ** 2 * (4 - T ** 2)
    bx.fill_between(T, 0, dens, color=GREY, alpha=0.35, lw=0, label=r'Haar law of $T$')
    bx.plot(T, w / 16, color=GREEN, lw=1.6, label=r'$T^2(4-T^2)/16$')
    for t0 in (-2, 0, 2):
        bx.plot(t0, 0, 'o', color=RED, ms=6, clip_on=False, zorder=5)
    bx.text(0.0, 0.045, 'exceptional set', color=RED, fontsize=8, ha='center')
    bx.set_xlim(-2.1, 2.1); bx.set_ylim(0, 0.42)
    bx.set_xlabel(r'$T=\mathrm{Tr}\,P_0$'); bx.set_yticks([])
    bx.legend(loc='upper right', fontsize=8, frameon=False)
    for sp in ('top', 'right', 'left'):
        bx.spines[sp].set_visible(False)
    fig.text(0.53, 0.92, '(b)', fontsize=11)
    save(fig, 'fig_bpneutral.png')


if __name__ == '__main__':
    fig_bpranges()
    fig_bpsectors()
    fig_endslinear()
    fig_bpneutral()
