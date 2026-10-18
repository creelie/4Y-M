"""Figures of the PRD paper: python3 make_figures.py writes figures/fig_*.png and .pdf."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from mpl_toolkits.mplot3d import proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.special import iv

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "axes.titlesize": 9.5, "axes.labelsize": 9, "legend.fontsize": 8,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "savefig.dpi": 200,
})
BLUE, RED, ORANGE, GREEN, GREY = "#1f5aa6", "#c0392b", "#e67e22", "#2e8b57", "#7f8c8d"
W = 7.0


class Arrow3D(FancyArrowPatch):
    def __init__(self, p, q, **kw):
        super().__init__((0, 0), (0, 0), **kw)
        self._p, self._q = p, q

    def do_3d_projection(self, renderer=None):
        xs, ys, zs = proj3d.proj_transform(*zip(self._p, self._q), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        return np.min(zs)


def arrow3(ax, p, q, color="k", lw=1.2, ms=9):
    ax.add_artist(Arrow3D(p, q, mutation_scale=ms, lw=lw, arrowstyle="-|>", color=color))


def clean3d(ax):
    ax.set_axis_off()
    ax.set_box_aspect((1, 1, 1))


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(HERE, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


# 1. thresholds -------------------------------------------------------------
def fig_thresholds():
    fig, ax = plt.subplots(figsize=(W, 2.5))
    rows = [
        ("earlier curvature bound", 1 / 12, GREY),
        ("earlier single-link bound", 1 / 6, GREY),
        ("Theorems B, D (curvature)", 1 / 8, BLUE),
        ("Theorem E (chordal Dobrushin)", 1 / (3 * np.sqrt(3)), BLUE),
        ("ceiling of single-link arguments", 2 / 9, RED),
    ]
    for i, (lab, v, c) in enumerate(rows):
        y = len(rows) - i
        ax.barh(y, v - 0.04, left=0.04, height=0.55, color=c, alpha=0.85 if c != RED else 0.35,
                hatch=None if c != RED else "///", edgecolor=c)
        ax.text(v * 1.06, y, f"{lab}: {v:.4f}" if c != RED else f"{lab}: $2/9$", va="center", fontsize=8)
    ax.axvspan(1.8, 2.6, color=ORANGE, alpha=0.3)
    ax.text(2.16, 3.0, "Monte Carlo\ncrossover", ha="center", va="center", fontsize=8, color="#8a4b0f")
    ax.set_xscale("log")
    ax.set_xlim(0.04, 4)
    ax.set_ylim(0.4, len(rows) + 0.6)
    ax.set_yticks([])
    ax.set_xlabel(r"Wilson coupling $\beta_W$ (logarithmic)")
    ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1, 2])
    ax.set_xticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2"])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    save(fig, "fig_thresholds")


# 2. plaquette ---------------------------------------------------------------
C35, S35 = np.cos(np.radians(40)), np.sin(np.radians(40))


def cab(p):
    """Cabinet projection: x to the right, y into the page at 35 degrees, z up."""
    x, y, z = p
    return np.array([x + 0.75 * y * C35, z + 0.75 * y * S35])


def arr2(ax, p, q, color="k", lw=1.4, ms=10, shrink=0.0):
    ax.add_patch(FancyArrowPatch(tuple(p), tuple(q), arrowstyle="-|>", mutation_scale=ms, lw=lw,
                                 color=color, shrinkA=shrink, shrinkB=shrink))


def fig_plaquette():
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.6))
    ax = axs[0]
    L = 2.0
    P3 = [(0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0)]
    P = [cab(p) for p in P3]
    ax.add_patch(plt.Polygon(P, closed=True, fc="#dbe7f5", ec="none"))
    for i in range(4):
        a, b = P[i], P[(i + 1) % 4]
        arr2(ax, a + 0.06 * (b - a), b - 0.06 * (b - a), lw=1.8)
    off = [(0, -0.22), (0.25, -0.05), (0, 0.2), (0.3, -0.12)]
    for i, o in enumerate(off):
        m = (P[i] + P[(i + 1) % 4]) / 2 + np.array(o)
        ax.text(*m, rf"$e_{i+1}$", ha="center", va="center", fontsize=10)
    sg = ["+", "$-$", "+", "$-$"]
    so = [(-0.22, -0.12), (0.2, -0.12), (0.2, 0.12), (-0.2, 0.12)]
    for i in range(4):
        ax.text(*(P[i] + np.array(so[i])), sg[i], ha="center", va="center", fontsize=12, color="#555")
    ax.text(*(P[0] + np.array([0.12, -0.22])), "$x$", fontsize=11)
    cols = [BLUE, GREEN, ORANGE, RED]
    for i, ang in enumerate([98, 116, 134, 152]):
        e = P[0] + 1.25 * np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
        arr2(ax, P[0], e, color=cols[i], lw=1.5)
        t = P[0] + 1.45 * np.array([np.cos(np.radians(ang)), np.sin(np.radians(ang))])
        ax.text(*t, rf"$Y_{i+1}$", color=cols[i], ha="center", va="center", fontsize=10)
    ax.set_xlim(-1.0, 3.6); ax.set_ylim(-0.5, 2.3)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("link vectors at $x$ and corner signs")

    ax = axs[1]
    shades = ["#dbe7f5", "#e3f1e4", "#fbe7d3", "#f6dcdc"]
    quads = [((0, 0), (1, 0), (1, 1), (0, 1)), ((0, 0), (0, 1), (-1, 1), (-1, 0)),
             ((0, 0), (-1, 0), (-1, -1), (0, -1)), ((0, 0), (0, -1), (1, -1), (1, 0))]
    for q, c in zip(quads, shades):
        ax.add_patch(plt.Polygon([cab((1.4 * a, 1.4 * b, 0)) for a, b in q], closed=True, fc=c, ec="#555", lw=0.8))
    for d in [(1, 0), (0, 1), (-1, 0), (0, -1)]:
        arr2(ax, cab((0, 0, 0)), cab((0.8 * d[0], 0.8 * d[1], 0)), color=BLUE, lw=2.0)
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            ax.plot(*cab((1.4 * a, 1.4 * b, 0)), "k.", ms=4)
    ax.text(*(cab((0, 0, 0)) + np.array([0.1, -0.3])), "$x$", fontsize=11)
    ax.annotate("$|D_xX|^2$: one corner\nfrom each plaquette", xy=cab((0.25, 0.25, 0)),
                xytext=(1.2, 1.55), fontsize=9, color=BLUE,
                arrowprops=dict(arrowstyle="-", color=BLUE, lw=0.7))
    ax.set_xlim(-2.3, 3.2); ax.set_ylim(-1.0, 2.1)
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("four plaquettes of one plane at $x$")
    fig.subplots_adjust(wspace=0.05)
    save(fig, "fig_plaquette")


# 3. extremal fields ---------------------------------------------------------
def fig_extremal():
    fig = plt.figure(figsize=(W, 3.1))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    eps = np.array([1.0, -0.4, -0.6])
    dirs = np.eye(3)
    for x in range(3):
        for y in range(3):
            for z in range(3):
                ax.plot([x], [y], [z], "o", color="#333", ms=2.2)
                s = (-1) ** (x + y + z)
                for mu in range(3):
                    p = np.array([x, y, z], float)
                    if p[mu] + 1 > 2:
                        continue
                    ax.plot(*zip(p, p + dirs[mu]), color="#bbb", lw=0.6)
                    mid = p + 0.5 * dirs[mu]
                    v = s * eps[mu] * 0.5 * np.array([0, 0, 1.0])
                    c = BLUE if s * eps[mu] > 0 else RED
                    arrow3(ax, mid - v / 2, mid + v / 2, color=c, lw=1.1, ms=6)
    ax.view_init(elev=22, azim=-58)
    clean3d(ax)
    ax.set_title(r"identity field: $X_{x,\mu}=(-1)^{|x|}\varepsilon_\mu T$", pad=-2, fontsize=9)

    ax = fig.add_subplot(1, 2, 2, projection="3d")
    for x in range(3):
        for y in range(3):
            for z in range(3):
                ax.plot([x], [y], [z], "o", color="#333", ms=2.2)
                p = np.array([x, y, z])
                for mu in range(3):
                    if p[mu] + 1 > 2:
                        continue
                    q = p + np.eye(3, dtype=int)[mu]
                    red = (sum(p[:mu]) % 2) == 1
                    ax.plot(*zip(p, q), color=RED if red else "#999", lw=1.8 if red else 0.7)
    for x in range(2):
        for y in range(2):
            ax.add_collection3d(Poly3DCollection([[(x, y, 0), (x + 1, y, 0), (x + 1, y + 1, 0), (x, y + 1, 0)]],
                                                 facecolor="#f6d6d2", edgecolor="none", alpha=0.55))
            ax.add_collection3d(Poly3DCollection([[(0, x, y), (0, x + 1, y), (0, x + 1, y + 1), (0, x, y + 1)]],
                                                 facecolor="#f6d6d2", edgecolor="none", alpha=0.35))
    ax.view_init(elev=22, azim=-58)
    clean3d(ax)
    ax.set_title(r"frustrated field: red links carry $z$", pad=-2, fontsize=9)
    fig.subplots_adjust(wspace=0.12)
    save(fig, "fig_extremal")


# 4. spectrum ----------------------------------------------------------------
def fig_spectrum():
    d = DATA
    rnd = np.loadtxt(os.path.join(d, "spec_random.dat"))
    idn = np.loadtxt(os.path.join(d, "spec_identity.dat"))
    fru = np.loadtxt(os.path.join(d, "spec_frustrated.dat"))
    fig = plt.figure(figsize=(W, 2.9))
    ax = fig.add_subplot(1, 2, 1)
    ax.hist(rnd, bins=np.arange(-17, 17.01, 0.5), color="#9db7d8", edgecolor="white", lw=0.3,
            label="two random fields")
    for arr, c, lab in ((idn, BLUE, "identity"), (fru, RED, "frustrated")):
        u, k = np.unique(np.round(arr), return_counts=True)
        m = u != 0
        ax.bar(u[m], k[m], width=0.45, color=c, alpha=0.85, label=lab)
    for s in (-1, 1):
        ax.axvline(12 * s, color="k", ls="--", lw=0.9)
        ax.axvline(16 * s, color="k", ls=":", lw=0.9)
    for x0, lab in ((12, r"$4d=12$"), (16, r"$8(d-1)=16$"), (-12, r"$-4d$"), (-16, r"$-8(d-1)$")):
        ax.text(x0 - 0.5, 75, lab, fontsize=7.5, ha="right", va="center", rotation=90)
    ax.set_xlim(-19, 27); ax.set_ylim(0, 128)
    ax.set_xlabel(r"eigenvalue of $\mathrm{Hess}\,S_0$"); ax.set_ylabel("count")
    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1.03, 1.0), fontsize=7.5)
    ax.set_title("spectra on $(\\mathbb{Z}/4\\mathbb{Z})^3$, zero modes omitted")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax = fig.add_subplot(1, 2, 2, projection="3d")
    bins = np.arange(-14, 14.01, 1.0)
    for j, (arr, c) in enumerate(((rnd, "#9db7d8"), (idn, BLUE), (fru, RED))):
        h, _ = np.histogram(arr[np.abs(arr) > 1e-6], bins=bins)
        if j == 0:
            h = h / 2
        ax.bar(bins[:-1] + 0.5, h, zs=j, zdir="y", width=0.8, color=c, alpha=0.85, edgecolor="none")
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels(["random", "identity", "frustrated"], fontsize=7)
    ax.set_xticks([-12, 0, 12])
    ax.set_xlabel("eigenvalue", labelpad=2); ax.set_zlabel("count per field", labelpad=2)
    ax.view_init(elev=24, azim=-60)
    ax.set_box_aspect((1.6, 1, 0.8))
    ax.set_title("the three spectra side by side", pad=-2)
    fig.subplots_adjust(wspace=0.12)
    save(fig, "fig_spectrum")


# 5. weights -----------------------------------------------------------------
def fig_weights():
    a = 1.6
    fig = plt.figure(figsize=(W, 3.0))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    xs, ys, hs, cs = [], [], [], []
    cmap = plt.get_cmap("viridis")
    for x in range(-4, 4):
        for y in range(-4, 5):
            n = abs(x) + abs(y)
            xs.append(x + 0.5); ys.append(y); hs.append(a ** (-n)); cs.append(cmap(n / 8))
    ax.bar3d(np.array(xs) - 0.3, np.array(ys) - 0.15, 0, 0.6, 0.3, hs, color=cs, shade=True)
    ax.set_xlabel("$x_1$", labelpad=-6); ax.set_ylabel("$x_2$", labelpad=-6); ax.set_zlabel("")
    ax.tick_params(pad=-2)
    ax.view_init(elev=28, azim=-52)
    ax.set_title(r"$w_e=a^{-|e|}$ on horizontal links, $a=1.6$", pad=-2)

    ax = fig.add_subplot(1, 2, 2)
    aa = np.linspace(1, 3, 300)
    price = (np.sqrt(aa) - 1) ** 2 / np.sqrt(aa)
    bw = 4.0 / (32.0 + 18.0 * price)
    ax.plot(aa, bw, color=BLUE, lw=1.7)
    ax.axhline(1 / 8, color=GREY, ls=":", lw=1)
    ax.text(2.2, 0.1255, r"$\beta_W=1/8$ at $a=1$", color=GREY, fontsize=8, va="bottom")
    ax.plot([a], [4.0 / (32.0 + 18.0 * (np.sqrt(a) - 1) ** 2 / np.sqrt(a))], "o", color=ORANGE, ms=4)
    ax.text(a + 0.05, 0.1195, "$a=1.6$", color=ORANGE, fontsize=8)
    ax.set_xlabel("weight base $a$"); ax.set_ylabel(r"largest $\beta_W$ with $K_a>0$")
    ax.set_ylim(0.105, 0.128)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.set_title(r"Theorem D for $\mathrm{SU}(2)$, $d=4$")
    fig.subplots_adjust(wspace=0.38)
    save(fig, "fig_weights")


# 6. sphere ------------------------------------------------------------------
def fig_sphere():
    fig = plt.figure(figsize=(W, 3.2))
    u, v = np.mgrid[0:2 * np.pi:120j, 0:np.pi:60j]
    X, Y, Z = np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v)
    kappa = 2.0
    dens = np.exp(kappa * Z)
    dens = (dens - dens.min()) / (dens.max() - dens.min())
    for k, (h, title) in enumerate(((np.array([1.0, 0, 0]), r"$h\perp b$: the reflection preserves $\nu_b$"),
                                    (np.array([0.6, 0, 0.8]), r"$h$ not orthogonal to $b$: the open case"))):
        ax = fig.add_subplot(1, 2, k + 1, projection="3d")
        ax.plot_surface(X, Y, Z, facecolors=plt.cm.YlOrRd(0.15 + 0.8 * dens), rstride=2, cstride=2,
                        lw=0, antialiased=True, alpha=0.55, shade=False)
        arrow3(ax, (0, 0, 0), (0, 0, 1.45), color="k", lw=1.6)
        ax.text(0.05, 0.05, 1.55, "$b$", fontsize=11)
        # plane h^perp
        e1 = np.cross(h, [0, 1.0, 0]); e1 /= np.linalg.norm(e1)
        e2 = np.cross(h, e1)
        s, t = np.meshgrid(np.linspace(-1.2, 1.2, 2), np.linspace(-1.2, 1.2, 2))
        P = s[..., None] * e1 + t[..., None] * e2
        ax.plot_surface(P[..., 0], P[..., 1], P[..., 2], color="#6a8fc7", alpha=0.18, lw=0)
        arrow3(ax, (0, 0, 0), tuple(1.25 * h), color=BLUE, lw=1.4)
        ax.text(*(1.38 * h), "$h$", color=BLUE, fontsize=11)
        x = np.array([0.55, -0.35, 0.76]); x /= np.linalg.norm(x)
        Rx = x - 2 * np.dot(h, x) * h
        ax.plot(*x, "o", color="k", ms=4); ax.plot(*Rx, "o", color="k", ms=4, mfc="white")
        ax.plot(*zip(x, Rx), color="k", lw=0.8, ls="--")
        ax.text(*(x * 1.12), "$x$", fontsize=10)
        ax.text(*(Rx * 1.15), "$Rx$", fontsize=10)
        ax.set_xlim(-1.1, 1.1); ax.set_ylim(-1.1, 1.1); ax.set_zlim(-1.1, 1.1)
        ax.view_init(elev=14, azim=-40)
        clean3d(ax)
        ax.set_title(title, pad=-6)
    fig.subplots_adjust(wspace=0.0)
    save(fig, "fig_sphere")


# 7. surface -----------------------------------------------------------------
def fig_surface():
    raw = [l.split() for l in open(os.path.join(DATA, "surf.dat")) if l.strip()]
    D = np.array(raw, float)
    ks = np.unique(D[:, 0]); ps = np.unique(D[:, 1])
    Lz = D[:, 2].reshape(len(ks), len(ps))
    K, P = np.meshgrid(ks, ps, indexing="ij")
    fig = plt.figure(figsize=(W, 3.3))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.plot_surface(K, P, Lz, cmap="viridis", edgecolor="#333", lw=0.25, alpha=0.95)
    ax.plot_surface(K, P, 0.25 + 0 * K, color="#999", alpha=0.25, lw=0)
    kk = np.linspace(0, 2, 60)
    ax.plot(kk, np.pi / 2 + 0.0 * kk, (3 * (1 + np.sqrt(kk ** 2 + 9))) ** -0.5, color=BLUE, lw=1.6)
    A = iv(2, kk[1:]) / iv(1, kk[1:])
    ax.plot(kk[1:], np.pi / 2 + 0 * kk[1:], A / kk[1:], color=RED, lw=1.6)
    ax.set_xlabel(r"$\kappa$", labelpad=-4); ax.set_ylabel(r"$\varphi$", labelpad=-4)
    ax.set_zlabel(""); ax.tick_params(pad=-2)
    ax.set_zlim(0.19, 0.30)
    ax.view_init(elev=20, azim=-50)
    ax.set_title(r"rate surface $L(\kappa,\varphi)$ and the plane $1/4$", pad=-2)

    ax = fig.add_subplot(1, 2, 2)
    for j, ph in enumerate(ps):
        ax.plot(ks, Lz[:, j], color=plt.cm.viridis(j / (len(ps) - 1)), lw=1.1, marker="o", ms=2.2)
    ax.plot(kk, (3 * (1 + np.sqrt(kk ** 2 + 9))) ** -0.5, color=BLUE, lw=1.4, label="proved bound")
    ax.plot(kk[1:], A / kk[1:], color=RED, lw=1.4, ls="--", label=r"$A(\kappa)/\kappa$ (edge $\varphi=\pi/2$)")
    ax.axhline(0.25, color=GREY, lw=1)
    ax.text(1.35, 0.2515, "conjectured $1/4$", color=GREY, fontsize=8)
    ax.set_xlabel(r"$\kappa$"); ax.set_ylabel(r"$L(\kappa,\varphi)$")
    ax.set_ylim(0.20, 0.295)
    ax.legend(frameon=False, loc="upper right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(r"sections at $\varphi=0,\pi/12,\dots,\pi/2$")
    fig.subplots_adjust(wspace=0.36)
    save(fig, "fig_surface")


# 8. box ---------------------------------------------------------------------
def cube_edges(ax, o, n, color, lw, alpha=1.0, ls="-"):
    o = np.array(o, float)
    for a in range(3):
        for b in (0, 1):
            for c in (0, 1):
                p = o.copy(); q = o.copy()
                idx = [i for i in range(3) if i != a]
                p[idx[0]] += b * n; p[idx[1]] += c * n
                q[:] = p; q[a] += n
                ax.plot(*zip(p, q), color=color, lw=lw, alpha=alpha, ls=ls)


def fig_box():
    fig = plt.figure(figsize=(W, 3.2))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    cube_edges(ax, (-1, -1, -1), 5, GREY, 0.8, ls="--")
    for i in range(4):
        for j in range(4):
            ax.plot([0, 3], [i, i], [j, j], color="#9db7d8", lw=0.6)
            ax.plot([i, i], [0, 3], [j, j], color="#9db7d8", lw=0.6)
            ax.plot([i, i], [j, j], [0, 3], color="#9db7d8", lw=0.6)
    cube_edges(ax, (0, 0, 0), 3, BLUE, 1.4)
    ax.plot([-1, 0], [1, 1], [2, 2], color=RED, lw=3)
    ax.text(-1.4, 1, 2.4, "$y$", color=RED, fontsize=11)
    ax.plot([1, 2], [2, 2], [1, 1], color=ORANGE, lw=3)
    ax.text(1.5, 2.0, 1.4, "$b$", color=ORANGE, fontsize=11)
    ax.text(4.3, -1.0, -1.9, "frozen shell", color=GREY, fontsize=8)
    ax.text(3.0, 3.2, 3.2, "$V$", color=BLUE, fontsize=11)
    ax.view_init(elev=20, azim=-62)
    clean3d(ax)
    ax.set_title(r"$\kappa_V(y,b)$: disagreement at $b$ caused by $y$", pad=-4)

    ax = fig.add_subplot(1, 2, 2, projection="3d")
    shifts = [(0, 0, 0), (-1, 0, 0), (0, -1, 0), (-1, -1, 0), (0, 0, -1)]
    cols = [BLUE, GREEN, "#8e44ad", ORANGE, GREY]
    for s, c in zip(shifts, cols):
        o = np.array(s, float) + 0.04 * np.array(s)
        cube_edges(ax, o, 2, c, 1.1, alpha=0.9)
    ax.plot([0.5, 1.5], [1, 1], [1, 1], color=RED, lw=3.2)
    ax.text(0.15, 1.0, 1.25, "$b$", color=RED, fontsize=11)
    ax.view_init(elev=22, azim=-58)
    clean3d(ax)
    ax.set_title("translates of $V$ containing $b$", pad=-4)
    fig.subplots_adjust(wspace=0.0)
    save(fig, "fig_box")


# 9. mass gap ----------------------------------------------------------------
def fig_massgap():
    fig = plt.figure(figsize=(W, 2.8))
    ax = fig.add_subplot(1, 2, 1)
    b1 = np.linspace(0.004, 1 / (3 * np.sqrt(3)) - 1e-4, 400)
    b2 = np.linspace(0.004, 2 / 9 - 1e-4, 400)
    ax.plot(b1, -np.log(3 * np.sqrt(3) * b1), color=BLUE, lw=1.7, label=r"Theorem F: $-\log(3\sqrt{3}\,\beta_W)$")
    ax.plot(b2, -np.log(4.5 * b2), color=RED, lw=1.4, ls="--", label=r"if the quarter holds: $-\log(\frac{9}{2}\beta_W)$")
    ax.axvline(1 / (3 * np.sqrt(3)), color=BLUE, lw=0.7, ls=":")
    ax.axvline(2 / 9, color=RED, lw=0.7, ls=":")
    ax.set_ylim(0, 5.5); ax.set_xlim(0, 0.24)
    ax.set_xlabel(r"$\beta_W$"); ax.set_ylabel("lattice mass gap $m$ (lower bound)")
    ax.legend(frameon=False, loc="upper right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("lower bounds for the gap in lattice units")

    ax = fig.add_subplot(1, 2, 2, projection="3d")
    bb = np.linspace(0.01, 1 / (3 * np.sqrt(3)) - 2e-3, 50)
    tt = np.linspace(0, 6, 40)
    B, T = np.meshgrid(bb, tt)
    Mz = np.exp(np.log(3 * np.sqrt(3) * B) * T)
    ax.plot_surface(B, T, Mz, cmap="Blues_r", lw=0, alpha=0.95)
    ax.set_xlabel(r"$\beta_W$", labelpad=-4); ax.set_ylabel("distance $t$", labelpad=-4)
    ax.set_zlabel("correlation bound", labelpad=-4); ax.tick_params(pad=-2)
    ax.view_init(elev=24, azim=-130)
    ax.set_title(r"decay $e^{-mt}=(3\sqrt{3}\,\beta_W)^t$", pad=-2)
    fig.subplots_adjust(wspace=0.2)
    save(fig, "fig_massgap")


# 10. scales ------------------------------------------------------------------
def fig_scales():
    fig, ax = plt.subplots(figsize=(W, 2.9))
    xc = 1.0

    def g2(x):
        sp = np.log1p(np.exp(1.2 * (xc - x))) / 1.2
        return 1.0 / (0.08 * sp + 0.125)

    xs = np.linspace(-16, 6, 2000)
    lo, hi = 1.6, 4.5
    xa = xs[np.argmin(np.abs(g2(xs) - lo))]
    xb = xs[np.argmin(np.abs(g2(xs) - hi))]
    ax.axvspan(-16, xa, color=GREEN, alpha=0.16, lw=0)
    ax.axvspan(xa, xb, color=ORANGE, alpha=0.28, lw=0)
    ax.axvspan(xb, 6, color=BLUE, alpha=0.16, lw=0)
    ax.text((-16 + xa) / 2, 8.6, "weak coupling: Ba\u0142aban's\nultraviolet stability", color=GREEN,
            ha="center", va="top", fontsize=8.5)
    ax.text((xa + xb) / 2, 8.6, "crossover\nwindow", color="#a0520c", ha="center", va="top", fontsize=8.5)
    ax.text((xb + 6) / 2, 5.9, "strong coupling:\nthis paper,\ncluster expansions", color=BLUE,
            ha="center", va="top", fontsize=8.5)
    shades = ["#111111", "#3b3b3b", "#666666", "#8f8f8f"]
    for x0, c, lab in zip([-6.5, -9.5, -12.5, -15.2], shades, ["$a_0$", "$a_0/8$", "$a_0/64$", "$a_0/512$"]):
        x = np.linspace(x0, 6, 600)
        ax.plot(x, g2(x), color=c, lw=1.4)
        ax.plot([x0], [g2(np.array([x0]))[0]], "o", color=c, ms=3.5)
        ax.text(x0, g2(np.array([x0]))[0] - 0.55, lab, color=c, ha="center", fontsize=8)
    ax.annotate("", xy=(xb, 0.35), xytext=(xa, 0.35),
                arrowprops=dict(arrowstyle="<->", color="#a0520c", lw=1.0))
    ax.text((xa + xb) / 2, 0.55, "same number of steps for every $a$", color="#a0520c", ha="center", fontsize=7.5)
    ax.set_xlim(-16, 6); ax.set_ylim(0, 8.8)
    ax.set_xlabel(r"block scale $\log_2(L a/\ell)$, with $\ell$ a fixed physical length")
    ax.set_ylabel("running coupling $g^2$")
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    save(fig, "fig_scales")


if __name__ == "__main__":
    for f in (fig_thresholds, fig_plaquette, fig_extremal, fig_spectrum, fig_weights,
              fig_sphere, fig_surface, fig_box, fig_massgap, fig_scales):
        f()
        print("done", f.__name__)
