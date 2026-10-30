"""Checks of the convex chart of the fibres of the cut-off average
(Theorem thm:fibrechart of the forest paper, Sec. ssec:smallfibres), for
G = SU(2) on the periodic lattice of side 4.

The backgrounds are those of small_field_fibre.py: block gauge transforms of
fields whose plaquettes have angle of order eps. For a set W of blocks, T is
the tangent space of the fibre F at U, and the chart is
Phi(eta) = U exp(eta + h(eta)), eta in T, h(eta) in the orthogonal complement
N = range DC^*, found by the iteration z -> z - G^{-1} C(eta + DC^* z) of the
proof. The program

  (f) solves for h at random points eta of the l^infinity ball K_s of T and
      checks that Phi(eta) lies on F, that h(eta) is orthogonal to T, that
      the iteration contracts, and records |h(eta)|_inf / |eta|_inf^2;
  (g) computes, by second differences, the Hessian of Q o Phi at such points
      in the direction of least Q_U on T and in random directions of T, and
      checks that it is at least 1/26;
  (h) computes the Hessian of psi = -sum_e log det Psi(ad X_e)
      + 1/2 log det(A_X A_X^*) - 1/2 log det(1 + Dh^* Dh), the logarithm of
      the Haar, coarea and area factors of the theorem, and records
      |psi''(zeta, zeta)| / |zeta|^2, which the theorem bounds by a constant
      c_6 that does not depend on lambda, W or the lattice;
  (j) checks the arithmetic of the constants of the proof exactly.

These are floating-point checks at fields far larger than those of the
theorem; they are not interval certificates.

Usage: python3 fibre_chart.py
"""

from fractions import Fraction as Fr

import numpy as np

import small_field_fibre as S

checks = []


def check(name, ok):
    checks.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name)


def setup(eps, W, seed):
    rng = np.random.default_rng(seed)
    U = S.make_field(eps, rng)
    kappa = max(S.dist(S.hol(U, p), S.I2) for p in S.PLAQ)
    Wset = set(W)
    var = [(x, mu) for x in S.SITES for mu in range(S.DIM)
           if not S.is_tree(x, mu) and S.blk(x) in Wset]
    vid = {l: i for i, l in enumerate(var)}
    cons = [(y, mu) for y in W for mu in range(S.DIM)]
    V0 = {k: S.vsharp(U, *k)[0] for k in cons}
    plist = [p for p in S.PLAQ if any(l in vid for l in S.plinks(p))]
    return dict(U=U, kappa=kappa, var=var, vid=vid, cons=cons, V0=V0,
                plist=plist, n=3 * len(var), rng=rng)


def moved(st, X):
    Uf = dict(st["U"])
    for l, i in st["vid"].items():
        Uf[l] = st["U"][l] @ S.expm(X[3 * i:3 * i + 3])
    return Uf


def Cval(st, X):
    Uf = moved(st, X)
    return np.concatenate([S.logm(st["V0"][k].conj().T @ S.vsharp(Uf, *k)[0])
                           for k in st["cons"]])


def Qval(st, X):
    Uf = moved(st, X)
    return sum(S.dist(S.hol(Uf, p), S.I2) ** 2 for p in st["plist"])


def jac(st, X, t=1e-5):
    """A_X = derivative of X -> C(U e^X), by central differences on the 24
    links of each constraint."""
    vid, cons = st["vid"], st["cons"]
    A = np.zeros((3 * len(cons), st["n"]))
    Uf = moved(st, X)
    for k, (y, mu) in enumerate(cons):
        Vk = st["V0"][(y, mu)]
        for l in S.vlinks(y, mu):
            if l not in vid:
                continue
            i = vid[l]
            base = Uf[l]
            for c in range(3):
                d = np.zeros(3)
                d[c] = t
                Uf[l] = st["U"][l] @ S.expm(X[3 * i:3 * i + 3] + d)
                Vp = S.vsharp(Uf, y, mu)[0]
                Uf[l] = st["U"][l] @ S.expm(X[3 * i:3 * i + 3] - d)
                Vm = S.vsharp(Uf, y, mu)[0]
                Uf[l] = base
                A[3 * k:3 * k + 3, 3 * i + c] = (S.logm(Vk.conj().T @ Vp) - S.logm(Vk.conj().T @ Vm)) / (2 * t)
    return A


def solve_h(st, eta, B, Ginv, iters=12):
    z = np.zeros(B.shape[0])
    steps = []
    for _ in range(iters):
        znew = z - Ginv @ Cval(st, eta + B.T @ z)
        steps.append(np.max(np.abs(znew - z)))
        z = znew
        if steps[-1] < 1e-14:
            break
    return B.T @ z, steps


def linf(X):
    return np.max(np.linalg.norm(X.reshape(-1, 3), axis=1))


def logdetPsi(c):
    """log det Psi(ad X) for X in su(2) with coordinates c (HS norm): the
    eigenvalues of ad X are 0 and +-i sqrt(2)|c|."""
    th = np.sqrt(2) * np.linalg.norm(c)
    if th < 1e-12:
        return 0.0
    # Psi(z) = (1 - e^{-z})/z at z = +-i th: |Psi|^2 = (2 - 2 cos th)/th^2
    return np.log((2 - 2 * np.cos(th)) / th ** 2)


def run(eps, W, seed, s, label):
    st = setup(eps, W, seed)
    rng = st["rng"]
    kappa = st["kappa"]
    tag = f"eps={eps:g}, |W|={len(W)}, kappa={kappa:.3e}, s={s:g}"
    print(f"-- {label}: {tag}")
    n = st["n"]
    B = jac(st, np.zeros(n))
    G = B @ B.T
    Ginv = np.linalg.inv(G)
    _, sv, Vt = np.linalg.svd(B)
    rank = int(np.sum(sv > 1e-8))
    Nb = Vt[rank:].T  # orthonormal basis of T

    # the covariant curvature form at U, for the least direction
    D = np.zeros((3 * len(st["plist"]), n))
    for r, p in enumerate(st["plist"]):
        A4 = S.dmaps(st["U"], p)
        for i, l in enumerate(S.plinks(p)):
            if l in st["vid"]:
                D[3 * r:3 * r + 3, 3 * st["vid"][l]:3 * st["vid"][l] + 3] += A4[i]
    QU = D.T @ D
    least = Nb @ np.linalg.eigh(Nb.T @ QU @ Nb)[1][:, 0]

    def sample_eta(frac):
        v = Nb @ rng.normal(size=Nb.shape[1])
        return frac * s * v / linf(v)

    def chart(eta):
        h, steps = solve_h(st, eta, B, Ginv)
        return eta + h, h, steps

    # (f) chart
    ok_on, ok_orth, ok_contr = True, True, True
    ratios = []
    pts = [sample_eta(f) for f in (0.5, 1.0, 1.0)]
    for eta in pts:
        X, h, steps = chart(eta)
        ok_on &= np.max(np.abs(Cval(st, X))) <= 1e-11
        ok_orth &= np.max(np.abs(Nb.T @ h)) <= 1e-12
        ok_contr &= all(steps[i + 1] <= 0.5 * steps[i] + 1e-15 for i in range(len(steps) - 1))
        ratios.append(linf(h) / linf(eta) ** 2)
    check(f"(f) Phi(eta) lies on F (|C| <= 1e-11) at three points of K_s ({tag})", ok_on)
    check(f"(f) h(eta) is orthogonal to T ({tag})", ok_orth)
    check(f"(f) the iteration for h contracts by at least 1/2 per step ({tag})", ok_contr)
    check(f"(f) |h|_inf/|eta|_inf^2 = {max(ratios):.3f} stays bounded (<= 50) ({tag})", max(ratios) <= 50)
    pk = max(S.dist(S.hol(moved(st, chart(eta)[0]), p), S.I2) for eta in pts[1:] for p in st["plist"])
    check(f"(f) plaquettes on Phi(K_s) at most kappa + 6 s ({pk:.4f} <= {kappa + 6 * s:.4f})",
          pk <= kappa + 6 * s)

    # (g) and (h): second differences along eta + t zeta
    def QPhi(eta):
        return Qval(st, chart(eta)[0])

    def psi(eta):
        X, _, _ = chart(eta)
        A = jac(st, X)
        haar = -sum(logdetPsi(X[3 * i:3 * i + 3]) for i in range(len(st["var"])))
        cor = 0.5 * np.linalg.slogdet(A @ A.T)[1]
        # Dh on T: -B^*(A B^*)^{-1} A
        Sm = -B.T @ np.linalg.solve(A @ B.T, A @ Nb)
        area = -0.5 * np.linalg.slogdet(np.eye(Nb.shape[1]) + Sm.T @ Sm)[1]
        return haar + cor + area, np.linalg.norm(Sm, 2)

    hq, hp, snorm = [], [], []
    base = pts[1]
    dirs = [least] + [Nb @ rng.normal(size=Nb.shape[1]) for _ in range(2)]
    for zeta in dirs:
        zeta = zeta / np.linalg.norm(zeta)
        t = 2e-3
        q0, qp, qm = QPhi(base), QPhi(base + t * zeta), QPhi(base - t * zeta)
        hq.append((qp - 2 * q0 + qm) / t ** 2)
        t2 = 3e-2
        p0, sn = psi(base)
        pp, _ = psi(base + t2 * zeta)
        pm, _ = psi(base - t2 * zeta)
        hp.append(abs(pp - 2 * p0 + pm) / t2 ** 2)
        snorm.append(sn)
    check(f"(g) Hessian of Q o Phi at a point of K_s, least and two random directions: "
          f"{', '.join(f'{v:.5f}' for v in hq)} (all >= 1/26) ({tag})", min(hq) >= 1 / 26)
    check(f"(h) |psi''(zeta,zeta)|/|zeta|^2 = {', '.join(f'{v:.3f}' for v in hp)} (<= 10) ({tag})",
          max(hp) <= 10)
    check(f"(h) |Dh| = {max(snorm):.4f} <= 1/6 on K_s ({tag})", max(snorm) <= 1 / 6)
    return max(hp)


def constants():
    c3 = 1  # every bound below is used with c_3 >= 1 and s <= 1/(6000 c_3)
    check("constants: 17*12*4 = 816 and 0.88*816 <= 720", 17 * 12 * 4 == 816 and Fr(88, 100) * 816 <= 720)
    check("constants: 17*48*0.88 <= 720 and 720/1440 = 1/2, 816/6000 <= 1",
          17 * 48 * Fr(88, 100) <= 720 and Fr(720, 1440) == Fr(1, 2) and Fr(816, 6000) <= 1)
    check("constants: 720 s^2 <= s/2 for s <= 1/1440", Fr(720, 1440) <= Fr(1, 2))
    check("constants: 2*0.44 = 0.88, 0.44*sqrt(24) <= 2.2, 0.88*17*2.2 <= 33, 17*2.2 <= 38, 34/40 < 1",
          2 * 0.44 == 0.88 and 0.44 * np.sqrt(24) <= 2.2 and 0.88 * 17 * 2.2 <= 33 and 17 * 2.2 <= 38)
    check("constants: sqrt(24)*1.5 <= 7.4, 7.4*sqrt(2) <= 10.5, 0.44*sqrt(2) <= 0.63",
          np.sqrt(24) * 1.5 <= 7.4 and 7.4 * np.sqrt(2) <= 10.5 and 0.44 * np.sqrt(2) <= 0.63)
    s = Fr(1, 6000)
    check("constants: 17*10.5*0.63 s <= 1/2 and 0.63*34*10.5 s <= 1/6 at s = 1/6000",
          17 * Fr(105, 10) * Fr(63, 100) * s <= Fr(1, 2) and Fr(63, 100) * 34 * Fr(105, 10) * s <= Fr(1, 6))
    check("constants: 24*0.44 <= 10.6, 2*7.4 = 14.8, 14.8*10.6 <= 157, 17*157 s <= 1/2",
          24 * Fr(44, 100) <= Fr(106, 10) and Fr(148, 10) * Fr(106, 10) <= 157 and 17 * 157 * s <= Fr(1, 2))
    check("constants: 1 + 1/36 <= 1.03 and 10.6*34*2*1.03 <= 744",
          1 + Fr(1, 36) <= Fr(103, 100) and Fr(106, 10) * 34 * 2 * Fr(103, 100) <= 744)
    check("constants: 3s*1.03 + (1 + 3s)*744 <= 750 at s = 1/6000",
          3 * s * Fr(103, 100) + (1 + 3 * s) * 744 <= 750)
    check("constants: 12*750 = 9000 and 6(5.4e4 + 9000 c3) + 1 <= 4e5 c3 for c3 >= 1",
          12 * 750 == 9000 and all(6 * (54000 + 9000 * c) + 1 <= 400000 * c for c in (1, 2, 10, 10 ** 6)))
    check("constants: (1 - 6s)/26 >= 1/26 - s/4", Fr(1, 26) * (1 - 6 * s) >= Fr(1, 26) - s / 4)
    check("constants: 1/sqrt(17) - 10.5 s >= 0.2, so |M^{-1}| <= 25",
          1 / np.sqrt(17) - 10.5 / 6000 >= 0.2)
    check("constants: 1/26 - 1/195 = 1/30 and 1/30 - 1/60 = 1/60",
          Fr(1, 26) - Fr(1, 195) == Fr(1, 30) and Fr(1, 30) - Fr(1, 60) == Fr(1, 60))
    check("constants: 0.63*34*(1 + 0.63*34) <= 481", Fr(63, 100) * 34 * (1 + Fr(63, 100) * 34) <= 481)


if __name__ == "__main__":
    constants()
    blocks = S.BLOCKS
    two = [blocks[0], blocks[1]]
    eight = [y for y in blocks if sum(y) % 2 == 0]
    c6 = []
    for eps, W, seed, s, lab in ((1e-3, two, 5, 1e-2, "two blocks"),
                                 (1e-3, eight, 6, 1e-2, "eight blocks"),
                                 (1e-2, eight, 7, 3e-2, "eight blocks"),
                                 (1e-3, blocks, 8, 1e-2, "sixteen blocks")):
        c6.append(run(eps, W, seed, s, lab))
    print(f"largest |psi''|/|zeta|^2 observed: {max(c6):.3f}")
    n_ok = sum(ok for _, ok in checks)
    print(f"{n_ok} of {len(checks)} checks passed")
