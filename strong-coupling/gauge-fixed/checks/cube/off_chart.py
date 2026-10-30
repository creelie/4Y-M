"""Checks of the mass off the chart (Theorem thm:offchart of the forest
paper, Sec. ssec:smallfibres), for G = SU(2) on the periodic lattice of side 4.

The backgrounds are those of small_field_fibre.py. For a set W of blocks the
program

  (k) checks the first step of the proof: if U and U' agree off A_W, have the
      same cut-off averages and all plaquettes at most kappa and kappa', then
      delta(U_e, U'_e) <= 92.3 max(kappa, kappa') on every link of A_W. The
      points U' are found by a search on the fibre: random moves of the links
      of A_W other than the path links, the path links then solved from the
      constraints by Newton's method, followed by descent of the action Q on
      the fibre; the inequality is checked at every iterate;
  (l) checks the covariance comparison |Cov_mu(F,G) - Cov_{mu(.|A)}(F,G)|
      <= 6 eps |F|_inf |G|_inf, eps = mu(A^c), in exact arithmetic on random
      finite laws;
  (m) checks the arithmetic of the constants of the theorem exactly, and
      the least beta_W with q_kappa(beta_W) < 1 at kappa = 1.4e-14.

The search of (k) is a floating-point computation at fields far larger than
those of the theorem; it is not a certificate.

Usage: python3 off_chart.py
"""

import math
import random
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
    cons = [(y, mu) for y in W for mu in range(S.DIM)]
    path = {k: (S.add(tuple(2 * c for c in k[0]), S.e(k[1])), k[1]) for k in cons}
    free = [l for l in var if l not in set(path.values())]
    V0 = {k: S.vsharp(U, *k)[0] for k in cons}
    varset = set(var)
    plist = [p for p in S.PLAQ if any(l in varset for l in S.plinks(p))]
    return dict(U=U, kappa=kappa, var=var, cons=cons, path=path, free=free,
                V0=V0, plist=plist, rng=rng)


def solve_paths(st, Uf):
    """Solve V^sharp_k(Uf) = V0_k by moving the path link d_k, which enters
    only V_k. Returns False if Newton's method fails or leaves the branch."""
    for k in st["cons"]:
        d = st["path"][k]
        for _ in range(30):
            try:
                c = S.logm(st["V0"][k].conj().T @ S.vsharp(Uf, *k)[0])
            except AssertionError:
                return False
            if np.max(np.abs(c)) < 1e-13:
                break
            J = np.zeros((3, 3))
            base = Uf[d]
            t = 1e-6
            for a in range(3):
                da = np.zeros(3)
                da[a] = t
                Uf[d] = base @ S.expm(da)
                cp = S.logm(st["V0"][k].conj().T @ S.vsharp(Uf, *k)[0])
                Uf[d] = base @ S.expm(-da)
                cm = S.logm(st["V0"][k].conj().T @ S.vsharp(Uf, *k)[0])
                J[:, a] = (cp - cm) / (2 * t)
            Uf[d] = base @ S.expm(-np.linalg.solve(J, c))
        else:
            return False
    return True


def Q(st, Uf):
    return sum(S.dist(S.hol(Uf, p), S.I2) ** 2 for p in st["plist"])


def gradQ(st, Uf):
    """Left-trivialised gradient of Q in the free links: dQ(xi) =
    -2 sum_p Re Tr(D_p(xi) U_p)."""
    g = {l: np.zeros(3) for l in st["free"]}
    for p in st["plist"]:
        Up = S.hol(Uf, p)
        w = np.array([-2 * np.real(np.trace(S.TAU[a] @ Up)) for a in range(3)])
        for i, l in enumerate(S.plinks(p)):
            if l in g:
                g[l] += S.dmaps(Uf, p)[i].T @ w
    return g


def measure(st, Uf):
    kp = max(S.dist(S.hol(Uf, p), S.I2) for p in S.PLAQ)
    dm = max(S.dist(st["U"][l], Uf[l]) for l in st["var"])
    return kp, dm


def run(eps, W, seed, sizes, steps, label):
    st = setup(eps, W, seed)
    rng = st["rng"]
    kappa = st["kappa"]
    tag = f"eps={eps:g}, |W|={len(W)}, kappa={kappa:.3e}"
    print(f"-- {label}: {tag}")
    worst, npts, final = 0.0, 0, []
    ok = True
    for a in sizes:
        Uf = dict(st["U"])
        for l in st["free"]:
            Uf[l] = st["U"][l] @ S.expm(rng.normal(size=3) * a / np.sqrt(3))
        if not solve_paths(st, Uf):
            print(f"   start of size {a}: outside the small-field branch, skipped")
            continue
        q = Q(st, Uf)
        for _ in range(steps):
            kp, dm = measure(st, Uf)
            npts += 1
            r = dm / max(kappa, kp)
            worst = max(worst, r)
            ok &= dm <= 92.3 * max(kappa, kp)
            g = gradQ(st, Uf)
            t = 0.05
            while t > 1e-6:
                Un = dict(Uf)
                for l in st["free"]:
                    Un[l] = Uf[l] @ S.expm(-t * g[l])
                if solve_paths(st, Un) and Q(st, Un) < q:
                    Uf, q = Un, Q(st, Un)
                    break
                t /= 2
            else:
                break
        kp, dm = measure(st, Uf)
        final.append((a, kp, dm))
        ok &= dm <= 92.3 * max(kappa, kp)
        npts += 1
        worst = max(worst, dm / max(kappa, kp))
    for a, kp, dm in final:
        print(f"   start of size {a}: final kappa' = {kp:.3e}, max delta(U_e, U'_e) = {dm:.3e}")
    check(f"(k) delta(U_e, U'_e) <= 92.3 max(kappa, kappa') at {npts} points of the fibre, "
          f"largest ratio {worst:.2f} ({tag})", ok and npts > 0)
    return worst


def covariance_check(trials=300, seed=11):
    rnd = random.Random(seed)
    ok = True
    for _ in range(trials):
        n = rnd.randint(3, 9)
        w = [Fr(rnd.randint(1, 50)) for _ in range(n)]
        Z = sum(w)
        p = [x / Z for x in w]
        F = [Fr(rnd.randint(-20, 20), 7) for _ in range(n)]
        G = [Fr(rnd.randint(-20, 20), 5) for _ in range(n)]
        A = [rnd.random() < 0.8 for _ in range(n)]
        if not any(A):
            continue
        pa = sum(pi for pi, inA in zip(p, A) if inA)
        e = 1 - pa
        E = lambda f, m: sum(mi * fi for mi, fi in zip(m, f))
        q = [pi / pa if inA else Fr(0) for pi, inA in zip(p, A)]
        cov = lambda m: E([f * g for f, g in zip(F, G)], m) - E(F, m) * E(G, m)
        lhs = abs(cov(p) - cov(q))
        ok &= lhs <= 6 * e * max(abs(x) for x in F) * max(abs(x) for x in G)
    check(f"(l) |Cov_mu - Cov_mu(.|A)| <= 6 mu(A^c) |F|_inf |G|_inf on {trials} random finite laws", ok)


def constants():
    check("(m) 16 pi <= 50.3, 3 + 3 = 6 <= 92.3 and 21 + 21 + 50.3 = 92.3",
          16 * math.pi <= 50.3 and 21 + 21 + Fr(503, 10) == Fr(923, 10))
    check("(m) (pi/3) 92.3 <= 96.7 and 1280 * 96.7 <= 1.24e5",
          math.pi / 3 * 92.3 <= 96.7 and 1280 * Fr(967, 10) <= 124000)
    check("(m) 92.3 kappa <= 1 and <= rho_G when kappa <= 1e-8 rho_G", 92.3e-8 <= 1)
    check("(m) blocks at distance at most two: 5^4 = 625, and 625 * 16 = 10^4", 5 ** 4 * 16 == 10 ** 4)
    check("(m) half the diameter of SU(2) is sqrt(2) pi / 2 <= 2.23 and 6 * 2.23^2 < 30",
          math.sqrt(2) * math.pi / 2 <= 2.23 and 6 * Fr(223, 100) ** 2 < 30)
    check("(m) kappa + 6 s <= 1e-8 with s = 1.24e5 kappa forces kappa <= 1.4e-14",
          Fr(1, 10 ** 8) / (1 + 6 * 124000) <= Fr(14, 10 ** 15))
    C0 = 0.5 + 2 / 3 * math.log(24 * math.pi ** 3)
    k = 1.4e-14
    b = 1.0
    while b * k * k <= 96 * (math.log(b) + C0):
        b *= 1.01
    check(f"(m) q_kappa(beta) < 1 at kappa = 1.4e-14 needs beta > {b:.2e} > 1e31", b > 1e31)


if __name__ == "__main__":
    constants()
    covariance_check()
    blocks = S.BLOCKS
    two = [blocks[0], blocks[1]]
    eight = [y for y in blocks if sum(y) % 2 == 0]
    worst = []
    for eps, W, seed, sizes, steps, lab in (
            (1e-3, two, 21, (0.02, 0.05, 0.1, 0.2), 25, "two blocks"),
            (1e-2, two, 22, (0.02, 0.05, 0.1, 0.2), 25, "two blocks"),
            (1e-3, eight, 23, (0.05, 0.1), 15, "eight blocks")):
        worst.append(run(eps, W, seed, sizes, steps, lab))
    print(f"largest ratio delta / max(kappa, kappa') observed: {max(worst):.2f}")
    n_ok = sum(ok for _, ok in checks)
    print(f"{n_ok} of {len(checks)} checks passed")
