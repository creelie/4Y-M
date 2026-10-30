"""Checks of the product chart and of the decay of the chart Hessian
(Theorem thm:chartdecay of the forest paper, Sec. ssec:smallfibres), for
G = SU(2) on the periodic lattice of side 12, with W a row of five blocks.

For every constraint k = (y, mu) the path link d_k = (2y + e_mu, mu) enters
V_{y,mu} with weight 1/8 and no other average. The section S writes an
element of T = ker B through its coordinates xi on the other links: the
coordinate at d_k solves the constraint k. The program

  (k) checks that d_k enters only the constraint k, that B S = 0, that S is
      injective with image T, records max_k |B_{k,d_k}^{-1}| and the largest
      l^infinity block row sum of S at the path links, which the proof
      bounds by 8.01 and 15.1 for kappa <= 1e-8;
  (l) checks that Pbar = P S, P = 1 - B^*(A B^*)^{-1} A, has Pbar^* Pbar >= 1
      at points of the chart;
  (m) records, against the l^infinity distance between link midpoints, the
      largest block norms of G^{-1} = (B B^*)^{-1}, and of Pbar at a point of
      the chart, and checks that they decay (at the base point Pbar = S has
      finite range, since there A = B and P S = S);
  (n) computes, by second differences of Q(U exp(X(xi))) with X(xi) =
      S xi + h(S xi), blocks of the Hessian in xi at that point for pairs
      of links at increasing distance, and checks that their norms decay;
  (o) checks the arithmetic of the constants of the proof exactly.

These are floating-point checks at fields far larger than those of the
theorem; they are not interval certificates.

Usage: python3 chart_decay.py
"""

import itertools
from fractions import Fraction as Fr

import numpy as np

import small_field_fibre as S

# a larger lattice, so that distances up to nine occur inside W
S.L = 12
S.SITES = list(itertools.product(range(S.L), repeat=S.DIM))
S.BLOCKS = list(itertools.product(range(S.L // 2), repeat=S.DIM))
S.PLAQ = [(x, mu, nu) for x in S.SITES for mu in range(S.DIM) for nu in range(mu + 1, S.DIM)]

import fibre_chart as FC  # noqa: E402

checks = FC.checks
check = FC.check


def mid(l):
    x, mu = l
    return np.array(x, float) + 0.5 * np.eye(S.DIM)[mu]


def dist(l1, l2):
    return float(np.max(np.abs(mid(l1) - mid(l2))))


def blocknorm(M, i, j):
    return np.linalg.norm(M[3 * i:3 * i + 3, 3 * j:3 * j + 3], 2)


def envelope(pairs):
    """pairs of (distance, norm) -> sorted list of (d, max norm at distance >= d)."""
    ds = sorted({round(d, 6) for d, _ in pairs})
    return [(d, max(v for dd, v in pairs if dd >= d - 1e-9)) for d in ds]


def run(eps, seed, s):
    W = [(j, 0, 0, 0) for j in range(5)]
    st = FC.setup(eps, W, seed)
    kappa = st["kappa"]
    tag = f"eps={eps:g}, |W|={len(W)}, kappa={kappa:.3e}"
    print(f"-- row of five blocks: {tag}")
    var, vid, cons, n = st["var"], st["vid"], st["cons"], st["n"]
    B = FC.jac(st, np.zeros(n))
    G = B @ B.T

    # (k) the path links and the section
    des = [(tuple((2 * y[i] + (i == mu)) % S.L for i in range(S.DIM)), mu) for (y, mu) in cons]
    owners = {d: [k for k, (y, mu) in enumerate(cons) if d in S.vlinks(y, mu)] for d in des}
    check(f"(k) each path link d_k is active, in W, and enters only V_k ({tag})",
          all(d in vid for d in des) and all(owners[d] == [k] for k, d in enumerate(des)))
    dset = set(vid[d] for d in des)
    free = [i for i in range(len(var)) if i not in dset]
    fidx = np.concatenate([np.arange(3 * i, 3 * i + 3) for i in free])
    Smat = np.zeros((n, len(fidx)))
    Smat[fidx, np.arange(len(fidx))] = 1.0
    binv = []
    for k, d in enumerate(des):
        r = slice(3 * k, 3 * k + 3)
        dc = slice(3 * vid[d], 3 * vid[d] + 3)
        Bkd = B[r, dc]
        Binv = np.linalg.inv(Bkd)
        binv.append(np.linalg.norm(Binv, 2))
        Bk = B[r, :].copy()
        Bk[:, dc] = 0
        Smat[dc, :] = -Binv @ (Bk @ Smat)
    check(f"(k) B S = 0 (|B S| = {np.max(np.abs(B @ Smat)):.1e}) ({tag})",
          np.max(np.abs(B @ Smat)) <= 1e-12)
    sv = np.linalg.svd(Smat, compute_uv=False)
    check(f"(k) S is injective (least singular value {sv[-1]:.3f} >= 1) and "
          f"dim T = {Smat.shape[1]} = n - 3|K| ({tag})",
          sv[-1] >= 1 - 1e-12 and Smat.shape[1] == n - 3 * len(cons))
    rows = [sum(np.linalg.norm(Smat[3 * vid[d]:3 * vid[d] + 3, 3 * j:3 * j + 3], 2)
                for j in range(len(free))) for d in des]
    check(f"(k) max_k |B_(k,d_k)^-1| = {max(binv):.3f}, path-link row sums of S "
          f"<= {max(rows):.3f} (<= 15.1) ({tag})", max(rows) <= 15.1)

    # (m) decay of G^{-1}
    Ginv = np.linalg.inv(G)
    gpairs = [(dist(des[a], des[b]), blocknorm(Ginv, a, b)) for a in range(len(cons))
              for b in range(len(cons))]
    env = envelope(gpairs)
    line = ", ".join(f"{d:g}:{v:.1e}" for d, v in env)
    far = [v for d, v in env if d >= 6]
    check(f"(m) largest block of G^-1 at distance >= d: {line} "
          f"(at d >= 6 at most a hundredth of d = 0) ({tag})", bool(far) and far[0] <= 0.01 * env[0][1])

    # (l) and (n): the chart in the coordinate xi
    Gi = np.linalg.inv(G)

    def X_of(xi):
        eta = Smat @ xi
        h, _ = FC.solve_h(st, eta, B, Gi)
        return eta + h

    rng = st["rng"]
    xi0 = rng.normal(size=len(fidx))
    xi0 *= 0.3 * s / FC.linf(Smat @ xi0)
    X0 = X_of(xi0)
    A = FC.jac(st, X0)
    Pb0 = (np.eye(n) - B.T @ np.linalg.solve(A @ B.T, A)) @ Smat
    lam = np.linalg.eigvalsh(Pb0.T @ Pb0)[0]
    check(f"(l) Pbar^* Pbar >= 1 at a point of the chart (least eigenvalue {lam:.4f}) ({tag})",
          lam >= 1 - 1e-10)
    ppairs = []
    for jcol, i in enumerate(free):
        for a in range(len(var)):
            ppairs.append((dist(var[a], var[i]), blocknorm(Pb0, a, jcol)))
    env = envelope(ppairs)
    line = ", ".join(f"{d:g}:{v:.1e}" for d, v in env if abs(d - round(d)) < 1e-9)
    near = [v for d, v in env if d >= 3]
    far = [v for d, v in env if d >= 8]
    check(f"(m) largest block of Pbar at that point at distance >= d: {line} "
          f"(beyond the range of S it decays: d >= 8 at most a fiftieth of d >= 3) ({tag})",
          bool(far) and far[0] <= 0.02 * near[0])

    def Q(xi):
        return FC.Qval(st, X_of(xi))

    # pairs of free links along the row, at increasing distance from the first
    first = min(free, key=lambda i: (var[i][0][0], var[i][1]))
    cand = {}
    for i in free:
        d = dist(var[first], var[i])
        if abs(d - round(d)) < 1e-9 and int(round(d)) not in cand and var[i][1] == var[first][1]:
            cand[int(round(d))] = i
    t = 2e-3
    hpairs = []
    for d in sorted(cand):
        i = cand[d]
        ci, cf = free.index(i), free.index(first)
        blk = np.zeros((3, 3))
        for a in range(3):
            for b in range(3):
                va = np.zeros(len(fidx))
                vb = np.zeros(len(fidx))
                va[3 * cf + a] = t
                vb[3 * ci + b] = t
                blk[a, b] = (Q(xi0 + va + vb) - Q(xi0 + va - vb) - Q(xi0 + vb - va) + Q(xi0 - va - vb)) / (4 * t * t)
        hpairs.append((d, np.linalg.norm(blk, 2)))
    line = ", ".join(f"{d}:{v:.1e}" for d, v in hpairs)
    big = [v for d, v in hpairs if d <= 1]
    far = [v for d, v in hpairs if d >= 4]
    check(f"(n) Hessian blocks of Q o X in xi at that point against distance: {line} "
          f"(at d >= 4 at most a hundredth of d <= 1) ({tag})",
          bool(far) and bool(big) and max(far) <= 0.01 * max(big))


def constants():
    k = Fr(1, 10 ** 8)
    check("constants: 1/(1/8 - 260 kappa) <= 8.01 for kappa <= 1e-8", 1 / (Fr(1, 8) - 260 * k) <= Fr(801, 100))
    check("constants: 7/8 + 16/16 = 15/8 and 8.01 (15/8 + sqrt(24) 260e-8) <= 15.1",
          Fr(7, 8) + 1 == Fr(15, 8) and 8.01 * (15 / 8 + np.sqrt(24) * 260e-8) <= 15.1)
    check("constants: 15.1/32 < 1/2 and 15.1 (7/4)/32 < 15/16, 1.75/32 < 1/16",
          Fr(151, 320) < Fr(1, 2) and Fr(151, 10) * Fr(7, 4) / 32 < Fr(15, 16) and Fr(175, 3200) < Fr(1, 16))
    check("constants: v'' = 1 on [5r/4, oo) gives v'(t)/t >= 1 - (5/4)/(3/2) = 1/6 for t >= 3r/2",
          1 - Fr(5, 4) / Fr(3, 2) == Fr(1, 6))
    check("constants: lambda/30 - c6 >= lambda/60 and 1/60 - 1/120 = 1/120, 2 * 60 = 120",
          Fr(1, 60) - Fr(1, 120) == Fr(1, 120) and 1 / (Fr(1, 120)) == 120)


if __name__ == "__main__":
    constants()
    for eps, seed, s in ((1e-3, 11, 1e-2), (1e-2, 12, 3e-2)):
        run(eps, seed, s)
    n_ok = sum(ok for _, ok in checks)
    print(f"{n_ok} of {len(checks)} checks passed")
