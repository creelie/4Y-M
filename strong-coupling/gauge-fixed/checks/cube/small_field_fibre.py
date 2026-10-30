"""Checks of the small-field estimates on the fibres of the cut-off average
(Theorem thm:smallfibre and Lemma lem:smallframes of the forest paper,
Sec. ssec:smallfibres), for G = SU(2) on the periodic lattice of side 4.

The background fields are U = h . exp(eps X) in the tree gauge: X is a
random field on the active links, the tree links are I, and h is a random
gauge transformation that is constant on every block, so the links between
blocks are far from I while every plaquette has angle of order eps. For
eps = 1e-3, 1e-2 and 3e-2, and for W = all sixteen blocks and W = eight of
them, the program

  (a) checks the local bounds delta(U_e, I) <= 3 kappa on the active links
      inside the blocks and delta(U_a, v) <= 21 kappa on the links leaving
      a block, kappa being the largest plaquette distance;
  (b) checks the frame estimates of the lemma: the covariant curvature of an
      inner plaquette differs from the flat one by at most 18 kappa times
      the sum of the link norms, that of a face plaquette, in the frame of
      the next block, by at most 48 kappa times it, and the right-trivialised
      derivative of the cut-off average differs from the block average of
      the framed field by at most 260 kappa N_{y,mu} (operator norms, the
      derivative computed by central differences);
  (c) computes the tangent space of the fibre and the least value of
      Q_U(xi)/|xi|^2 on it, which the theorem bounds below by
      1/(52(1 + 1.4e6 kappa)) for kappa <= 1e-8 (the block gap of the
      quadratic approximation is 0.1075);
  (d) assembles the Hessian of Q along the geodesics U exp(t xi) from the
      formula of the proof, compares it with second differences, and checks
      Hess Q >= 2 Q_U - 84 kappa |xi|^2 for every xi (least eigenvalue);
  (e) checks the gradient bound |grad_e Q| <= 12 kappa, the Gram bound
      DC DC^* >= 1/16 - 700 kappa, and evaluates the Hessian of Q on the
      fibre, Hess Q + dQ(eta) with eta the normal acceleration, at random
      tangent vectors and at the least direction of (c), checking that the
      fibre curve U exp(t xi + t^2 eta / 2) stays on the fibre to third
      order.

These are floating-point checks of the inequalities of the proof at fields
far larger than the kappa <= 1e-8 of the theorem; they are not interval
certificates.

Usage: python3 small_field_fibre.py
"""

import itertools

import numpy as np

checks = []


def check(name, ok):
    checks.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name)


L = 4
DIM = 4
SIG = [np.array([[0, 1], [1, 0]], complex),
       np.array([[0, -1j], [1j, 0]], complex),
       np.array([[1, 0], [0, -1]], complex)]
TAU = [1j * s / np.sqrt(2) for s in SIG]  # orthonormal in the HS norm
I2 = np.eye(2, dtype=complex)


def alg(c):
    return sum(c[k] * TAU[k] for k in range(3))


def coords(X):
    return np.array([np.real(np.trace(TAU[k].conj().T @ X)) for k in range(3)])


def expm(c):
    a = np.asarray(c) / np.sqrt(2)
    th = np.linalg.norm(a)
    if th < 1e-300:
        return I2.copy()
    n = a / th
    return np.cos(th) * I2 + 1j * np.sin(th) * sum(n[k] * SIG[k] for k in range(3))


def logm(g):
    a0 = np.real(np.trace(g)) / 2
    b = np.array([np.real(np.trace(-1j * SIG[k] @ g)) / 2 for k in range(3)])
    nb = np.linalg.norm(b)
    th = np.arctan2(nb, a0)
    if nb < 1e-300:
        return np.zeros(3)
    return np.sqrt(2) * th * b / nb  # coordinates in TAU


def dist(a, b):
    return np.linalg.norm(a - b)


def ad(g):
    R = np.zeros((3, 3))
    for l in range(3):
        R[:, l] = coords(g @ TAU[l] @ g.conj().T)
    return R


def e(mu):
    v = [0] * DIM
    v[mu] = 1
    return tuple(v)


def add(x, y):
    return tuple((x[i] + y[i]) % L for i in range(DIM))


def sub(x, y):
    return tuple((x[i] - y[i]) % L for i in range(DIM))


SITES = list(itertools.product(range(L), repeat=DIM))
BLOCKS = list(itertools.product(range(L // 2), repeat=DIM))


def off(x):
    return tuple(c % 2 for c in x)


def blk(x):
    return tuple(c // 2 for c in x)


def is_tree(x, mu):
    o = off(x)
    return o[mu] == 0 and all(o[n] == 0 for n in range(mu + 1, DIM))


def make_field(eps, rng):
    U = {}
    h = {y: expm(rng.normal(size=3) * 3) for y in BLOCKS}
    for x in SITES:
        for mu in range(DIM):
            if is_tree(x, mu):
                g = I2.copy()
            else:
                g = expm(rng.normal(size=3) * eps)
            U[(x, mu)] = h[blk(x)] @ g @ h[blk(add(x, e(mu)))].conj().T
    return U


PLAQ = [(x, mu, nu) for x in SITES for mu in range(DIM) for nu in range(mu + 1, DIM)]


def plinks(p):
    x, mu, nu = p
    return [(x, mu), (add(x, e(mu)), nu), (add(x, e(nu)), mu), (x, nu)]


def hol(U, p):
    l1, l2, l3, l4 = (U[l] for l in plinks(p))
    return l1 @ l2 @ l3.conj().T @ l4.conj().T


def dmaps(U, p):
    """3x3 blocks A_i with D_p = sum_i A_i xi_i (coordinates)."""
    l1, l2, l3, l4 = (U[l] for l in plinks(p))
    g1, g2, g4 = l1, l1 @ l2, l1 @ l2 @ l3.conj().T
    return [ad(g1), ad(g2), -ad(g2), -ad(g4)]


def vsharp(U, y, mu):
    root = tuple(2 * c for c in y)
    P = U[(root, mu)] @ U[(add(root, e(mu)), mu)]
    Ms = []
    mx = 0.0
    for o in itertools.product((0, 1), repeat=DIM):
        x = add(root, o)
        Px = U[(x, mu)] @ U[(add(x, e(mu)), mu)]
        mx = max(mx, dist(Px, P))
        Ms.append(logm(P.conj().T @ Px))
    assert mx <= 0.5, "outside the small-field branch"
    return P @ expm(np.mean(Ms, axis=0)), P


def vlinks(y, mu):
    """The 24 active-or-tree links entering V_{y,mu} with their weights."""
    root = tuple(2 * c for c in y)
    out = {}
    for o in itertools.product((0, 1), repeat=DIM):
        x = add(root, o)
        for z in (x, add(x, e(mu))):
            out[(z, mu)] = out.get((z, mu), 0) + 1 / 16
    return out


def run(eps, Wfrac, seed):
    rng = np.random.default_rng(seed)
    U = make_field(eps, rng)
    kappa = max(dist(hol(U, p), I2) for p in PLAQ)
    W = BLOCKS if Wfrac == 1 else [y for y in BLOCKS if sum(y) % 2 == 0]
    Wset = set(W)
    tag = f"eps={eps:g}, |W|={len(W)}, kappa={kappa:.3e}"
    print(f"-- {tag}")

    # (a) local bounds
    inner = [(x, mu) for x in SITES for mu in range(DIM)
             if not is_tree(x, mu) and off(x)[mu] == 0]
    ok_in = max(dist(U[l], I2) for l in inner) <= 3 * kappa
    ok_out = True
    for y in BLOCKS:
        root = tuple(2 * c for c in y)
        for mu in range(DIM):
            v = U[(add(root, e(mu)), mu)]
            for a in itertools.product((0, 1), repeat=DIM):
                if a[mu] == 1:
                    ok_out &= dist(U[(add(root, a), mu)], v) <= 21 * kappa
    check(f"(a) inner links <= 3 kappa, leaving links within 21 kappa of v ({tag})",
          ok_in and ok_out)

    # variables
    var = [(x, mu) for x in SITES for mu in range(DIM)
           if not is_tree(x, mu) and blk(x) in Wset]
    vid = {l: i for i, l in enumerate(var)}
    n = 3 * len(var)

    # (b) frame estimates for plaquettes; assemble D
    rows = []
    ok_inner, ok_face = True, True
    worst_in, worst_face = 0.0, 0.0
    plist = [p for p in PLAQ if any(l in vid for l in plinks(p))]
    D = np.zeros((3 * len(plist), n))
    for r, p in enumerate(plist):
        A = dmaps(U, p)
        ls = plinks(p)
        for i, l in enumerate(ls):
            if l in vid:
                D[3 * r:3 * r + 3, 3 * vid[l]:3 * vid[l] + 3] += A[i]
        x, mu, nu = p
        bl = [blk(l[0]) for l in ls] + [blk(add(l[0], e(l[1]))) for l in ls]
        if len(set(bl)) == 1:
            sgn = [1, 1, -1, -1]
            err = max(np.linalg.norm(A[i] - sgn[i] * np.eye(3), 2) for i in range(4))
            worst_in = max(worst_in, err / kappa)
            ok_inner &= err <= 18 * kappa + 1e-12
        else:
            # face plaquette: the two links in the direction that leaves
            for lead, tran in ((mu, nu), (nu, mu)):
                if off(x)[lead] == 1 and blk(add(x, e(lead))) != blk(x) \
                        and blk(add(x, e(tran))) == blk(x):
                    y = blk(x)
                    root = tuple(2 * c for c in y)
                    v = U[(add(root, e(lead)), lead)]
                    Rv = ad(v.conj().T)
                    Bk = [Rv @ Ai for Ai in A]
                    if lead == mu:
                        want = [np.eye(3), np.eye(3), -np.eye(3), -Rv]
                    else:  # orientation reversed: compare with the inverse
                        Up = hol(U, p)
                        Bk = [-ad(v.conj().T @ Up.conj().T) @ Ai for Ai in A]
                        # links: (x,nu) in B(y), (x+nu,mu) leaving, (x+mu,nu) next block, (x,mu) leaving
                        want = [-Rv, -np.eye(3), np.eye(3), np.eye(3)]
                    err = max(np.linalg.norm(Bk[i] - want[i], 2) for i in range(4))
                    worst_face = max(worst_face, err / kappa)
                    ok_face &= err <= 48 * kappa + 1e-12
    check(f"(b) inner plaquettes within 18 kappa (worst {worst_in:.2f} kappa) ({tag})", ok_inner)
    check(f"(b) face plaquettes within 48 kappa in the next frame (worst {worst_face:.2f} kappa) ({tag})", ok_face)

    # constraint derivatives
    cons = [(y, mu) for y in W for mu in range(DIM)]
    DC = np.zeros((3 * len(cons), n))
    t = 1e-6
    worst_c = 0.0
    hfun = []
    for k, (y, mu) in enumerate(cons):
        V0, v = vsharp(U, y, mu)
        wts = vlinks(y, mu)
        Bt = np.zeros((3, n))
        Rv = ad(v.conj().T)
        for l, wgt in wts.items():
            if l not in vid:
                continue
            (z, _) = l
            frame = Rv if blk(z) == y and off(z)[mu] == 0 else np.eye(3)
            Bt[:, 3 * vid[l]:3 * vid[l] + 3] += wgt * frame
            for c in range(3):
                old = U[l]
                d = np.zeros(3)
                d[c] = t
                U[l] = old @ expm(d)
                Vp, _ = vsharp(U, y, mu)
                U[l] = old @ expm(-d)
                Vm, _ = vsharp(U, y, mu)
                U[l] = old
                DC[3 * k:3 * k + 3, 3 * vid[l] + c] = (logm(V0.conj().T @ Vp) - logm(V0.conj().T @ Vm)) / (2 * t)
        cols = [3 * vid[l] + c for l in wts if l in vid for c in range(3)]
        R = (DC[3 * k:3 * k + 3] - Bt)[:, cols]
        err = np.linalg.norm(R, 2)
        worst_c = max(worst_c, err / kappa)
    check(f"(b) linearised constraint within 260 kappa of the block average (worst {worst_c:.2f} kappa) ({tag})",
          worst_c <= 260)

    # (c) tangent space and covariant gap
    _, s, Vt = np.linalg.svd(DC)
    rank = int(np.sum(s > 1e-8))
    Nb = Vt[rank:].T
    QU = D.T @ D
    lam = np.linalg.eigvalsh(Nb.T @ QU @ Nb)[0]
    check(f"(c) DC has full rank {rank} = {3 * len(cons)}", rank == 3 * len(cons))
    check(f"(c) least Q_U/|xi|^2 on the tangent space = {lam:.6f} >= 1/52", lam >= 1 / 52)

    # (d) Hessian of Q along U exp(t xi)
    H = np.zeros((n, n))
    grad = np.zeros(n)
    ok_g = True
    gl = np.zeros((len(var), 3))
    for r, p in enumerate(plist):
        ls = plinks(p)
        Up = hol(U, p)
        A = dmaps(U, p)
        idx = [vid.get(l) for l in ls]
        w = np.array([-2 * np.real(np.trace(TAU[k] @ (Up - I2))) for k in range(3)])
        for i in range(4):
            if idx[i] is not None:
                g = A[i].T @ w
                grad[3 * idx[i]:3 * idx[i] + 3] += g
                gl[idx[i]] += g

        def q(c12):
            Bs = [alg(A[i] @ c12[3 * i:3 * i + 3]) for i in range(4)]
            Dm = sum(Bs)
            C = sum(Bs[i] @ Bs[j] - Bs[j] @ Bs[i] for i in range(4) for j in range(i + 1, 4))
            return -2 * np.real(np.trace((Dm @ Dm + C) @ Up))
        act = [3 * i + c for i in range(4) if idx[i] is not None for c in range(3)]
        glob = [3 * idx[i] + c for i in range(4) if idx[i] is not None for c in range(3)]
        E = np.eye(12)
        qd = {a: q(E[a]) for a in act}
        for ia, a in enumerate(act):
            for ib, b in enumerate(act):
                if ib < ia:
                    continue
                val = qd[a] if a == b else (q(E[a] + E[b]) - qd[a] - qd[b]) / 2
                H[glob[ia], glob[ib]] += val
                if ia != ib:
                    H[glob[ib], glob[ia]] += val
    ok_g = np.max(np.linalg.norm(gl, axis=1)) <= 12 * kappa
    check(f"(e) gradient bound |grad_e Q| <= 12 kappa ({tag})", ok_g)

    def Qval(Uf):
        return sum(dist(hol(Uf, p), I2) ** 2 for p in plist)

    def moved(X):
        Uf = dict(U)
        for l, i in vid.items():
            Uf[l] = U[l] @ expm(X[3 * i:3 * i + 3])
        return Uf
    xi = rng.normal(size=n)
    xi /= np.linalg.norm(xi)
    h = 1e-3
    fd = (Qval(moved(h * xi)) - 2 * Qval(U) + Qval(moved(-h * xi))) / h ** 2
    fd1 = (Qval(moved(h * xi)) - Qval(moved(-h * xi))) / (2 * h)
    check(f"(d) Hessian formula agrees with second differences ({xi @ H @ xi:.6f} vs {fd:.6f})",
          abs(xi @ H @ xi - fd) <= 1e-4)
    check(f"(e) gradient formula agrees with first differences ({grad @ xi:.3e} vs {fd1:.3e})",
          abs(grad @ xi - fd1) <= 1e-6 + 1e-4 * abs(fd1))
    mineig = np.linalg.eigvalsh(H - 2 * QU + 84 * kappa * np.eye(n))[0]
    check(f"(d) Hess Q >= 2 Q_U - 84 kappa |xi|^2 (least eigenvalue {mineig:.3e}) ({tag})", mineig >= -1e-9)
    hmin = np.linalg.eigvalsh(Nb.T @ H @ Nb)[0]
    check(f"(d) least ambient Hessian on the tangent space {hmin:.5f} >= 1/26 - 5.4e4 kappa",
          hmin >= 1 / 26 - 5.4e4 * kappa)

    # (e) Gram matrix and fibre Hessian
    G = DC @ DC.T
    gmin = np.linalg.eigvalsh(G)[0]
    check(f"(e) Gram matrix >= 1/16 - 700 kappa (least eigenvalue {gmin:.5f})", gmin >= 1 / 16 - 700 * kappa)

    def Cvals(X):
        Uf = moved(X)
        out = []
        for (y, mu) in cons:
            V0, _ = vsharp(U, y, mu)
            V1, _ = vsharp(Uf, y, mu)
            out.append(logm(V0.conj().T @ V1))
        return np.concatenate(out)
    w_least = Nb @ np.linalg.eigh(Nb.T @ QU @ Nb)[1][:, 0]
    samples = [w_least] + [Nb @ rng.normal(size=Nb.shape[1]) for _ in range(3)]
    ratios = []
    ok_curve = True
    for s_ in samples:
        s_ = s_ / np.linalg.norm(s_)
        hh = 1e-3
        hvec = (Cvals(hh * s_) + Cvals(-hh * s_)) / hh ** 2
        eta = -DC.T @ np.linalg.solve(G, hvec)
        hf = s_ @ H @ s_ + grad @ eta
        ratios.append(hf)
        tt = 1e-2
        c3 = np.linalg.norm(Cvals(tt * s_ + tt ** 2 * eta / 2))
        ok_curve &= c3 <= 50 * tt ** 3
    check(f"(e) fibre curves stay on the fibre to third order ({tag})", ok_curve)
    check(f"(e) fibre Hessian of Q at the least and three random tangent directions: "
          f"{', '.join(f'{r:.5f}' for r in ratios)} (all >= 1/26 - 5.4e4 kappa - 1e3 kappa)",
          min(ratios) >= 1 / 26 - 5.5e4 * kappa)
    return kappa, lam


def constants():
    from fractions import Fraction as F
    k = F(1, 10 ** 8)
    check("constants: 3888*8 = 31104", 3888 * 8 == 31104)
    check("constants: 8*260^2 = 540800, 2*540800 + 3*4608 = 1095424 <= 1.1e6",
          8 * 260 ** 2 == 540800 and 2 * 540800 + 3 * 4608 == 1095424 <= 1100000)
    check("constants: (13/2)*31104 + 1.1e6 <= 1.31e6", F(13, 2) * 31104 + 1100000 <= 1310000)
    check("constants: 52(1+k)^2/(1 - 1.31e6 k (1+k)^2) <= 52(1 + 1.4e6 k) at k = 1e-8",
          (1 + k) ** 2 / (1 - F(131, 100) * 10 ** 6 * k * (1 + k) ** 2) <= 1 + F(14, 10) * 10 ** 6 * k)
    check("constants: 2/(52(1+1.4e6 k)) - 84 k >= 1/26 - 5.4e4 k for k <= 1e-8",
          F(2, 52) / (1 + F(14, 10) * 10 ** 6 * k) - 84 * k >= F(1, 26) - 54000 * k)
    check("constants: 202 + 103*0.434 + 42/4 <= 260", 202 + 103 * 0.434 + 42 / 4 <= 260)
    check("constants: (1/16 + 454 k)/(1/8 - 227 k) <= 0.51 and 1/((1/8 - 227k)(0.49)) <= 17",
          (F(1, 16) + 454 * k) / (F(1, 8) - 227 * k) <= F(51, 100)
          and 1 / ((F(1, 8) - 227 * k) * F(49, 100)) <= 17)
    check("constants: 16 pi <= 50.3, e^(50.3k) - 1 <= 51 k, 51(1 + 50.3k) <= 51.6",
          16 * np.pi <= 50.3 and np.expm1(50.3e-8) <= 51e-8 and 51 * (1 + 50.3e-8) <= 51.6)
    check("constants: 51 + 51.6 + 99 <= 202, 51.6 + 51 <= 103, sqrt(48)/16 <= 0.434",
          51 + 51.6 + 99 <= 202 and 51.6 + 51 <= 103 and np.sqrt(48) / 16 <= 0.434)
    check("constants: sqrt(3/16) <= 0.434, 2*0.434*260 + 260^2 k <= 227",
          np.sqrt(3 / 16) <= 0.434 and 2 * 0.434 * 260 + 260 ** 2 * 1e-8 <= 227)
    check("constants: 1/8 - 227k - 2(1/32 + 227k) >= 1/16 - 700k",
          F(1, 8) - 227 * k - 2 * (F(1, 32) + 227 * k) >= F(1, 16) - 700 * k)
    check("constants: (0.434 + 260k) 12 sqrt(24) <= 26, 17*26 = 442, 2*442 = 884",
          (0.434 + 260e-8) * 12 * np.sqrt(24) <= 26 and 17 * 26 == 442)


if __name__ == "__main__":
    constants()
    for eps, wf, seed in ((1e-3, 1, 1), (1e-3, 0, 2), (1e-2, 1, 3), (3e-2, 0, 4)):
        run(eps, wf, seed)
    n_ok = sum(ok for _, ok in checks)
    print(f"{n_ok} of {len(checks)} checks passed")
