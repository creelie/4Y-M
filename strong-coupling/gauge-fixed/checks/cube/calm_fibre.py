"""Checks of the density bound for one cut-off average (lemma averagedensity)
and of the Peierls bound on a calm fibre (proposition calmfibre), for SU(2).

Paths of one coarse link: P = l_0 = g, P_x = u_a l_a, P_x' = l_a w_a.  In the
region R_1 (maximum at a site x_a) the coordinates are h_a = g^{-1} u_a l_a,
and v_x' = g^{-1} u_a^{-1} g h_a w_a.  The script
(A) checks the constants: M = 1 + 2(1 - 7 Lambda/8)^{-3} < 2.4e4,
    1 + 7 Lambda/8 < 1.957, (7 Lambda/8) 46 < 44, (M^4 1.957^12)^{1/45} < 2.93,
    14 188^3 < 9.31e7, and h(B_rho) = (2 theta - sin 2 theta)/(2 pi) >= rho^3/14
    for rho <= 0.1, theta = rho/sqrt 2, together with the Monte Carlo value of
    the Haar measure of a small ball;
(B) on random points of R_1 and R_2: the determinant of the derivative of
    g -> V(g, h) is at least (1 - 7 Lambda/8)^3, and the fixed point iteration
    g -> c exp(-chi N(g)) started at many points of S_h returns to the same g;
(C) on calm backgrounds: the iteration l -> c exp(-N(l)) with perturbed inner
    and private links stays in the ball of radius 44 rho about the path link,
    converges, and |det dV/dl_0| <= 1.957^3.
Double precision; the proofs in the paper are by hand.
Usage: python3 calm_fibre.py
"""
import numpy as np
from scipy.linalg import expm, logm

checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)


S = [np.array(m, complex) for m in ([[0, 1], [-1, 0]], [[0, 1j], [1j, 0]], [[1j, 0], [0, -1j]])]
B = [s / np.sqrt(2) for s in S]          # orthonormal for Re tr(X^* Y)


def alg(c):
    return sum(ci * b for ci, b in zip(c, B))


def coords(X):
    return np.array([np.real(np.trace(b.conj().T @ X)) for b in B])


def haar(rng):
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    return np.array([[q[0] + 1j * q[1], q[2] + 1j * q[3]], [-q[2] + 1j * q[3], q[0] - 1j * q[1]]])


def lg(g):
    X = logm(g)
    return (X - X.conj().T) / 2


def dist(a, b):
    return np.linalg.norm(lg(a.conj().T @ b))


def delta(v):
    return np.linalg.norm(v - np.eye(2))


def phi(t):
    return min(1.0, max(0.0, 2 - 2 * t))


ALPHA = 2 * np.arcsin(1 / (2 * np.sqrt(2)))
LAM = ALPHA / np.sin(ALPHA)
K7 = 7 * LAM / 8
MARGIN = 1 - K7
M = 1 + 2 * MARGIN ** -3
check(f"M = 1 + 2(1 - 7 Lambda/8)^-3 = {M:.1f} < 2.4e4", M < 2.4e4)
check(f"1 + 7 Lambda/8 = {1 + K7:.6f} < 1.957", 1 + K7 < 1.957)
check(f"(7 Lambda/8) 46 = {46 * K7:.3f} < 44", 46 * K7 < 44)
check(f"(M^4 1.957^12)^(1/45) = {(M ** 4 * 1.957 ** 12) ** (1 / 45):.4f} < 2.93",
      (M ** 4 * 1.957 ** 12) ** (1 / 45) < 2.93)
check(f"14 188^3 = {14 * 188 ** 3} < 9.31e7", 14 * 188 ** 3 < 9.31e7)
rhos = np.linspace(1e-6, 0.1, 100001)
th = rhos / np.sqrt(2)
check("h(B_rho) >= rho^3/14 for rho <= 0.1",
      np.all((2 * th - np.sin(2 * th)) / (2 * np.pi) >= rhos ** 3 / 14))
rng = np.random.default_rng(7)
q = rng.normal(size=(2_000_000, 4)); q /= np.linalg.norm(q, axis=1)[:, None]
r0 = 0.6
t0 = r0 / np.sqrt(2)
# the ball of radius r0 about I: angle arccos(q_0) <= t0 (q_0 > 0 side only)
frac = np.mean(np.arccos(np.clip(q[:, 0], -1, 1)) <= t0)
check(f"Haar measure of the ball of radius 0.6: Monte Carlo {frac:.5f}, "
      f"formula {(2 * t0 - np.sin(2 * t0)) / (2 * np.pi):.5f}",
      abs(frac - (2 * t0 - np.sin(2 * t0)) / (2 * np.pi)) < 5e-4)


def near(rng, s):
    return expm(alg(rng.normal(size=3) * s / np.sqrt(3)))


def V_region(g, h, u, w, kind):
    """cut-off average in the coordinates of R_1 (kind 1) or R_2 (kind 2)."""
    vs = []
    for a in range(7):
        if kind == 1:
            vs += [h[a], g.conj().T @ u[a].conj().T @ g @ h[a] @ w[a]]
        else:
            vs += [g.conj().T @ u[a] @ g @ h[a] @ w[a].conj().T, h[a]]
    m = max(delta(v) for v in vs)
    mh = max(delta(v) for v in h)
    N = sum(lg(v) for v in vs) / 16
    return g @ expm(phi(mh) * N), m, mh, N


def in_S(g, h, u, w, kind):
    _, m, mh, _ = V_region(g, h, u, w, kind)
    return mh < 1 and m <= mh + 1e-12


least_det, ok_unique, tested, solved = 9.0, True, 0, 0
for kind in (1, 2):
    for t in range(60):
        g = haar(rng)
        sc = rng.uniform(0.1, 0.6)
        h = [near(rng, sc) for _ in range(7)]
        mh = max(delta(v) for v in h)
        u = [haar(rng) for _ in range(7)]
        # choose the fixed links w_a so that the other seven sites are inside the ball of radius m_h
        w = []
        for a in range(7):
            s_a = near(rng, rng.uniform(0, 0.9) * mh / 1.5)
            if kind == 1:
                w.append(np.linalg.inv(g.conj().T @ u[a].conj().T @ g @ h[a]) @ s_a)
            else:
                w.append((np.linalg.inv(g.conj().T @ u[a] @ g @ h[a]) @ s_a).conj().T)
        if not in_S(g, h, u, w, kind):
            continue
        tested += 1
        V0 = V_region(g, h, u, w, kind)[0]
        cols = []
        for k in range(3):
            e = 1e-6
            Vp = V_region(expm(e * B[k]) @ g, h, u, w, kind)[0]
            Vm = V_region(expm(-e * B[k]) @ g, h, u, w, kind)[0]
            cols.append(coords(lg(Vp @ V0.conj().T) - lg(Vm @ V0.conj().T)) / (2 * e))
        least_det = min(least_det, abs(np.linalg.det(np.array(cols).T)))
        c = V0
        for start in range(4):
            x = g @ near(rng, rng.uniform(0.05, 0.5))
            if not in_S(x, h, u, w, kind):
                continue
            for it in range(300):
                _, _, mh2, N = V_region(x, h, u, w, kind)
                x = c @ expm(-phi(mh2) * N)
            if in_S(x, h, u, w, kind) and np.linalg.norm(V_region(x, h, u, w, kind)[0] - c) < 1e-9:
                solved += 1
                ok_unique &= dist(x, g) < 1e-7
print(f"regions R_1, R_2: {tested} sampled points, least |det| {least_det:.4f}, {solved} solutions found from other starts")
check(f"|det dV/dg| >= (1 - 7 Lambda/8)^3 = {MARGIN ** 3:.2e} on R_1 and R_2",
      tested > 50 and least_det >= MARGIN ** 3)
check("every solution of V(g, h) = c in S_h found by the iteration is the sampled g", ok_unique and solved > 20)

worst_ratio, max_det, conv = 0.0, 0.0, True
for t in range(200):
    kap = 10 ** rng.uniform(-6, -3)
    l0 = haar(rng)
    u = [near(rng, kap) for _ in range(7)]
    w = [near(rng, kap) for _ in range(7)]
    l = [l0 @ near(rng, 10 * kap) for _ in range(7)]

    def Nk(x0, u, l, w):
        vs = []
        for a in range(7):
            vs += [x0.conj().T @ u[a] @ l[a], x0.conj().T @ l[a] @ w[a]]
        return sum(lg(v) for v in vs) / 16

    c = l0 @ expm(Nk(l0, u, l, w))
    rho = 10 ** rng.uniform(-6, -3)
    def pert(g):
        z = rng.normal(size=3)
        return g @ expm(alg(z / np.linalg.norm(z)) * rho * rng.uniform())
    u2, l2, w2 = [pert(x) for x in u], [pert(x) for x in l], [pert(x) for x in w]
    x = l0
    for it in range(200):
        x = c @ expm(-Nk(x, u2, l2, w2))
        worst_ratio = max(worst_ratio, dist(x, l0) / rho)
    V = x @ expm(Nk(x, u2, l2, w2))
    conv &= np.linalg.norm(V - c) < 1e-10
    cols = []
    for k in range(3):
        e = 1e-7
        Vp = expm(e * B[k]) @ x @ expm(Nk(expm(e * B[k]) @ x, u2, l2, w2))
        Vm = expm(-e * B[k]) @ x @ expm(Nk(expm(-e * B[k]) @ x, u2, l2, w2))
        cols.append(coords(lg(Vp @ V.conj().T) - lg(Vm @ V.conj().T)) / (2 * e))
    max_det = max(max_det, abs(np.linalg.det(np.array(cols).T)))
print(f"calm backgrounds: largest d(l, l_0)/rho along the iteration {worst_ratio:.3f}, "
      f"largest |det| {max_det:.4f}")
check("the iteration stays within 44 rho of the path link and solves the constraint",
      conv and worst_ratio <= 44)
check(f"|det dV/dl_0| <= 1.957^3 = {1.957 ** 3:.3f} on calm backgrounds", max_det <= 1.957 ** 3)
print(f"{sum(checks)} of {len(checks)} checks passed")
