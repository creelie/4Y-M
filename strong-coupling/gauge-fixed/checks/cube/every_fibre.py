"""Checks of the aligned field (lemma alignedfield) and of the Peierls bound on
every fibre (proposition everyfibre), for SU(2).

A run of n blocks y_1, ..., y_n in direction mu, with y_n + e_mu outside W:
each pair has its path link l_0 = g, the links u_a = (x_a, mu) inside its
block (tree links, equal to I, when a_nu = 0 for every nu > mu), the private
links l_a, and w_a = u_a of the next block, fixed and arbitrary after y_n.
The script
(A) builds the aligned field for random values of the averages and random
    fixed links, in the four directions, and checks that every v_x equals I
    and that every average takes its prescribed value;
(B) on regular backgrounds with Haar random links inside the blocks and the
    sites of the averages within 1/4 of I, perturbs the links other than the
    path links by at most rho and checks that the iteration l -> c exp(-N(l))
    stays within 44 rho of the path link, keeps every site within 1/2 of I
    and solves the constraint, and that |det dV/dl_0| <= 1.957^3.
Double precision; the proofs in the paper are by hand.
Usage: python3 every_fibre.py
"""
import itertools
import numpy as np
from scipy.linalg import expm, logm

checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)


S = [np.array(m, complex) for m in ([[0, 1], [-1, 0]], [[0, 1j], [1j, 0]], [[1j, 0], [0, -1j]])]
B = [s / np.sqrt(2) for s in S]
I2 = np.eye(2, dtype=complex)


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
    return np.linalg.norm(v - I2)


def near(rng, s):
    z = rng.normal(size=3)
    return expm(alg(z / np.linalg.norm(z)) * s * rng.uniform())


ALPHA = 2 * np.arcsin(1 / (2 * np.sqrt(2)))
LAM = ALPHA / np.sin(ALPHA)
K7 = 7 * LAM / 8
OFFS = list(itertools.product((0, 1), repeat=3))   # transverse offsets a


def is_tree(a, mu):
    """u_a = (x_a, mu) with o_mu = 0 is a tree link iff o_nu = 0 for nu > mu."""
    others = [nu for nu in range(4) if nu != mu]
    return all(a[i] == 0 for i, nu in enumerate(others) if nu > mu)


def sites(g, u, l, w):
    """the sixteen v_x of one pair; u, l, w indexed by the offsets a."""
    vs = []
    for a in OFFS:
        vs += [g.conj().T @ u[a] @ l[a], g.conj().T @ l[a] @ w[a]]
    return vs


def average(g, u, l, w):
    vs = sites(g, u, l, w)
    m = max(delta(v) for v in vs)
    chi = min(1.0, max(0.0, 2 - 2 * m))
    N = sum(lg(v) for v in vs) / 16 if chi > 0 else 0 * I2
    return g @ expm(chi * N), vs


# (A) the aligned field
worst_v, worst_c = 0.0, 0.0
rng = np.random.default_rng(11)
for mu in range(4):
    for t in range(25):
        n = rng.integers(1, 5)
        c = [haar(rng) for _ in range(n)]
        w_end = {a: (I2 if is_tree(a, mu) else haar(rng)) for a in OFFS}
        u, l = [None] * n, [None] * n
        w_next = w_end
        for j in range(n - 1, -1, -1):
            g = c[j]
            u[j], l[j] = {}, {}
            for a in OFFS:
                if is_tree(a, mu):
                    u[j][a], l[j][a] = I2, g
                else:
                    l[j][a] = g @ w_next[a].conj().T
                    u[j][a] = g @ w_next[a] @ g.conj().T
            w_next = u[j]
        for j in range(n):
            w = u[j + 1] if j + 1 < n else w_end
            V, vs = average(l[j][(0, 0, 0)], u[j], l[j], w)
            worst_v = max(worst_v, max(delta(v) for v in vs))
            worst_c = max(worst_c, np.linalg.norm(V - c[j]))
            assert np.allclose(u[j][(0, 0, 0)], I2)
print(f"aligned fields: largest delta(v_x, I) {worst_v:.1e}, largest |V - c| {worst_c:.1e}")
check("the aligned field has v_x = I at every site and the prescribed averages",
      worst_v < 1e-12 and worst_c < 1e-12)

# (B) regular backgrounds with large links inside the blocks
worst_ratio, worst_site, max_det, conv, tested = 0.0, 0.0, 0.0, True, 0
mu = 0
for t in range(200):
    g0 = haar(rng)
    u = {a: (I2 if is_tree(a, mu) else haar(rng)) for a in OFFS}
    w = {a: (I2 if is_tree(a, mu) else haar(rng)) for a in OFFS}
    # private links with v_{x_a}, v_{x'_a} within 1/4 of I
    l = {}
    for a in OFFS:
        if a == (0, 0, 0):
            l[a] = g0
        else:
            l[a] = u[a].conj().T @ g0 @ near(rng, 0.12)
    for a in OFFS:
        if a != (0, 0, 0):
            # adjust w so that v_{x'_a} is also near I
            w[a] = l[a].conj().T @ g0 @ near(rng, 0.12) if not is_tree(a, mu) else I2
    vs0 = sites(g0, u, l, w)
    if max(delta(v) for v in vs0) > 0.25:
        continue
    tested += 1
    c, _ = average(g0, u, l, w)
    rho = 10 ** rng.uniform(-6, -3)
    def pert(x):
        return x @ near(rng, rho)
    u2 = {a: (I2 if is_tree(a, mu) else pert(u[a])) for a in OFFS}
    w2 = {a: (I2 if is_tree(a, mu) else pert(w[a])) for a in OFFS}
    l2 = {a: (pert(l[a]) if a != (0, 0, 0) else None) for a in OFFS}
    x = g0
    for it in range(300):
        l2[(0, 0, 0)] = x
        vs = sites(x, u2, l2, w2)
        worst_site = max(worst_site, max(delta(v) for v in vs))
        N = sum(lg(v) for v in vs) / 16
        x = c @ expm(-N)
        worst_ratio = max(worst_ratio, dist(x, g0) / rho)
    l2[(0, 0, 0)] = x
    V, vs = average(x, u2, l2, w2)
    conv &= np.linalg.norm(V - c) < 1e-10
    cols = []
    for k in range(3):
        e = 1e-7
        Vs = []
        for s in (e, -e):
            l3 = dict(l2); l3[(0, 0, 0)] = expm(s * B[k]) @ x
            Vs.append(average(l3[(0, 0, 0)], u2, l3, w2)[0])
        cols.append(coords(lg(Vs[0] @ V.conj().T) - lg(Vs[1] @ V.conj().T)) / (2 * e))
    max_det = max(max_det, abs(np.linalg.det(np.array(cols).T)))
print(f"regular backgrounds: {tested} tested, largest d(l, l_0)/rho {worst_ratio:.3f}, "
      f"largest site {worst_site:.3f}, largest |det| {max_det:.4f}")
check("the iteration stays within 44 rho of the path link and solves the constraint",
      tested > 100 and conv and worst_ratio <= 44)
check("the sites stay within 1/2 of I, so the cut-off equals 1", worst_site < 0.5)
check(f"|det dV/dl_0| <= 1.957^3 = {1.957 ** 3:.3f}", max_det <= 1.957 ** 3)
check(f"1/4 + 46/1076 < 1/2 and 7 Lambda/8 * 46 = {46 * K7:.3f} < 44", 0.25 + 46 / 1076 < 0.5 and 46 * K7 < 44)
P = 4 * 269
check(f"2.93 * 14 * 1076^3 = {2.93 * 14 * P ** 3:.3e} < 5.2e10", 2.93 * 14 * P ** 3 < 5.2e10)
check("94 * 2 sqrt 2 + 47^2 * 1e-3 < 269", 94 * 2 * np.sqrt(2) + 47 ** 2 * 1e-3 < 269)
print(f"{sum(checks)} of {len(checks)} checks passed")
