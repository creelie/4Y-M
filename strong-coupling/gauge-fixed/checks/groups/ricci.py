"""Ricci constants of the compact simple groups in a faithful representation.

For a compact simple group G with a faithful unitary representation rho of
dimension n, the metric <X,Y> = -Tr(rho(X) rho(Y)) is bi-invariant, and the
Ricci curvature is Ric = c/4, where B = c * (trace form) and B is the Killing
form. On a Cartan element H, B(H,H) = sum over roots a of a(H)^2 and the trace
form is sum over weights m of rho (with multiplicity) of m(H)^2. Simplicity
makes the ratio independent of H; this is checked on several H.

The script works in exact rational arithmetic and prints, for every group,
Ric, the dual Coxeter number, the Dynkin index, and the ranges
n|beta| < 2 Ric/29 and beta_W = n^2 beta < 2 n Ric/29 of the forest-gauge
curvature theorem. Run: python3 ricci.py
"""
from fractions import Fraction as Fr
from itertools import combinations, product
import random

def vec(*xs):
    return tuple(Fr(x) for x in xs)

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def unit(n, i, s=1):
    return tuple(Fr(s) if k == i else Fr(0) for k in range(n))

def add(a, b, s=1):
    return tuple(x + s * y for x, y in zip(a, b))

def roots_A(n):  # A_n in R^{n+1}
    m = n + 1
    return [add(unit(m, i), unit(m, j), -1) for i in range(m) for j in range(m) if i != j]

def weights_A(n):
    m = n + 1
    avg = tuple(Fr(1, m) for _ in range(m))
    return [add(unit(m, i), avg, -1) for i in range(m)]

def pm_pairs(n):
    out = []
    for i, j in combinations(range(n), 2):
        for s, t in product((1, -1), repeat=2):
            out.append(add(unit(n, i, s), unit(n, j, t)))
    return out

def roots_B(n):
    return pm_pairs(n) + [unit(n, i, s) for i in range(n) for s in (1, -1)]

def weights_B(n):  # vector of SO(2n+1)
    return [unit(n, i, s) for i in range(n) for s in (1, -1)] + [tuple(Fr(0) for _ in range(n))]

def roots_C(n):
    return pm_pairs(n) + [unit(n, i, 2 * s) for i in range(n) for s in (1, -1)]

def weights_C(n):  # fundamental of Sp(n), dimension 2n
    return [unit(n, i, s) for i in range(n) for s in (1, -1)]

def roots_D(n):
    return pm_pairs(n)

def weights_D(n):  # vector of SO(2n)
    return [unit(n, i, s) for i in range(n) for s in (1, -1)]

def roots_G2():
    e = [unit(3, i) for i in range(3)]
    long_ = [add(e[i], e[j], -1) for i in range(3) for j in range(3) if i != j]
    short = []
    for i in range(3):
        j, k = [x for x in range(3) if x != i]
        v = tuple((2 * e[i][t] - e[j][t] - e[k][t]) / 3 for t in range(3))
        short += [v, tuple(-x for x in v)]
    return long_ + short, short

def weights_G2():
    _, short = roots_G2()
    return short + [tuple(Fr(0) for _ in range(3))]

def roots_F4():
    short = [unit(4, i, s) for i in range(4) for s in (1, -1)]
    short += [tuple(Fr(s, 2) for s in signs) for signs in product((1, -1), repeat=4)]
    return pm_pairs(4) + short, short

def weights_F4():
    _, short = roots_F4()
    return short + [tuple(Fr(0) for _ in range(4))] * 2

def roots_E8():
    r = pm_pairs(8)
    for signs in product((1, -1), repeat=8):
        if signs.count(-1) % 2 == 0:
            r.append(tuple(Fr(s, 2) for s in signs))
    return r

def project(v, basis):
    """Orthogonal projection of v onto the complement of span(basis) (basis orthogonalised)."""
    ob = []
    for b in basis:
        w = b
        for o in ob:
            w = add(w, o, -dot(w, o) / dot(o, o))
        ob.append(w)
    for o in ob:
        v = add(v, o, -dot(v, o) / dot(o, o))
    return v

def E7_data():
    R = roots_E8()
    r = R[0]
    roots = [a for a in R if dot(a, r) == 0]
    wts = [project(a, [r]) for a in R if dot(a, r) == 1]
    return roots, wts, [r]

def E6_data():
    R = roots_E8()
    r1 = R[0]
    r2 = next(a for a in R if dot(a, r1) == -1)  # r1, r2 simple roots of an A2
    roots = [a for a in R if dot(a, r1) == 0 and dot(a, r2) == 0]
    wts = [project(a, [r1, r2]) for a in R if (dot(a, r1), dot(a, r2)) == (1, 0)]
    return roots, wts, [r1, r2]

def ratio(roots, wts, dim, perp=()):
    vals = set()
    rng = random.Random(7)
    for _ in range(4):
        H = tuple(Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(dim))
        if perp:
            H = project(H, list(perp))
        if all(x == 0 for x in H):
            continue
        num = sum(dot(a, H) ** 2 for a in roots)
        den = sum(dot(m, H) ** 2 for m in wts)
        vals.add(num / den)
    assert len(vals) == 1, vals
    return vals.pop()

def sumzero(H):
    m = len(H)
    s = sum(H) / m
    return tuple(x - s for x in H)

def main():
    rows = []
    for n in range(1, 8):
        rows.append((f"SU({n+1})", n + 1, roots_A(n), weights_A(n), n + 1, "A"))
    for n in range(2, 7):
        rows.append((f"Sp({n})", 2 * n, roots_C(n), weights_C(n), n, ""))
    for n in range(1, 6):
        rows.append((f"SO({2*n+1})", 2 * n + 1, roots_B(n), weights_B(n), n, ""))
    for n in range(3, 7):
        rows.append((f"SO({2*n})", 2 * n, roots_D(n), weights_D(n), n, ""))
    g2r, _ = roots_G2()
    rows.append(("G2", 7, g2r, weights_G2(), 3, "A"))
    f4r, _ = roots_F4()
    rows.append(("F4", 26, f4r, weights_F4(), 4, ""))
    e6r, e6w, e6p = E6_data()
    e7r, e7w, e7p = E7_data()
    e8r = roots_E8()
    print(f"{'group':8} {'n':>4} {'#roots':>6} {'#wts':>5} {'Ric':>8} {'h':>4} {'index':>6} {'n|beta| <':>10} {'beta_W <':>10}")
    out = []
    for name, nrep, roots, wts, dim, kind in rows + [("E6", 27, e6r, e6w, 8, "E6"), ("E7", 56, e7r, e7w, 8, "E7"),
                                                     ("E8", 248, e8r, e8r + [tuple(Fr(0) for _ in range(8))] * 8, 8, "")]:
        assert len(wts) == nrep, (name, len(wts))
        perp = ()
        if kind == "E6":
            perp = e6p
        elif kind == "E7":
            perp = e7p
        if kind in ("A",):
            # Cartan of A_n and of G2 lies in the sum-zero hyperplane
            vals = set()
            rng = random.Random(3)
            for _ in range(4):
                H = sumzero(tuple(Fr(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(dim)))
                vals.add(sum(dot(a, H) ** 2 for a in roots) / sum(dot(m, H) ** 2 for m in wts))
            assert len(vals) == 1, (name, vals)
            c = vals.pop()
        else:
            c = ratio(roots, wts, dim, perp)
        ric = c / 4
        # dual Coxeter number h and Dynkin index l: l_adj = 2h and c = l_adj / l_rho
        # h from the highest root theta: 2h = (1/2) sum_a <a, theta^vee>^2
        long2 = max(dot(a, a) for a in roots)
        theta = max(roots, key=lambda a: (dot(a, a), a))
        thv = tuple(2 * x / dot(theta, theta) for x in theta)
        two_h = sum(dot(a, thv) ** 2 for a in roots) / 2
        index = sum(dot(m, thv) ** 2 for m in wts) / 2
        assert two_h / index == c, (name, two_h, index, c)
        h = two_h / 2
        out.append((name, nrep, ric, h, index))
        print(f"{name:8} {nrep:>4} {len(roots):>6} {len(wts):>5} {str(ric):>8} {str(h):>4} {str(index):>6} {str(2*ric/29):>10} {str(2*nrep*ric/29):>10}")
    # checks against the closed forms
    known = {"G2": (Fr(1), 4, 2), "F4": (Fr(3, 4), 9, 6), "E6": (Fr(1), 12, 6), "E7": (Fr(3, 4), 18, 12), "E8": (Fr(1, 4), 30, 60)}
    for name, nrep, ric, h, index in out:
        if name.startswith("SU("):
            N = nrep
            assert ric == Fr(N, 2) and h == N and index == 1
        elif name.startswith("Sp("):
            k = nrep // 2
            assert ric == Fr(k + 1, 2) and h == k + 1 and index == 1
        elif name.startswith("SO("):
            N = nrep
            assert ric == Fr(N - 2, 4), (name, ric)
            if N >= 5:
                assert h == N - 2 and index == 2, (name, h, index)
        else:
            assert (ric, h, index) == known[name], (name, ric, h, index)
    print("all closed forms confirmed")

if __name__ == "__main__":
    main()
