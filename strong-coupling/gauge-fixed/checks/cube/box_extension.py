"""Checks of the box extension lemma and the conditional density bound
(Lemma lem:boxextension and Proposition prop:conddensity of the forest
paper, Sec. ssec:largefields).

For boxes [0,l]^4 with l = 2, 3, 4 and random closed boundary data with
plaquette angles at most rho_s on the boundary sphere, this program

  (a) builds the real lift alpha on the boundary facet by facet in the
      order f1+, f2+, f2-, f3+, f3-, f4+, f4-, f1-, with the axial gauge on
      each facet and the McShane extension of the potential of the
      difference, and checks d alpha = theta, the agreement on earlier
      facets, and |alpha| <= 13118 l rho_s;
  (b) checks that the path metric of the union Gamma of earlier facets is
      at most three times the l^1 metric of the facet;
  (c) recovers the given torus links from alpha up to a gauge;
  (d) evaluates the radial Whitney extension by quadrature and checks that
      it reproduces alpha on the boundary links and that every plaquette of
      the box has angle at most C rho_s, C = sqrt(3)(24*13118+9);
  (e) checks the constants: the recursion M_8 = 13118, C, the Haar bound
      of a ball in SU(2), and the exponent of the density bound.

Usage: python3 box_extension.py
"""

import itertools
import math
import random
from collections import deque

import numpy as np

checks = []


def check(name, ok):
    checks.append((name, bool(ok)))
    print(("ok   " if ok else "FAIL ") + name)


def unit(m):
    v = [0, 0, 0, 0]
    v[m] = 1
    return tuple(v)


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def box(l):
    verts = list(itertools.product(range(l + 1), repeat=4))
    links = [(x, m) for x in verts for m in range(4) if x[m] < l]
    plaqs = [(x, m, n) for x in verts for m in range(4) for n in range(m + 1, 4)
             if x[m] < l and x[n] < l]
    return verts, links, plaqs


def on_bdry(x, l):
    return any(c in (0, l) for c in x)


def link_in_face(e, i, s, l):
    x, m = e
    return m != i and x[i] == s * l


def plaq_in_face(p, i, s, l):
    x, m, n = p
    return m != i and n != i and x[i] == s * l


def bdry_link(e, l):
    x, m = e
    return any(x[k] in (0, l) for k in range(4) if k != m)


def bdry_plaq(p, l):
    x, m, n = p
    return any(x[k] in (0, l) for k in range(4) if k not in (m, n))


def plaq_links(p):
    x, m, n = p
    return [((x, m), 1), ((add(x, unit(m)), n), 1),
            ((add(x, unit(n)), m), -1), ((x, n), -1)]


def curl(a, p):
    return sum(s * a[e] for e, s in plaq_links(p))


def wrap(t):
    return (t + math.pi) % (2 * math.pi) - math.pi


ORDER = [(0, 1), (1, 1), (1, 0), (2, 1), (2, 0), (3, 1), (3, 0), (0, 0)]


def axial(l, i, s, theta):
    d1, d2, d3 = [k for k in range(4) if k != i]
    a = {}
    face = [x for x in itertools.product(range(l + 1), repeat=4) if x[i] == s * l]
    for x in face:
        if x[d3] < l:
            a[(x, d3)] = 0.0
    for x in sorted(face, key=lambda y: y[d3]):
        if x[d2] < l:
            if x[d3] == 0:
                a[(x, d2)] = 0.0
            if x[d3] < l:
                a[(add(x, unit(d3)), d2)] = a[(x, d2)] - theta[(x, d2, d3)]
    for x in sorted(face, key=lambda y: (y[d3], y[d2])):
        if x[d1] < l:
            if x[d2] == 0 and x[d3] == 0:
                a[(x, d1)] = 0.0
            if x[d3] == 0 and x[d2] < l:
                a[(add(x, unit(d2)), d1)] = a[(x, d1)] - theta[(x, d1, d2)]
            if x[d3] < l:
                a[(add(x, unit(d3)), d1)] = a[(x, d1)] - theta[(x, d1, d3)]
    return a


def lift(l, theta):
    """Facet-by-facet lift; returns alpha and the list of M_k/(l rho)."""
    rho = max(abs(theta[p]) for p in theta)
    alpha = {}
    Ms = []
    gamma_ratio = 0.0
    for k, (i, s) in enumerate(ORDER):
        a0 = axial(l, i, s, theta)
        flinks = list(a0)
        fpl = [p for p in theta if plaq_in_face(p, i, s, l)]
        assert all(abs(curl(a0, p) - theta[p]) < 1e-12 for p in fpl)
        assert max(abs(v) for v in a0.values()) <= 2 * l * rho + 1e-12
        gl = [e for e in flinks if e in alpha]
        if not gl:
            for e in flinks:
                alpha[e] = a0[e]
        else:
            om = {e: alpha[e] - a0[e] for e in gl}
            adj = {}
            for (x, m) in gl:
                y = add(x, unit(m))
                adj.setdefault(x, []).append((y, om[(x, m)]))
                adj.setdefault(y, []).append((x, -om[(x, m)]))
            gv = list(adj)
            lam = {gv[0]: 0.0}
            dq = deque([gv[0]])
            while dq:
                x = dq.popleft()
                for y, w in adj[x]:
                    if y not in lam:
                        lam[y] = lam[x] + w
                        dq.append(y)
                    else:
                        assert abs(lam[y] - lam[x] - w) < 1e-9, "omega not exact"
            assert len(lam) == len(gv)
            L = 3 * max(abs(v) for v in om.values())
            # path metric of Gamma against the l^1 metric
            for x in gv:
                dist = {x: 0}
                dq = deque([x])
                while dq:
                    u = dq.popleft()
                    for v, _ in adj[u]:
                        if v not in dist:
                            dist[v] = dist[u] + 1
                            dq.append(v)
                for v, dv in dist.items():
                    d1 = sum(abs(a - b) for a, b in zip(x, v))
                    if d1:
                        gamma_ratio = max(gamma_ratio, dv / d1)
            fverts = {x for (x, m) in flinks} | {add(x, unit(m)) for (x, m) in flinks}
            lt = {x: min(lam[y] + L * sum(abs(a - b) for a, b in zip(x, y)) for y in gv)
                  for x in fverts}
            for (x, m) in flinks:
                val = a0[(x, m)] + lt[add(x, unit(m))] - lt[x]
                if (x, m) in alpha:
                    assert abs(val - alpha[(x, m)]) < 1e-9
                else:
                    alpha[(x, m)] = val
        Ms.append(max(abs(v) for v in alpha.values()) / (l * rho))
    return alpha, Ms, gamma_ratio


def whitney_extension(l, alpha, npts=12, nsub=24):
    """B_e = integral over e of (r/R0)^2 phi^* hat(alpha)."""
    c = np.full(4, l / 2.0)
    R0 = l / 2.0
    gx, gw = np.polynomial.legendre.leggauss(npts)
    ts, ws = [], []
    for j in range(nsub):
        ts.extend((j + (gx + 1) / 2) / nsub)
        ws.extend(gw / (2 * nsub))
    ts, ws = np.array(ts), np.array(ws)

    def form(y, i, s, v):
        """hat(alpha) at the point y of facet (i, s) applied to v."""
        dirs = [k for k in range(4) if k != i]
        u = np.array([y[k] for k in dirs])
        q = np.minimum(np.floor(u), l - 1).astype(int)
        f = u - q
        tot = 0.0
        for a, k in enumerate(dirs):
            others = [b for b in range(3) if b != a]
            val = 0.0
            for e0, e1 in itertools.product((0, 1), repeat=2):
                corner = q.copy()
                corner[others[0]] += e0
                corner[others[1]] += e1
                w = ((f[others[0]] if e0 else 1 - f[others[0]])
                     * (f[others[1]] if e1 else 1 - f[others[1]]))
                x = [0, 0, 0, 0]
                x[i] = s * l
                for b, kk in enumerate(dirs):
                    x[kk] = int(corner[b])
                val += w * alpha[(tuple(x), k)]
            tot += val * v[k]
        return tot

    B = {}
    for (x, m) in alpha_links(l):
        tot = 0.0
        for t, w in zip(ts, ws):
            p = np.array(x, dtype=float)
            p[m] += t
            z = p - c
            r = np.max(np.abs(z))
            if r < 1e-14:
                continue
            i = int(np.argmax(np.abs(z)))
            s = 1 if z[i] > 0 else 0
            phi = c + R0 * z / r
            phi[i] = s * l
            v = np.zeros(4)
            v[m] = 1.0
            grad = np.zeros(4)
            grad[i] = 1.0 if s else -1.0
            dphi = (R0 / r) * (v - z * (grad @ v) / r)
            tot += w * (r / R0) ** 2 * form(phi, i, s, dphi)
        B[(x, m)] = tot
    return B


def alpha_links(l):
    return box(l)[1]


def run(l, seed, ext=True):
    rng = random.Random(seed)
    verts, links, plaqs = box(l)
    bl = [e for e in links if bdry_link(e, l)]
    bp = [p for p in plaqs if bdry_plaq(p, l)]
    # torus data: a small curvature part, a large gauge part, reduced mod 2 pi
    small = {e: rng.uniform(-1, 1) for e in bl}
    mu = {x: rng.uniform(-math.pi, math.pi) for x in verts if on_bdry(x, l)}
    raw = {e: small[e] for e in bl}
    th0 = {p: curl(raw, p) for p in bp}
    rho_s = 0.3
    sc = rho_s / max(abs(v) for v in th0.values())
    given = {(x, m): wrap(sc * small[(x, m)] + mu[add(x, unit(m))] - mu[x])
             for (x, m) in bl}
    theta = {p: wrap(curl(given, p)) for p in bp}
    rho = max(abs(v) for v in theta.values())
    check(f"l={l} seed={seed}: boundary angles at most rho_s < pi/3",
          rho <= rho_s + 1e-12 and rho_s < math.pi / 3)
    alpha, Ms, gr = lift(l, theta)
    check(f"l={l} seed={seed}: d alpha = theta on every boundary plaquette",
          max(abs(curl(alpha, p) - theta[p]) for p in bp) < 1e-9)
    bound = [2]
    for _ in range(7):
        bound.append(3 * bound[-1] + 8)
    check(f"l={l} seed={seed}: M_k/(l rho) {[round(v, 2) for v in Ms]} within "
          f"{bound}", all(a <= b + 1e-9 for a, b in zip(Ms, bound)))
    check(f"l={l} seed={seed}: path metric of Gamma at most {gr:.3f} <= 3 "
          f"times the l^1 metric", gr <= 3 + 1e-12)
    # gauge: given - alpha is d mu mod 2 pi on the boundary
    gam = {e: wrap(given[e] - alpha[e]) for e in bl}
    check(f"l={l} seed={seed}: given - alpha is flat mod 2 pi",
          max(abs(wrap(curl(gam, p))) for p in bp) < 1e-9)
    if not ext:
        return
    B = whitney_extension(l, alpha)
    check(f"l={l} seed={seed}: extension reproduces alpha on boundary links",
          max(abs(B[e] - alpha[e]) for e in bl) < 1e-9)
    C = math.sqrt(3) * (24 * 13118 + 9)
    m = max(abs(curl(B, p)) for p in plaqs)
    check(f"l={l} seed={seed}: largest box plaquette angle {m / rho:.2f} rho_s "
          f"<= C rho_s, C = {C:.0f}", m <= C * rho)


print("(a)-(d) lift and extension")
for l in (2, 3, 4):
    for seed in range(3):
        run(l, seed, ext=(l <= 3 or seed == 0))

print("(e) constants")
M = 2
for _ in range(7):
    M = 3 * M + 8
check(f"M_8 = {M} = 13118", M == 13118)
C = math.sqrt(3) * (24 * 13118 + 9)
check(f"C = sqrt(3)(24*13118+9) = {C:.1f} < 545400", C < 545400)
# Haar measure of the ball of angle t in SU(2): (2/pi) int_0^t sin^2
for t in np.linspace(1e-3, math.pi / 2, 200):
    v = (t - math.sin(t) * math.cos(t)) / math.pi
    assert v >= 8 * t ** 3 / (3 * math.pi ** 3) - 1e-15
check("Haar ball of angle t <= pi/2 has measure >= 8 t^3/(3 pi^3)", True)
C0 = 0.5 + 2 / 3 * math.log(24 * math.pi ** 3)
for l in range(2, 40):
    nP = 6 * l * l * (l + 1) ** 2
    nL = 4 * l * (l + 1) ** 3
    assert nP <= 6 * (l + 1) ** 4 and nL <= 4 * (l + 1) ** 4
    # n_P + n_L log(24 pi^3 beta^{3/2}) <= 6(l+1)^4 (log beta + C0 + 1/2)
    for beta in (1.0, 10.0, 1e3, 1e8):
        lhs = nP + nL * math.log(24 * math.pi ** 3 * beta ** 1.5)
        rhs = 6 * (l + 1) ** 4 * (math.log(beta) + C0 + 0.5)
        assert lhs <= rhs + 1e-9
check("counts of plaquettes and links and the exponent of the bound", True)

npass = sum(ok for _, ok in checks)
print(f"{npass} checks passed, {len(checks) - npass} failed")
