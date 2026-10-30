"""Checks of the nonabelian filling lemma (Lemma lem:nonabelianfilling of
the forest paper, Sec. ssec:largefields).

For boxes [0,l]^4 and SU(2) boundary fields that are gauge transforms, by
arbitrary (large) vertex gauges, of fields with small plaquettes, this
program

  (a) builds the coordinate-path gauge of each facet and checks that its
      links have angle at most 2 l rho_s (l = 2, 3, 4);
  (b) glues the facet gauges in the order f1+, f2+, f2-, f3+, f3-, f4+,
      f4-, f1- through principal logarithms and McShane extensions, and
      checks the agreement on the earlier facets and the bounds
      M_k <= C_k l rho_s (l = 2, 3, 4);
  (c) checks the curl bound |(dA)_p| <= rho_s + 8 M^2 for the logarithms
      of the gauged boundary links;
  (d) for l = 2, 3 evaluates the radial cubical extension B, integrates its
      transports along every link of the box piece by piece between the
      cone interfaces, undoes the gauge, and checks that the boundary
      links are reproduced and that every plaquette of the box has angle
      at most K rho_s, K = 42 C_8 + 20;
  (e) checks the constants C_k, K, eps_0 and the numerical inequalities of
      the proof in exact rational arithmetic.

The boundary fields satisfy l^2 rho_s <= eps_0, so rho_s is of order
1e-21; the computation uses mpmath with 60 digits. These are numerical
checks of the construction, not interval certificates; the proof in the
paper covers every l >= 2.

Usage: python3 nonabelian_filling.py
"""

import itertools
import random
from fractions import Fraction as Fr

import mpmath as mp

mp.mp.dps = 60

checks = []


def check(name, ok):
    checks.append((name, bool(ok)))
    print(("ok   " if ok else "FAIL ") + name)


# ---------------------------------------------------------------- quaternions

def qmul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0)


def qinv(a):
    return (a[0], -a[1], -a[2], -a[3])


ONE = (mp.mpf(1), mp.mpf(0), mp.mpf(0), mp.mpf(0))


def vnorm(v):
    return mp.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)


def angle(q):
    return mp.atan2(vnorm(q[1:]), q[0])


def qexp(v):
    t = vnorm(v)
    if t == 0:
        return ONE
    s = mp.sin(t) / t
    return (mp.cos(t), v[0] * s, v[1] * s, v[2] * s)


def qlog(q):
    t = angle(q)
    n = vnorm(q[1:])
    if n == 0:
        return (mp.mpf(0),) * 3
    return tuple(c * t / n for c in q[1:])


def vadd(u, w):
    return tuple(a + b for a, b in zip(u, w))


def vscale(c, u):
    return tuple(c * a for a in u)


def haar(rng):
    v = [mp.mpf(rng.gauss(0, 1)) for _ in range(4)]
    n = mp.sqrt(sum(x * x for x in v))
    return tuple(x / n for x in v)


def small(rng, amp):
    return tuple(mp.mpf(rng.uniform(-1, 1)) * amp for _ in range(3))


# ----------------------------------------------------------------- constants

def constants():
    C = [None, 2]
    for _ in range(2, 9):
        C.append(11 * C[-1] + 24)
    return C


C = constants()
C8 = C[8]
K = 42 * C8 + 20
EPS0 = Fr(1, 100 * C8 * C8)
SQRT3_UP = Fr(17320509, 10 ** 7)  # > sqrt(3)


def exact_checks():
    check("C_8 = 85743550", C8 == 85743550)
    check("K = 42 C_8 + 20 = 3601229120", K == 3601229120)
    check("eps_0 = 1/735195636660250000",
          EPS0 == Fr(1, 735195636660250000))
    check("sqrt(3) < 1.7320509", SQRT3_UP ** 2 > 3)
    check("6 sqrt(3) <= 11", 6 * SQRT3_UP <= 11)
    check("12 sqrt(3) + 2 <= 24", 12 * SQRT3_UP + 2 <= 24)
    check("9 (C_8 + 2) eps_0 < 1/2", 9 * (C8 + 2) * EPS0 < Fr(1, 2))
    check("C_8 eps_0 / 2 < 1/2 (M < 1/2 for l >= 2)",
          C8 * EPS0 / 2 < Fr(1, 2))
    check("C_8^2 eps_0 = 1/100 (M^2 <= rho_s/100)",
          C8 * C8 * EPS0 == Fr(1, 100))
    lhs = 24 * SQRT3_UP * C8 + 9 * SQRT3_UP + (72 * SQRT3_UP + 216) / 100
    check("24 sqrt3 C_8 + 9 sqrt3 + (72 sqrt3 + 216)/100 < 42 C_8 + 20",
          lhs < K)
    check("(1/2)/sin(1/2) < 2", mp.mpf("0.5") / mp.sin(mp.mpf("0.5")) < 2)
    ok = all(2 * mp.asin(x / 2) <= x * (1 + x * x)
             for x in [mp.mpf(k) / 1000 for k in range(1, 1001)])
    check("2 arcsin(x/2) <= x (1 + x^2) on (0,1]", ok)
    check("(4 h 3 sqrt3 M)^2 / 2 = 216 M^2 h^2", Fr(1, 2) * 16 * 9 * 3 == 216)


# ----------------------------------------------------------------- lattice

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


def bdry_link(x, m, l):
    y = add(x, unit(m))
    return any(i != m and x[i] in (0, l) for i in range(4)) and on_bdry(y, l)


def bdry_plaq(x, m, n, l):
    return any(i not in (m, n) and x[i] in (0, l) for i in range(4))


def holonomy(U, x, m, n):
    a = U[(x, m)]
    b = U[(add(x, unit(m)), n)]
    c = qinv(U[(add(x, unit(n)), m)])
    d = qinv(U[(x, n)])
    return qmul(qmul(a, b), qmul(c, d))


ORDER = [(0, 1), (1, 1), (1, 0), (2, 1), (2, 0), (3, 1), (3, 0), (0, 0)]


def facet_verts(l, i, s):
    return [x for x in itertools.product(range(l + 1), repeat=4)
            if x[i] == s * l]


def facet_links(l, i, s):
    return [(x, m) for x in facet_verts(l, i, s) for m in range(4)
            if m != i and x[m] < l]


def path_hol(eta, l, i, s, x):
    """Holonomy of the path from the facet corner to x, moving first in the
    lowest free direction, then the next, then the last."""
    free = [m for m in range(4) if m != i]
    y = [0, 0, 0, 0]
    y[i] = s * l
    h = ONE
    for m in free:
        for _ in range(x[m]):
            h = qmul(h, eta[(tuple(y), m)])
            y[m] += 1
    return h


def gauged(eta, g, x, m):
    return qmul(qmul(g[x], eta[(x, m)]), qinv(g[add(x, unit(m))]))


def d1(x, y):
    return sum(abs(a - b) for a, b in zip(x, y))


def boundary_gauge(eta, l, rho):
    """Facet gauges glued as in the proof. Returns g and the list of
    observed maxima per stage."""
    g = {}
    done_links = []
    stats = []
    for k, (i, s) in enumerate(ORDER, start=1):
        V = facet_verts(l, i, s)
        g0 = {x: path_hol(eta, l, i, s, x) for x in V}
        loc = max(angle(gauged(eta, g0, x, m)) for (x, m) in facet_links(l, i, s))
        check(f"l={l} facet {k}: coordinate-path gauge <= 2 l rho_s",
              loc <= 2 * l * rho * (1 + mp.mpf(10) ** -30))
        Gam = [x for x in V if x in g]
        if k == 1:
            for x in V:
                g[x] = g0[x]
        else:
            t = {x: qmul(g[x], qinv(g0[x])) for x in Gam}
            x0 = Gam[0]
            a = {x: qlog(qmul(qinv(t[x0]), t[x])) for x in Gam}
            L = 6 * (C[k - 1] * l * rho + 2 * l * rho)
            lip = max((vnorm(vadd(a[x], vscale(-1, a[y]))) / d1(x, y)
                       for x in Gam for y in Gam if x != y), default=0)
            check(f"l={l} facet {k}: log of transition is L-Lipschitz",
                  lip <= L * (1 + mp.mpf(10) ** -30))
            ta = {}
            for x in V:
                ta[x] = tuple(min(a[y][j] + L * d1(x, y) for y in Gam)
                              for j in range(3))
            agree = max(vnorm(vadd(ta[x], vscale(-1, a[x]))) for x in Gam)
            check(f"l={l} facet {k}: McShane extension agrees on Gamma",
                  agree < mp.mpf(10) ** -50 * max(rho, mp.mpf(10) ** -40))
            for x in V:
                if x not in g:
                    g[x] = qmul(qmul(t[x0], qexp(ta[x])), g0[x])
        done_links += facet_links(l, i, s)
        Mk = max(angle(gauged(eta, g, x, m)) for (x, m) in done_links)
        stats.append(Mk / (l * rho))
        check(f"l={l} after facet {k}: links <= C_{k} l rho_s "
              f"(observed {float(Mk / (l * rho)):.3g} l rho_s)",
              Mk <= C[k] * l * rho * (1 + mp.mpf(10) ** -30))
    return g, stats


def make_boundary(l, rng, amp):
    verts, links, plaqs = box(l)
    gg = {x: haar(rng) for x in verts if on_bdry(x, l)}
    eta = {}
    for (x, m) in links:
        if bdry_link(x, m, l):
            y = add(x, unit(m))
            eta[(x, m)] = qmul(qmul(gg[x], qexp(small(rng, amp))), qinv(gg[y]))
    rho = max(angle(holonomy(eta, x, m, n)) for (x, m, n) in plaqs
              if bdry_plaq(x, m, n, l))
    return eta, rho


# ----------------------------------------------------------- the extension

def make_B(A, l):
    c = mp.mpf(l) / 2
    R0 = mp.mpf(l) / 2

    def Ahat(i, base, y):
        """Whitney form on the unit cube of facet i with lower corner base,
        at the boundary point y; returns {direction: coefficient}."""
        free = [m for m in range(4) if m != i]
        t = {m: y[m] - base[m] for m in free}
        out = {}
        for j in free:
            o1, o2 = [m for m in free if m != j]
            acc = (mp.mpf(0),) * 3
            for b1 in (0, 1):
                for b2 in (0, 1):
                    e = list(base)
                    e[o1] += b1
                    e[o2] += b2
                    lam = ((t[o1] if b1 else 1 - t[o1]) *
                           (t[o2] if b2 else 1 - t[o2]))
                    acc = vadd(acc, vscale(lam, A[(tuple(e), j)]))
            out[j] = acc
        return out

    def B(cell, x, m):
        """B(e_m) at x, with the formula of the cone cell = (i, sgn, base)."""
        i, sgn, base = cell
        r = sgn * (x[i] - c)
        if r == 0:
            return (mp.mpf(0),) * 3
        y = [c + R0 * (xj - c) / r for xj in x]
        dphi = [R0 * (1 if j == m else 0) / r for j in range(4)]
        if m == i:
            dphi = [dphi[j] - R0 * (x[j] - c) * sgn / r ** 2 for j in range(4)]
        a = Ahat(i, base, y)
        out = (mp.mpf(0),) * 3
        for j, coef in a.items():
            out = vadd(out, vscale(dphi[j], coef))
        return vscale((r / R0) ** 2, out)

    def cell_at(x):
        d = [abs(x[j] - c) for j in range(4)]
        i = max(range(4), key=lambda j: d[j])
        sgn = 1 if x[i] > c else -1
        r = d[i]
        y = [c + R0 * (xj - c) / r for xj in x]
        base = []
        for j in range(4):
            if j == i:
                base.append(l if sgn > 0 else 0)
            else:
                base.append(min(int(mp.floor(y[j])), l - 1))
        return (i, sgn, tuple(base))

    def breakpoints(x0, m):
        cand = set()
        xm = mp.mpf(x0[m])
        for i in range(4):
            if i == m:
                continue
            di = abs(x0[i] - c)
            for sg in (1, -1):
                cand.add(c + sg * di - xm)
            for k in range(l + 1):
                cand.add(c + (k - c) * di / R0 - xm)
        for j in range(4):
            if j == m or x0[j] == c:
                continue
            for k in range(l + 1):
                if k == c:
                    continue
                rt = R0 * (x0[j] - c) / (k - c)
                if rt > 0:
                    for sg in (1, -1):
                        cand.add(c + sg * rt - xm)
        pts = sorted(s for s in cand if 0 < s < 1)
        out = [mp.mpf(0)]
        for s in pts:
            if s - out[-1] > mp.mpf(10) ** -40:
                out.append(s)
        if 1 - out[-1] > mp.mpf(10) ** -40:
            out.append(mp.mpf(1))
        else:
            out[-1] = mp.mpf(1)
        return out

    def transport(x0, m, steps=12):
        bp = breakpoints(x0, m)
        T = ONE
        for s0, s1 in zip(bp, bp[1:]):
            mid = [mp.mpf(v) for v in x0]
            mid[m] += (s0 + s1) / 2
            if all(abs(mid[j] - c) == 0 for j in range(4)):
                continue
            cell = cell_at(mid)
            h = (s1 - s0) / steps

            def f(s):
                p = [mp.mpf(v) for v in x0]
                p[m] += s
                return (mp.mpf(0),) + B(cell, p, m)

            s = s0
            for _ in range(steps):
                k1 = qmul(T, f(s))
                T2 = tuple(a + h / 2 * b for a, b in zip(T, k1))
                k2 = qmul(T2, f(s + h / 2))
                T3 = tuple(a + h / 2 * b for a, b in zip(T, k2))
                k3 = qmul(T3, f(s + h / 2))
                T4 = tuple(a + h * b for a, b in zip(T, k3))
                k4 = qmul(T4, f(s + h))
                T = tuple(a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
                          for a, b1, b2, b3, b4 in zip(T, k1, k2, k3, k4))
                s += h
        n = mp.sqrt(sum(a * a for a in T))
        return tuple(a / n for a in T)

    return transport


def full_filling(l, rng, amp):
    verts, links, plaqs = box(l)
    eta, rho = make_boundary(l, rng, amp)
    check(f"l={l} full: hypothesis l^2 rho_s <= eps_0 "
          f"(rho_s = {mp.nstr(rho, 3)})", l * l * rho <= mp.mpf(EPS0.numerator) / EPS0.denominator)
    g, _ = boundary_gauge(eta, l, rho)
    A = {}
    for (x, m) in links:
        if bdry_link(x, m, l):
            A[(x, m)] = qlog(gauged(eta, g, x, m))
    M = max(vnorm(a) for a in A.values())
    check(f"l={l} full: |A| <= C_8 l rho_s", M <= C8 * l * rho)
    curl = 0
    for (x, m, n) in plaqs:
        if bdry_plaq(x, m, n, l):
            v = vadd(vadd(A[(x, m)], A[(add(x, unit(m)), n)]),
                     vscale(-1, vadd(A[(add(x, unit(n)), m)], A[(x, n)])))
            curl = max(curl, vnorm(v))
    check(f"l={l} full: |dA| <= rho_s + 8 M^2", curl <= rho + 8 * M * M)
    transport = make_B(A, l)
    V = {(x, m): transport(x, m) for (x, m) in links}
    err = max(angle(qmul(qinv(V[e]), qexp(A[e]))) for e in A)
    check(f"l={l} full: boundary transports equal exp(A_e) "
          f"(error {mp.nstr(err, 3)})", err <= rho * mp.mpf(10) ** -12)
    gext = {x: g.get(x, ONE) for x in verts}
    U = {(x, m): qmul(qmul(qinv(gext[x]), V[(x, m)]), gext[add(x, unit(m))])
         for (x, m) in links}
    berr = max(angle(qmul(qinv(U[e]), eta[e])) for e in eta)
    check(f"l={l} full: filled field equals eta on boundary links",
          berr <= rho * mp.mpf(10) ** -12)
    worst = max(angle(holonomy(U, x, m, n)) for (x, m, n) in plaqs)
    check(f"l={l} full: every plaquette angle <= K rho_s "
          f"(observed {float(worst / rho):.4g} rho_s)", worst <= K * rho)
    inner = [p for p in plaqs if not bdry_plaq(*p, l)]
    check(f"l={l} full: {len(inner)} inner plaquettes checked", len(inner) > 0)


def main():
    exact_checks()
    rng = random.Random(20261030)
    for l in (2, 3, 4):
        for trial in range(2):
            amp = mp.mpf(EPS0.numerator) / EPS0.denominator / (8 * l * l)
            eta, rho = make_boundary(l, rng, amp)
            check(f"l={l} trial {trial}: hypothesis l^2 rho_s <= eps_0",
                  l * l * rho <= mp.mpf(EPS0.numerator) / EPS0.denominator)
            big = max(angle(e) for e in eta.values())
            check(f"l={l} trial {trial}: boundary links are not small "
                  f"(max angle {float(big):.3g})", big > 1)
            boundary_gauge(eta, l, rho)
    for l in (2, 3):
        amp = mp.mpf(EPS0.numerator) / EPS0.denominator / (8 * l * l)
        full_filling(l, rng, amp)
    n_ok = sum(ok for _, ok in checks)
    print(f"{n_ok} of {len(checks)} checks passed")
    raise SystemExit(0 if n_ok == len(checks) else 1)


if __name__ == "__main__":
    main()
