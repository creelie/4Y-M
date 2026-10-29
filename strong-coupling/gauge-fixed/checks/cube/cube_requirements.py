"""Finite counts behind the subsection on what the cube estimate requires.

The line field sigma of the residue table has period two. The script checks:
  (a) the table: every plaquette contains exactly one fixed link and every
      active link has twelve active partners;
  (b) the eight lines modulo two, and the colouring of the theorem on the
      cube estimate: residues 0010, 0011, 0100, 0110, 1001, 1011, 1100, 1101
      form four of the eight pairs, and exactly 32 of the 48 active classes
      join residues of different colour (32 is the maximum over all 2^8
      colourings of the pairs);
  (c) the census of (a_c, b_c) over the 48 classes in the boundary count
      a_c m^3 - b_c m^2, and the count itself against direct enumeration of
      boxes of side n = 2, 4, 6, 8, 10;
  (d) the slow line: for each of the 48 line types, the four transverse
      plaquettes at each link have two active links each, off the line, not
      twists, and pairwise distinct over a stretch of the line;
  (e) the exact arithmetic of the proofs;
  (f) the proposition on the cube size: every active link of the cube is the
      holonomy of a disc of at most m(m+1) plaquettes whose other boundary
      links lie in F or outside the cube, each plaquette lies in at most
      4n^2+2n discs (n = 2, 4, 6, 8), and the lattice Stokes inequality on
      random configurations (n = 2, 4);
  (g) the flux configuration of the remark on loops, in exact angles;
  (h) the arithmetic of the lower bound on eps_V;
  (i) the Gaussian benchmark: K_c = alpha_n(t) N_n for every residue class,
      alpha_n(0) = 1, and the bound for t > 0, in exact arithmetic;
  (j) the joint-cube transport route: a legal configuration where the weighted
      curvature is negative once beta_W > 1/3;
  (k) the four colours of the SU(2) proposition: the sets X_1, ..., X_4 are
      unions of two pairs, every active link joins different sets, X_3 u X_4
      is the colour set of (b), and the tetrahedral quaternions z_k have
      <z_k, z_l> = -1/3 for k != l, in exact arithmetic;
  (l) the four colours for every group: the quotients z_k z_l^{-1} have real
      part -1/3, their images in SU(3) (2x2 block), SU(4) (q tensor I_2) and
      SO(3) have equal Hilbert-Schmidt distance to the identity, and for SU(3)
      the diameter squared is 9 against the Haar mean 6 of the squared distance;
  (m) the cube size for SU(N): every element of SU(N) is exp X with X traceless
      and eigenvalue angles at most 3 pi; the Jacobian of exp in Hilbert-Schmidt
      coordinates is the product of (sin(phi/2)/(phi/2))^2 over pairs of angles
      (finite differences, N = 2, 3, 4); the lower bound
      (11/12)^(N(N-1)) (r/(3 pi sqrt N))^(N^2-1) for the Haar measure of the
      ball of radius r <= 1 against Monte Carlo (N = 2, 3); and the arithmetic
      of L_N(b);
  (n) the proposition on positive reflected norms in a box: the coefficients
      of the SU(2) characters chi_j in exp((b/2) Re chi) are positive for
      2j <= 4 at b = 0.5 and 2j <= 10 at b = 2, 6 (Weyl integration formula);
  (o) the remarks on the flowed ratio: the free massless covariance with
      spatial cutoff Lambda at time t and at coincidence (radial quadrature),
      and the bounds of the proposition on shifted renormalisation on random
      discrete spectral measures (tail, moments, Lipschitz, Jensen, slab tail);
  (w) two-dimensional SU(N): the holonomy identity and the expansion of the
      character coefficients;
  (x) the proposition on block averages of the free field: the exact
      conditional variance E Var(F | block averages) of a quadratic packet
      against 2 sigma^4 sum alpha_B^2, and Var F against 2 d^2 sum f^2, on
      periodic boxes in two and three dimensions;
  (y) the proposition on perturbed actions: on random finite metric spaces,
      W1 of two measures tilted by the same weight e^u is at most
      e^Omega (1 + 2 Delta l) times W1 of the measures (l the Lipschitz
      constant, Omega the oscillation of u, Delta the diameter), and W1 of one
      measure tilted by e^u and e^u' is at most 2 Delta sup|u - u'|, by
      linear programming; and the constants 48, 96 of the proposition.
Run: python3 cube_requirements.py
"""
from collections import Counter
from fractions import Fraction as Fr
from itertools import combinations, product

TABLE = """0000 0 0 2 3 1
1000 0 0 3 1 2
0100 2 3 0 2 1
1100 3 2 0 1 3
0010 3 1 2 0 3
1010 1 3 1 0 2
0110 2 1 3 2 0
1110 1 2 1 3 0
0001 1 2 1 3 0
1001 2 1 3 2 0
0101 1 3 1 0 2
1101 3 1 2 0 3
0011 3 2 0 1 3
1011 2 3 0 2 1
0111 0 0 3 1 2
1111 0 0 2 3 1"""

SIG, NEXT = {}, {}
for row in TABLE.splitlines():
    f = row.split()
    r = tuple(int(c) for c in f[0])
    SIG[r] = int(f[1])
    NEXT[r] = [int(c) for c in f[2:]]

def res(x):
    return tuple(t % 2 for t in x)

def sigma(x):
    return SIG[res(x)]

def move(x, mu, s=1):
    return tuple(v + (s if k == mu else 0) for k, v in enumerate(x))

def active(link):
    return sigma(link[0]) != link[1]

def plaquette(z, mu, nu):
    return [(z, mu), (move(z, mu), nu), (move(z, nu), mu), (z, nu)]

def plaquettes_of(link):
    y, mu = link
    out = []
    for nu in range(4):
        if nu != mu:
            for s in (1, -1):
                z = y if s == 1 else move(y, nu, -1)
                out.append((nu, plaquette(z, mu, nu)))
    return out

def partners(link):
    return {l for _, p in plaquettes_of(link) for l in p if l != link and active(l)}

checks = 0
def ok(cond, msg):
    global checks
    if not cond:
        raise SystemExit("FAILED: " + msg)
    checks += 1

# (a) the table
for r in SIG:
    for mu in range(4):
        ok(sigma(move(r, mu)) == NEXT[r][mu], "table column")
    for mu, nu in combinations(range(4), 2):
        ok(sum(not active(l) for l in plaquette(r, mu, nu)) == 1, "one fixed link")
CLASSES = [(r, mu) for r in SIG for mu in range(4) if SIG[r] != mu]
ok(len(CLASSES) == 48, "48 classes")
for c in CLASSES:
    ok(len(partners(c)) == 12, "twelve partners")

# (b) lines modulo two and the colouring
PAIR = {}
for r in SIG:
    r2 = res(move(r, SIG[r]))
    ok(SIG[r2] == SIG[r], "straight lines")
    PAIR[r] = frozenset([r, r2])
PAIRS = sorted(set(PAIR.values()), key=sorted)
ok(len(PAIRS) == 8, "eight lines modulo two")
W = {tuple(int(c) for c in s) for s in
     "0010 0011 0100 0110 1001 1011 1100 1101".split()}
ok(all(PAIR[r] <= W for r in W), "W is a union of pairs")
def cut(colour):
    return sum(colour(r) != colour(res(move(r, mu))) for r, mu in CLASSES)
ok(cut(lambda r: r in W) == 32, "32 classes cut")
best = max(cut(lambda r, b=b: b[PAIRS.index(PAIR[r])])
           for b in product((0, 1), repeat=8))
ok(best == 32, "32 is the maximum")

# (c) the boundary count
census = Counter()
for r, mu in CLASSES:
    a = b = 0
    for base, _ in partners((r, mu)):
        d = tuple(u - v for u, v in zip(base, r))
        q = sum((r[l] == 1 and d[l] == 1) or (r[l] == 0 and d[l] == -1)
                for l in range(4))
        a += q
        b += q == 2
    census[(a, b)] += 1
ok(sorted(census.items()) == [((2, 0), 3), ((3, 0), 6), ((4, 0), 6),
   ((5, 0), 6), ((6, 0), 6), ((7, 0), 3), ((7, 1), 3), ((8, 1), 6),
   ((9, 1), 3), ((9, 2), 3), ((10, 2), 3)], "census")
for (a, b) in census:
    ok(a <= 10 and 8 - a + b >= 0, "dominated by (10, 2)")
for n in (2, 4, 6, 8, 10):
    m = n // 2
    inside = lambda y: all(0 <= t < n for t in y)
    P = Counter()
    for y in product(range(-1, n + 1), repeat=4):
        if inside(y):
            continue
        for mu in range(4):
            if active((y, mu)):
                k = sum(inside(l[0]) for l in partners((y, mu)))
                if k:
                    P[(res(y), mu)] += k
    formula = {}
    for r, mu in CLASSES:
        a = b = 0
        for base, _ in partners((r, mu)):
            d = tuple(u - v for u, v in zip(base, r))
            q = sum((r[l] == 1 and d[l] == 1) or (r[l] == 0 and d[l] == -1)
                    for l in range(4))
            a += q
            b += q == 2
        formula[(r, mu)] = a * m**3 - b * m**2
    ok(all(P[c] == formula[c] for c in CLASSES), f"box count n={n}")
    ok(max(P.values()) == 10 * m**3 - 2 * m**2, f"maximum n={n}")

# (d) the slow line
for x0 in SIG:
    s = SIG[x0]
    for mu in range(4):
        if mu == s:
            continue
        line = [(tuple(v + (k if i == s else 0) for i, v in enumerate(x0)), mu)
                for k in range(30)]
        on = set(line)
        ok(all(active(l) for l in line), "line active")
        twists, others = set(), []
        for link in line[2:-2]:
            for nu, p in plaquettes_of(link):
                rest = [l for l in p if l != link and active(l)]
                ok(sum(not active(l) for l in p) == 1, "fixed link")
                if nu == s:
                    nb = [l for l in rest if l in on]
                    tw = [l for l in rest if l not in on]
                    ok(len(nb) == 1 and len(tw) == 1 and tw[0][1] == s, "twist")
                    twists.add(tw[0])
                else:
                    ok(len(rest) == 2 and not set(rest) & on, "transverse")
                    others += rest
        ok(len(others) == len(set(others)), "transverse links distinct")
        ok(not set(others) & twists, "transverse links are not twists")

# (e) exact arithmetic
ok(Fr(17, 64) + 3 * Fr(17, 64)**2 + Fr(1, 64) == Fr(2019, 4096), "barrier")
ok(Fr(2019, 4096) < Fr(3, 4), "barrier contradiction")
ok(Fr(1, 8) * (2 + Fr(1, 8)) == Fr(17, 64), "eta")
ok(Fr(8, 3 * 4**3) == Fr(1, 24), "ball volume bound with pi < 4")
bw, n = Fr(3, 10), 16
ok(bw / (1 - 3 * bw) * (Fr(5, n) - Fr(2, n * n)) == Fr(117, 128), "alpha_16")
r = [(1, 0, 0, 0), (-1, 0, 0, 0), (0, 1, 0, 0), (0, -1, 0, 0), (0, 0, 1, 0),
     (0, 0, -1, 0), (0, 0, 1, 0), (0, -1, 0, 0), (0, 1, 0, 0), (-1, 0, 0, 0)]
ok(all(r[k] == r[(10 - k) % 10] for k in range(10)), "palindrome")
for k in range(10):
    h = [Fr(2 * r[k][i] - r[k - 1][i] - r[(k + 1) % 10][i], 2) for i in range(4)]
    ok(sum(t * t for t in h) <= 4, "field at most 2b")

# (f) Proposition on the cube size: links as loops, multiplicity, Stokes
import math, random

def qm(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    return (a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
            a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0)

def qc(a):
    return (a[0], -a[1], -a[2], -a[3])

def chord(a):  # delta(1, a) = sqrt(2 - Re Tr a), Re Tr a = 2 a_0
    return math.sqrt(max(0.0, 2 - 2 * a[0]))

ONE = (1.0, 0.0, 0.0, 0.0)

def disc(e, n):
    """The disc R(e) of the proof: its plaquettes (base, plane) and its boundary
    as a list of (link, orientation) starting at x."""
    x, mu = e
    m = n // 2
    y = move(x, mu)
    s_, t_ = sigma(x), sigma(y)
    def reach(z, d):
        up, down = n - z[d], z[d] + 1
        return (up, 1) if up <= m else (down, -1)
    s, eps = reach(x, s_)
    t, eps2 = reach(y, t_)
    xp = move(x, s_, s * eps)
    path, pts = [], [x]
    def walk(d, k, sgn):
        for _ in range(k):
            z = pts[-1]
            if sgn == 1:
                path.append(((z, d), 1))
            else:
                path.append(((move(z, d, -1), d), -1))
            pts.append(move(z, d, sgn))
    walk(s_, s, eps)          # along the line of F through x
    walk(mu, 1, 1)            # to x' + e_mu
    walk(t_, t, eps2)         # along e_tau
    walk(s_, s, -eps)         # along e_sigma
    walk(t_, t, -eps2)        # back along the line of F through y
    walk(mu, 1, -1)           # over e
    plaq = set()
    for a in range(s):
        z = move(x, s_, eps * a if eps == 1 else -(a + 1))
        plaq.add((z, frozenset((s_, mu))))
        for c in range(t):
            w = move(move(y, s_, eps * a if eps == 1 else -(a + 1)), t_,
                     c if eps2 == 1 else -(c + 1))
            plaq.add((w, frozenset((s_, t_))))
    return s, t, path, pts, plaq

for n in (2, 4, 6, 8):
    m = n // 2
    inside = lambda z: all(0 <= v < n for v in z)
    V = [(x, mu) for x in product(range(n), repeat=4) for mu in range(4)
         if active((x, mu))]
    Vset = set(V)
    ok(len(V) == 48 * m**4, f"|V| = 48 N_n, n={n}")
    mult = Counter()
    for e in V:
        x, mu = e
        ok(sigma(move(x, mu)) not in (mu, sigma(x)), "sigma, tau, mu distinct")
        s, t, path, pts, plaq = disc(e, n)
        ok(s <= m and t <= m and len(plaq) == s * (t + 1) <= m * (m + 1), "disc size")
        ok(pts[-1] == x, "closed loop")
        others = [l for l, _ in path[:-1]]
        ok(path[-1][0] == e and e not in others, "e once")
        ok(all(not active(l) or l not in Vset for l in others), "boundary off V")
        for p in plaq:
            mult[p] += 1
    ok(max(mult.values()) <= 4 * n * n + 2 * n, f"multiplicity n={n}")
    if n <= 4:
        rng = random.Random(n)
        U = {}
        for e in V:
            g = [rng.gauss(0, 1) for _ in range(4)]
            r = math.sqrt(sum(v * v for v in g))
            U[e] = tuple(v / r for v in g)
        link = lambda l: U.get(l, ONE)
        def pl(z, a, b):
            return qm(qm(link((z, a)), link((move(z, a), b))),
                      qm(qc(link((move(z, b), a))), qc(link((z, b)))))
        for e in V:
            s, t, path, pts, plaq = disc(e, n)
            h = ONE
            for l, o in path:
                h = qm(h, link(l) if o == 1 else qc(link(l)))
            ok(abs(chord(h) - chord(U[e])) < 1e-9, "holonomy is Q_e")
            rhs = sum(chord(pl(z, *sorted(P))) for z, P in plaq)
            ok(chord(h) <= rhs + 1e-9, "lattice Stokes")
for n in range(2, 60, 2):
    m = n // 2
    ok(m * (m + 1) * (4 * n * n + 2 * n) * 2 == n * n * (n + 2) * (2 * n + 1),
       "m(m+1)(4n^2+2n) = n^2(n+2)(2n+1)/2")
    ok((n + 2) * (2 * n + 1) <= 5 * n * n, "(n+2)(2n+1) <= 5n^2")

# (g) the flux configuration of the remark: angles in units of 2 pi / L^2
for L in (4, 6, 8, 10, 12):
    N = L * L
    phi = 2
    def ang(z, d):
        z = tuple(v % L for v in z)
        if d == 1:
            return (phi * z[0]) % N
        if d == 0 and z[0] == L - 1:
            return (-phi * L * z[1]) % N
        return 0
    def pa(z, a, b):
        return (ang(z, a) + ang(move(z, a), b) - ang(move(z, b), a) - ang(z, b)) % N
    for z in product(range(L), repeat=2):
        for a, b in combinations(range(4), 2):
            zz = (z[0], z[1], 0, 0)
            ok(pa(zz, a, b) == (phi if (a, b) == (0, 1) else 0), "flux plaquettes")
    h, k = 0, L // 2
    for i in range(k):
        h += ang((i, 0, 0, 0), 0) + ang((k, i, 0, 0), 1) \
            - ang((i, k, 0, 0), 0) - ang((0, i, 0, 0), 1)
    ok(h % N == N // 2, "holonomy -1 on the square")
    ok(phi * k * k == N // 2, "disc bound pi")

# (h) the lower bound on eps_V: 2 e^(1/2) pi^2 (2/9)^3 < 1/2, with e < 2.7183
#     (so e^(1/2) < 1.6488) and pi < 3.1416
ok(Fr(16488, 10000) ** 2 > Fr(27183, 10000), "sqrt(e) < 1.6488")
ok(2 * Fr(16488, 10000) * Fr(31416, 10000) ** 2 * Fr(2, 9) ** 3 < Fr(1, 2), "small ball")

# (i) the Gaussian benchmark: K_c = alpha_n(t) N_n for every class, alpha_n(0) = 1,
#     and the lower bound on 1 - alpha_n(t), in exact arithmetic
def solve(M, rhs):
    N = len(M)
    M = [row[:] + [r] for row, r in zip(M, rhs)]
    for c in range(N):
        p = next(r for r in range(c, N) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        inv = 1 / M[c][c]
        M[c] = [v * inv for v in M[c]]
        for r in range(N):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [M[r][N] for r in range(N)]

for d, n in ((1, 2), (1, 6), (2, 2), (2, 4), (2, 6), (3, 2), (3, 4), (4, 2)):
    sites = list(product(range(n), repeat=d))
    index = {x: i for i, x in enumerate(sites)}
    def nbrs(x):
        for i in range(d):
            for s in (1, -1):
                yield tuple(v + (s if k == i else 0) for k, v in enumerate(x))
    ext = sorted({y for x in sites for y in nbrs(x) if y not in index})
    for t in (Fr(0), Fr(1, 3), Fr(2)):
        A = [[Fr(0)] * len(sites) for _ in sites]
        for x in sites:
            A[index[x]][index[x]] = 2 * d + t
            for y in nbrs(x):
                if y in index:
                    A[index[x]][index[y]] = Fr(-1)
        Gone = solve(A, [Fr(1)] * len(sites))      # G 1 (A is symmetric)
        alpha = 1 - t * sum(Gone) / n**d
        K = Counter()
        for j in ext:
            C = sum(Gone[index[x]] for x in nbrs(j) if x in index)  # 1^T G B_j
            K[tuple(v % 2 for v in j)] += C
        ok(len(K) == 2**d and all(k == alpha * Fr(n // 2) ** d for k in K.values()),
           f"K_c = alpha N_n, d={d} n={n} t={t}")
        if t == 0:
            ok(alpha == 1, "massless alpha = 1")
        else:
            ok(1 - alpha >= Fr(5, 6) ** d * t * n * n / (t * n * n + 10 * d), "mass bound")
            ok(all(g > 0 for g in Gone), "G 1 > 0")

# (j) the joint-cube transport route: at the legal configuration U_j = -1, all
#     other links 1, the six plaquettes through j have holonomy -1, and along a
#     geodesic of SU(2) = S^3 through U_j the action -beta_W sum <Q_p, 1> has
#     second derivative -6 beta_W, against Ricci curvature 2 of the unit S^3:
#     the weighted curvature is negative once beta_W > 1/3, at every cube size
for X in ((0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (0, 0.6, 0, 0.8)):
    for beta in (0.2, 1 / 3, 0.5, 2.3):
        def act(s):
            u = qm((-1, 0, 0, 0), (math.cos(s),) + tuple(math.sin(s) * c for c in X[1:]))
            return -beta * 6 * u[0]                  # six plaquettes, <Q_p, 1> = Re part
        h = 1e-4
        second = (act(h) - 2 * act(0) + act(-h)) / h**2
        ok(abs(second + 6 * beta) < 1e-5, "second derivative -6 beta_W")
        ok((2 + second < 0) == (beta > 1 / 3 + 1e-12), "curvature sign at 1/3")
ok(Fr(4, 29) < Fr(1, 3), "Bakry-Emery range 4/29 below the threshold 1/3")

# (k) four colours for SU(2)
XS = [{tuple(int(c) for c in w) for w in ws.split()} for ws in (
    "0000 1000 0111 1111", "0001 0101 1010 1110",
    "0010 0011 1100 1101", "0100 0110 1001 1011")]
ok(sorted(r for X in XS for r in X) == sorted(SIG), "four sets partition the residues")
colour = {r: k for k, X in enumerate(XS) for r in X}
for X in XS:
    ok(all(res(move(r, SIG[r])) in X for r in X), "union of pairs")
for r in SIG:
    for mu in range(4):
        if mu != SIG[r]:
            ok(colour[r] != colour[res(move(r, mu))], "active link joins different sets")
ok(XS[2] | XS[3] == {tuple(int(c) for c in w) for w in
   "0010 0011 0100 0110 1001 1011 1100 1101".split()}, "X3 u X4 is the colour set of (b)")
Z = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
for k, l in combinations(range(4), 2):
    ok(Fr(sum(a * b for a, b in zip(Z[k], Z[l])), 3) == Fr(-1, 3), "<z_k, z_l> = -1/3")
ok(2 - 2 * Fr(-1, 3) == Fr(8, 3), "d_*^2 = 8/3")

# (l) four colours for every group
import numpy as np
def su2(q):
    a, b, c, d = q
    return np.array([[a + 1j * b, c + 1j * d], [-c + 1j * d, a - 1j * b]])
ZQ = [(0,) + tuple(v / math.sqrt(3) for v in z) for z in Z]
def so3(q):
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
embed = {"SU(3)": lambda q: np.block([[su2(q), np.zeros((2, 1))], [np.zeros((1, 2)), np.eye(1)]]),
         "SU(4)": lambda q: np.kron(su2(q), np.eye(2)),
         "SO(3)": so3}
for name, phi in embed.items():
    dist = []
    for k in range(4):
        for l in range(4):
            if k != l:
                quo = qm(ZQ[k], qc(ZQ[l]))
                ok(abs(quo[0] + 1 / 3) < 1e-12, "real part -1/3")
                M = phi(ZQ[k]) @ np.conj(phi(ZQ[l])).T
                dist.append(np.linalg.norm(M - np.eye(M.shape[0])))
    ok(max(dist) - min(dist) < 1e-12 and min(dist) > 0.5, "equal positive distances " + name)
best = min(math.cos(a) + math.cos(b) + math.cos(a + b)
           for a in np.linspace(0, 2 * math.pi, 721) for b in np.linspace(0, 2 * math.pi, 721))
ok(abs(best + 1.5) < 1e-9, "min Re tr on SU(3) is -3/2")
ok(6 + 2 * Fr(3, 2) == 9 and Fr(6, 3) == Fr(2, 3) * 3, "SU(3): D^2 = 9, 2d/D = (2/3) D")

# (m) the cube size for SU(N)
rng = np.random.default_rng(20261028)
def haar_su(N, k):
    Zm = (rng.standard_normal((k, N, N)) + 1j * rng.standard_normal((k, N, N))) / math.sqrt(2)
    Qm, Rm = np.linalg.qr(Zm)
    dg = np.diagonal(Rm, axis1=1, axis2=2)
    Qm = Qm * (dg / abs(dg))[:, None, :]
    return Qm / (np.linalg.det(Qm) ** (1 / N))[:, None, None]
def expm_h(X):  # X skew-Hermitian
    w, U = np.linalg.eigh(-1j * X)
    return (U * np.exp(1j * w)) @ np.conj(U).T
def basis_su(N):  # Hilbert-Schmidt orthonormal basis of su(N)
    B = []
    for j in range(N):
        for k in range(j + 1, N):
            E = np.zeros((N, N), complex); E[j, k] = 1; E[k, j] = -1; B.append(E / math.sqrt(2))
            E = np.zeros((N, N), complex); E[j, k] = 1j; E[k, j] = 1j; B.append(E / math.sqrt(2))
    for j in range(1, N):
        E = np.zeros((N, N), complex)
        for k in range(j): E[k, k] = 1j
        E[j, j] = -1j * j; B.append(E / np.linalg.norm(E))
    return B
for N in (2, 3, 4):
    for q in haar_su(N, 20):
        w, U = np.linalg.eig(q)
        th = np.angle(w)
        th = np.where(th >= math.pi, th - 2 * math.pi, th)
        kk = int(round(th.sum() / (2 * math.pi)))
        order = np.argsort(th)
        if kk > 0: th[order[-kk:]] -= 2 * math.pi
        if kk < 0: th[order[:-kk]] += 2 * math.pi
        ok(abs(th.sum()) < 1e-9 and max(abs(th)) <= 3 * math.pi, "angles traceless, at most 3 pi")
        ok(np.linalg.norm((U * np.exp(1j * th)) @ np.linalg.inv(U) - q) < 1e-8, "exp X = q")
    B = basis_su(N)
    G = np.array([[np.trace(np.conj(a).T @ b).real for b in B] for a in B])
    ok(np.allclose(G, np.eye(N * N - 1)), "orthonormal basis of su(N)")
    for _ in range(5):
        c = rng.standard_normal(N * N - 1); c *= rng.uniform(0.1, 2.5) / np.linalg.norm(c)
        X = sum(a * b for a, b in zip(c, B))
        E0inv = np.conj(expm_h(X)).T
        h = 1e-6
        J = np.array([[np.trace(np.conj(Bb).T @ (E0inv @ (expm_h(X + h * Ba) - expm_h(X - h * Ba)) / (2 * h))).real
                       for Ba in B] for Bb in B])
        ang = np.linalg.eigvalsh(-1j * X)
        prod = 1.0
        for a, b in combinations(ang, 2):
            x = (a - b) / 2
            prod *= (math.sin(x) / x) ** 2
        ok(abs(abs(np.linalg.det(J)) - prod) < 1e-5, f"Jacobian of exp, N = {N}")
        if np.linalg.norm(X) <= 1:
            ok(prod >= (11 / 12) ** (N * (N - 1)), "Jacobian at least (11/12)^(N(N-1))")
for N, r in ((2, 0.6), (3, 0.9)):
    dim = N * N - 1
    lower = (11 / 12) ** (N * (N - 1)) * (r / (3 * math.pi * math.sqrt(N))) ** dim
    hits = total = 0
    for _ in range(20):
        Qs = haar_su(N, 50000)
        d2 = np.sum(abs(Qs - np.eye(N)) ** 2, axis=(1, 2))
        hits += int(np.sum(d2 <= r * r)); total += len(d2)
    ok(hits / total > lower, f"ball measure above the bound, N = {N}")
for N in (2, 3, 4, 5):
    for bb in (N / 9, 1, 6, 100):
        r = math.sqrt(N / (9 * bb))
        lhs = (N * N - 1) / 2 * math.log(9 * bb / N) + (N * N - 1) * math.log(3 * math.pi * math.sqrt(N))
        ok(abs(lhs - (N * N - 1) / 2 * math.log(81 * math.pi ** 2 * bb)) < 1e-9, "L_N arithmetic")
        ok(abs(math.log(3 * math.pi * math.sqrt(N) / r) - 0.5 * math.log(81 * math.pi ** 2 * bb)) < 1e-9, "radius")
ok(all(math.sin(x) / x >= 1 - x * x / 6 for x in np.linspace(1e-6, 1 / math.sqrt(2), 1000)), "sin x / x >= 1 - x^2/6")

# (n) positive character coefficients of the Wilson weight for SU(2)
th = np.linspace(0, math.pi, 20001)
wts = (2 / math.pi) * np.sin(th) ** 2
trap = getattr(np, "trapezoid", None) or getattr(np, "trapz")
for bb, top in ((0.5, 4), (2.0, 10), (6.0, 10)):
    dens = np.exp(bb * np.cos(th))  # (b/N) Re chi = (b/2) 2 cos theta
    for twoj in range(top + 1):
        chi = np.sin((twoj + 1) * th[1:-1]) / np.sin(th[1:-1])
        c = trap(dens[1:-1] * chi * wts[1:-1], th[1:-1])
        ok(c > 0, f"positive character coefficient, b = {bb}, 2j = {twoj}")

# (o) the flowed ratio
for Lam, t in ((1.0, 0.5), (10.0, 0.3), (100.0, 1.0)):
    p = np.linspace(0, Lam, 400001)
    num = trap(p * np.exp(-p * t), p) / (4 * math.pi ** 2)
    ok(abs(num - (1 - (1 + Lam * t) * math.exp(-Lam * t)) / (4 * math.pi ** 2 * t * t)) < 1e-8 * max(1, num), "C(t,0)")
    ok(abs(trap(p, p) / (4 * math.pi ** 2) - Lam ** 2 / (8 * math.pi ** 2)) < 1e-8 * Lam ** 2, "C(0,0)")
ok(2 * (1 - (1 + 1e4) * math.exp(-1e4)) / 1e8 < 1e-7, "contact ratio tends to zero")
for trial in range(200):
    E = rng.exponential(rng.uniform(0.5, 20), 30); w = rng.dirichlet(np.ones(30))
    sig = rng.uniform(0.2, 2); R = float(np.sum(w * np.exp(-sig * E)))
    wp = w * np.exp(-sig * E) / R
    for Lam in (0.5, 2, 8):
        ok(np.sum(wp[E >= Lam]) <= np.exp(-sig * Lam) / R + 1e-12, "tail")
    for m in (1, 2, 3):
        ok(np.sum(wp * E ** m) <= (m / (math.e * sig)) ** m / R + 1e-9, "moments")
    ts = rng.uniform(0, 3, 2)
    ok(abs(np.sum(wp * (np.exp(-ts[0] * E) - np.exp(-ts[1] * E)))) <= abs(ts[0] - ts[1]) / (math.e * sig * R) + 1e-12, "Lipschitz")
    a = rng.uniform(0.001, 0.05); d = float(np.sum(w * (1 - np.exp(-a * E)))); C = d / a
    for n in (1, 5, 40):
        if C * a <= 0.5:
            ok(np.sum(w * np.exp(-a * n * E)) >= math.exp(-2 * C * a * n) - 1e-12, "Jensen")
    for Lam in (0.5, 2, 8):
        ok(np.sum(w[E >= Lam]) <= C * (a + 1 / Lam) + 1e-12, "slab tail")

# (p) packet ratio against slice correlators: D(2 a n1 + s)/D(2 a n0) <= C(s)/C(0)
# <= D(2 a n0 + s)/D(2 a n1), log-convexity, and the extension from one pair
def lse(x):
    return float(np.logaddexp.reduce(np.asarray(x, dtype=float)))
for trial in range(400):
    a = 10 ** rng.uniform(-3, -1); K = int(rng.integers(1, 40))
    ll = -a * 10 ** rng.uniform(-1, 3.5, K); lnu = np.log(rng.exponential(size=K))
    t0 = rng.uniform(0.05, 1); t1 = t0 + rng.uniform(0.01, 1)
    n0 = int(np.ceil(t0 / a)); n1 = int(np.floor(t1 / a))
    if n1 <= n0:
        continue
    n = np.arange(n0, n1 + 1); lw = np.log(rng.uniform(1e-3, 1, len(n)) * a)
    lW = np.array([lse(lw + l * n) for l in ll])
    lD = lambda j: lse(lnu + ll * j)
    lC = lambda m: lse(lnu + 2 * lW + ll * m)
    m = int(rng.uniform(0.05, 2) / a); r = lC(m) - lC(0)
    ok(lD(2 * n1 + m) - lD(2 * n0) - 1e-9 <= r <= lD(2 * n0 + m) - lD(2 * n1) + 1e-9, "slice bounds")
    j = int(rng.integers(1, 200)); ok(2 * lD(j) <= lD(j - 1) + lD(j + 1) + 1e-9, "log-convex D")
    ok(2 * lC(j) <= lC(j - 1) + lC(j + 1) + 1e-9, "log-convex C")
    i0, i1 = sorted(rng.choice(300, 2, replace=False)); lc = lD(i1) - lD(i0)
    al = int(rng.integers(i0, 400)); be = al + int(rng.integers(1, 400))
    ok(lD(be) - lD(al) >= lc * max(1, (be - al) / (i1 - i0)) - 1e-9, "extension from one pair")

# (q) the free lattice field: D_k(tau')/D_k(tau) converges to the continuum ratio
def slice_corr(a, tau, tf=0.05, N=128):
    L = np.pi / a; P = min(L, 8.0); ps = np.linspace(-P, P, N, endpoint=False) + P / N
    px, py, pz = np.meshgrid(ps, ps, ps, indexing='ij')
    ph2 = sum((2 / a * np.sin(a * q / 2)) ** 2 for q in (px, py, pz))
    w = np.arccosh(1 + a * a * ph2 / 2) / a
    hh = np.zeros_like(px)
    for kx in (-1, 0, 1):
        for ky in (-1, 0, 1):
            for kz in (-1, 0, 1):
                hh += np.exp(-((px + 2 * L * kx) ** 2 + (py + 2 * L * ky) ** 2 + (pz + 2 * L * kz) ** 2) / 4)
    return np.sum(hh ** 2 * np.exp(-2 * tf * ph2) * np.exp(-tau * w) / (2 * np.sinh(a * w) / a))
x = np.linspace(1e-6, 30, 300001)
cont = [trap(x * x * np.exp(-x * x / 2 - 2 * 0.05 * x * x - t * x) / (2 * x), x) for t in (0.5, 1.0)]
rc = cont[1] / cont[0]
ok(abs(rc - 0.63365) < 1e-4, "continuum free ratio")
for a in (0.1, 0.0125, 0.001):
    ok(abs(slice_corr(a, 1.0) / slice_corr(a, 0.5) - rc) < 2e-3, f"lattice free ratio, a = {a}")
th = np.linspace(1e-9, 6, 100001)
ok(np.all(np.arccosh(1 + th) >= math.acosh(7) / math.sqrt(12) * np.sqrt(2 * th) - 1e-12), "arccosh bound")

# (r) the three-state chain: two-point e^{-mt}, four-point q1 q3 (q2/eps + 1 - q2)
from fractions import Fraction as Fr
for eps in (Fr(1, 2), Fr(1, 10), Fr(1, 1000)):
    pi = {-1: eps / 2, 0: 1 - eps, 1: eps / 2}
    S = (-1, 0, 1)
    def P(q, s, t):
        return q * (s == t) + (1 - q) * pi[t]
    for q1, q2, q3 in ((Fr(1, 2), Fr(1, 2), Fr(1, 2)), (Fr(1, 3), Fr(3, 4), Fr(2, 5)), (Fr(9, 10), Fr(1, 7), Fr(1, 1))):
        ok(all(pi[s] * P(q1, s, t) == pi[t] * P(q1, t, s) for s in S for t in S), "detailed balance")
        two = sum(pi[s] * s * P(q1, s, t) * t for s in S for t in S) / eps
        ok(two == q1, "two-point")
        four = sum(pi[s0] * s0 * P(q1, s0, s1) * s1 * P(q2, s1, s2) * s2 * P(q3, s2, s3) * s3
                   for s0 in S for s1 in S for s2 in S for s3 in S) / eps ** 2
        ok(four == q1 * q3 * (q2 / eps + 1 - q2), "four-point")
        ok(four - q1 * q3 - 2 * q1 * q2 ** 2 * q3 == q1 * q3 * (q2 / eps - q2 - 2 * q2 ** 2), "truncated four-point")

# (s) the range of the slice-ratio extension, and the later-pair counterexample
D = lambda t: (Fr(1, 2) ** t + Fr(1, 4) ** t) / 2
c = D(3) / D(1)
ok(c == Fr(3, 16) and D(2) / D(1) == Fr(5, 12), "slice values")
ok((D(2) / D(1)) ** 2 < c, "power bound fails inside the pair")
ok(D(2) / D(1) >= c, "monotone bound inside the pair")
for al, be in ((1, 4), (2, 5), (1, 7), (3, 9)):
    ok(D(be) / D(al) >= c ** Fr(max(2, be - al), 2) if (be - al) % 2 == 0 else
       (D(be) / D(al)) ** 2 >= c ** max(2, be - al), "power bound to the right")
for k in (10, 20, 40):
    Dk = lambda t: (1 - Fr(1, 4 ** k)) * Fr(1, 2 ** (k * t)) + Fr(1, 4 ** k) * Fr(1, 2 ** t)
    ok(abs(Dk(4) / Dk(3) - Fr(1, 2)) < Fr(1, 2 ** (k - 2)), "later pair ratio near 1/2")
    ok(Dk(2) / Dk(1) < Fr(2, 2 ** k), "earlier pair ratio tends to zero")
# (t) the amplitude scales cancel: r^n (c r^2 V)^(-n/2) V^(n-1) = c^(-n/2) V^(n/2-1)
for r in (Fr(1, 7), Fr(3, 2), Fr(10 ** -6)):
    for V in (Fr(1, 3), Fr(5, 1)):
        for n in (2, 4, 6):
            lhs = r ** n * V ** (n - 1) / (r * r * V) ** (n // 2)
            ok(lhs == V ** (n // 2 - 1), "amplitude cancels")

# (u) tree bounds force a closing gap: C(s)/C(0) <= e^{-m s/a} for spectral measures on
# Omega^perp, and the product-form tail ratio 2 e^{-mu R'}/(1 + e^{-mu}) <= 2 e^{-mu R}
for trial in range(300):
    a = 10 ** rng.uniform(-3, -1); m = rng.uniform(0, 3); K = int(rng.integers(1, 30))
    lam = np.exp(-m) * rng.uniform(0, 1, K); w = rng.exponential(size=K); w /= w.sum()
    s = a * int(rng.integers(1, 400))
    ok(float(np.sum(w * lam ** round(s / a))) <= math.exp(-m * s / a) + 1e-15, "gap bound")
for mu in (Fr(1, 3), Fr(1, 1), Fr(5, 2)):
    q = math.exp(-float(mu))
    for R in (0.5, 3.0, 7.25, 20.0):
        Rp = math.ceil(R); J = 4000
        tail = sum(q ** abs(j) for j in range(-J, J + 1) if abs(j) >= Rp)
        full = sum(q ** abs(j) for j in range(-J, J + 1))
        ok(abs(tail / full - 2 * q ** Rp / (1 + q)) < 1e-12, "tail fraction")
        ok(tail / full <= 2 * math.exp(-float(mu) * R) + 1e-15, "tail bound")
        ok((1 + q) / (1 - q) >= 2 / float(mu), "coth bound")

# (v) energy control from the order-two tree bounds, uniform local estimates, blocking:
# finite differences, 1 - e^{-x} >= x/(1+x), log(1-x) >= -2x, the SU(2) one-link
# constants, the Bessel majorant, and |Tr(T^2 M)| <= 1 for unit T and unitary M
from scipy.special import iv, ive
from scipy.integrate import quad
for m in (1, 2, 3):
    for a in (0.1, 0.01):
        xs = np.linspace(0.3, 2.0, 50)
        f = lambda x: np.exp(-1 / (x * (2.3 - x))) if 0 < x < 2.3 else 0.0
        D = [sum((-1) ** j * math.comb(m, j) * f(x - j * a) for j in range(m + 1)) for x in xs]
        h = 1e-3; g = np.linspace(0.01, 2.29, 4000)
        der = [max(abs(sum((-1) ** j * math.comb(m, j) * f(y + (m / 2 - j) * h) for j in range(m + 1)) / h ** m) for y in g)]
        ok(max(abs(d) for d in D) <= a ** m * der[0] * 1.05 + 1e-12, "finite difference")
for x in np.linspace(0, 50, 501):
    ok(1 - math.exp(-x) >= x / (1 + x) - 1e-15, "1 - e^-x")
for x in np.linspace(0, 0.5, 101):
    ok(math.log1p(-x) >= -2 * x - 1e-15, "log(1-x)")
ok(math.e ** 2 * math.sqrt(math.pi / 2) < 14, "Z_T bound")
ok(4 * 224 * 336 < 2 ** 19 and 1792 * 1 <= 2 ** 19, "variance constants")
ok(math.exp(-0.5) / 8 >= 1 / 16 and math.exp(-1.5) / 4 >= 1 / 24, "interval masses")
for r in np.concatenate([np.linspace(0, 10, 41), [30.0, 100.0, 600.0]]):
    n = 1 + r
    br = [k / n for k in (1, 4, 16, 64) if k / n < 2]
    Z = quad(lambda t: math.exp(-r * t) * math.sqrt(t * (2 - t)), 0, 2, limit=200, points=br)[0]
    ok(Z <= 14 * n ** -1.5, "normalisation")
    pI = quad(lambda t: math.exp(-r * t) * math.sqrt(t * (2 - t)), 1 / (4 * n), 1 / (2 * n))[0] / Z
    pJ = quad(lambda t: math.exp(-r * t) * math.sqrt(t * (2 - t)), 1 / n, 1.5 / n)[0] / Z
    ok(pI >= 1 / 224 and pJ >= 1 / 336, "P(I), P(J)")
    m1 = quad(lambda t: t * math.exp(-r * t) * math.sqrt(t * (2 - t)), 0, 2, limit=200, points=br)[0] / Z
    m2 = quad(lambda t: t * t * math.exp(-r * t) * math.sqrt(t * (2 - t)), 0, 2, limit=200, points=br)[0] / Z
    ok(m2 - m1 ** 2 >= 1 / (2 ** 19 * n * n), "Var X")
    EX = 1 - m1
    if r > 0:
        ok(abs(EX - iv(2, r) / iv(1, r)) < 1e-7, "E X = I2/I1")
        ok(abs(3 * EX - r * (2 * m1 - m2)) < 1e-8 * (1 + r), "sphere identity")
        ok(EX <= (math.sqrt(9 + 4 * r * r) - 3) / (2 * r) + 1e-12 <= 2 * r / (2 * r + 3) + 1e-12, "mean bound")
        ok(all(iv(l + 1, r) <= iv(2, r) * (1 + 1e-12) for l in range(1, 8)), "Bessel order")
for b in (0.5, 1.0, 5.0, 50.0):
    ok(2 * 6 * b / (2 * 6 * b + 3) <= 4 * b / (4 * b + 1) + 1e-15 and 3 / (12 * b + 3) >= 1 / (4 * b + 1) - 1e-15, "b bounds")
    for B in (4, 10, 100):
        ok(math.exp(3 * b) * (4 * b / (4 * b + 1)) ** (B - 3) <= math.exp(3 * b - (B - 3) / (4 * b + 1)) * (1 + 1e-12), "blocking bound")
for N in (2, 3, 4):
    for trial in range(40):
        X = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)); T = X - X.conj().T
        T -= np.trace(T) / N * np.eye(N); T /= math.sqrt(-np.trace(T @ T).real)
        M = np.linalg.qr(rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N)))[0]
        ok(abs(np.trace(T @ T @ M)) <= 1 + 1e-12, "|Tr(T^2 M)| <= 1")

# (w) two-dimensional SU(N): the holonomy identity, chat_lambda = 1 - N C2/beta + O(beta^-2)
# (SU(2) by Bessel functions, SU(3) fundamental and adjoint by Weyl integration), and the
# energies -(ell/a^2) log chat_F -> g^2 ell C2/2 at beta = 2N/(g^2 a^2)
def haar(N):
    Z = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z); Q = Q * (np.diag(R) / abs(np.diag(R)))
    return Q / np.linalg.det(Q) ** (1 / N)
for N in (2, 3):
    for L in (2, 3, 5):
        V = [haar(N) for _ in range(L)]; W = [haar(N) for _ in range(L)]
        lhs = np.eye(N); rhs = np.eye(N); P = np.eye(N)
        for x in range(L):
            lhs = lhs @ (W[x] @ V[x]); rhs = rhs @ (P @ W[x] @ np.linalg.inv(P)); P = P @ V[x]
        ok(np.allclose(lhs, rhs @ P), "holonomy identity")
for beta in (400.0, 1600.0, 6400.0):
    c = ive(2, beta) / ive(1, beta)
    ok(abs(beta * (1 - c) / 2 - 0.75) < 1.0 / beta, "SU(2) chat")
th = np.linspace(-np.pi, np.pi, 1600, endpoint=False)
t1, t2 = np.meshgrid(th, th, indexing="ij"); t3 = -(t1 + t2)
z = [np.exp(1j * t1), np.exp(1j * t2), np.exp(1j * t3)]
vd = abs(z[0] - z[1]) ** 2 * abs(z[0] - z[2]) ** 2 * abs(z[1] - z[2]) ** 2
chi = z[0] + z[1] + z[2]; retr = np.cos(t1) + np.cos(t2) + np.cos(t3)
for beta in (800.0, 3200.0):
    w = vd * np.exp(beta / 3 * (retr - 3)); w /= w.sum()
    cF = float(np.sum(w * chi).real) / 3; cA = float(np.sum(w * (abs(chi) ** 2 - 1))) / 8
    ok(abs(beta * (1 - cF) / 3 - 4 / 3) < 40 / beta and abs(beta * (1 - cA) / 3 - 3) < 80 / beta, "SU(3) chat")
    ok(0 < cA < cF < 1, "SU(3) chat order")
for a in (0.02, 0.01, 0.005):
    beta = 4 / a ** 2
    E = -math.log(ive(2, beta) / ive(1, beta)) / a ** 2
    ok(abs(E - 0.375) < 0.4 * a, "SU(2) energy")

# (x) block averages of the free field: E Var(F|B) >= 2 sigma^4 sum alpha_B^2 and
# Var F <= 2 d^2 sum f^2, exactly by Gaussian traces
def free_block_check(shape, mu, fvals, blocks):
    d = len(shape); sites = list(product(*[range(n) for n in shape])); idx = {x: i for i, x in enumerate(sites)}
    n = len(sites)
    def mv(x, nu, s):
        y = list(x); y[nu] = (y[nu] + s) % shape[nu]; return tuple(y)
    P = np.zeros((n, n))
    for x in sites:
        P[idx[x], idx[x]] += 2 * d + mu ** 2
        for nu in range(d):
            for s_ in (1, -1):
                P[idx[x], idx[mv(x, nu, s_)]] -= 1
    C = np.linalg.inv(P)
    A = np.zeros((n, n))
    for x in sites:
        for nu in range(d):
            g = np.zeros(n); g[idx[mv(x, nu, 1)]] += 1; g[idx[x]] -= 1
            A += fvals[idx[x]] * np.outer(g, g)
    K = np.zeros((len(blocks), n))
    for b, z in enumerate(blocks):
        for e in product((0, 1), repeat=d):
            K[b, idx[tuple((2 * z[i] + e[i]) % shape[i] for i in range(d))]] = 2.0 ** -d
    Cc = C - C @ K.T @ np.linalg.solve(K @ C @ K.T, K @ C)
    evar = 2 * np.trace(A @ Cc @ A @ Cc) + 4 * np.trace(A @ Cc @ A @ (C - Cc))
    var = 2 * np.trace(A @ C @ A @ C)
    sig2 = 2 / (2 * d + mu ** 2); al = 0.0
    for z in blocks:
        xb = tuple(2 * zi for zi in z); yb = tuple(2 * z[i] + (1 if i < 2 else 0) for i in range(d))
        t = sum(d * fvals[idx[w]] + sum(fvals[idx[mv(w, nu, -1)]] for nu in range(d)) for w in (xb, yb))
        al += (t / 4) ** 2
    ok(evar >= 2 * sig2 ** 2 * al * (1 - 1e-9), "block-average lower bound")
    ok(var <= 2 * d * d * float(np.sum(fvals ** 2)) * (1 + 1e-9), "variance upper bound")
    ok(evar <= var * (1 + 1e-9), "conditioning reduces variance")
for shape in ((4, 4), (6, 6), (8, 4), (4, 4, 4)):
    allb = list(product(*[range(n // 2) for n in shape]))
    for mu in (0.1, 0.7, 2.0):
        for trial in range(3):
            fv = rng.random(int(np.prod(shape))) * (rng.random(int(np.prod(shape))) < 0.6)
            sub = [b for b in allb if rng.random() < 0.7] or allb[:1]
            free_block_check(shape, mu, fv, sub)

# (y) perturbed actions: the tilt inequalities behind the openness of D(n, alpha)
from scipy.optimize import linprog
def w1(mu, nu, Dm):
    k = len(mu)
    A_eq = []; b_eq = []
    for i in range(k):
        row = np.zeros(k * k); row[i * k:(i + 1) * k] = 1; A_eq.append(row); b_eq.append(mu[i])
    for j in range(k):
        row = np.zeros(k * k); row[j::k] = 1; A_eq.append(row); b_eq.append(nu[j])
    r = linprog(Dm.ravel(), A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=(0, None), method="highs")
    return r.fun
def tilt(mu, u):
    h = mu * np.exp(u - u.max()); return h / h.sum()
yrng = np.random.default_rng(29102026)
for trial in range(60):
    k = int(yrng.integers(4, 9))
    pts = yrng.normal(size=(k, 4)); pts /= np.linalg.norm(pts, axis=1)[:, None]
    Dm = np.arccos(np.clip(pts @ pts.T, -1, 1))          # geodesic metric on S^3 = SU(2)
    Delta = Dm.max()
    ell = float(yrng.uniform(0.05, 2.0))
    cen = yrng.integers(0, k, size=3); off = yrng.uniform(0, 1, size=3)
    u = ell * np.min(off[None, :] + Dm[:, cen], axis=1)    # ell-Lipschitz
    Om = u.max() - u.min()
    mu = yrng.dirichlet(np.ones(k)); nu = yrng.dirichlet(np.ones(k))
    lhs = w1(tilt(mu, u), tilt(nu, u), Dm); rhs = np.exp(Om) * (1 + 2 * Delta * ell) * w1(mu, nu, Dm)
    ok(lhs <= rhs + 1e-9, "tilt by a common weight")
    up = u + yrng.uniform(-1, 1, size=k) * float(yrng.uniform(0, 0.5))
    eta = np.abs(up - u).max()
    ok(w1(tilt(mu, u), tilt(mu, up), Dm) <= 2 * Delta * eta + 1e-9, "two tilts of one measure")
ok(48 * 2 == 96 and 2 * 48 * 48 == 96 * 48, "constants of the perturbation proposition")

print(f"all {checks} checks passed")
