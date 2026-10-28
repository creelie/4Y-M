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
      the diameter squared is 9 against the Haar mean 6 of the squared distance.
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

print(f"all {checks} checks passed")
