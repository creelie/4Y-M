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
  (e) the exact arithmetic of the proofs.
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

print(f"all {checks} checks passed")
