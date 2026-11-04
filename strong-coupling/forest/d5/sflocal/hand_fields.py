"""Checks of the hand-constructed fields phi_1 (thm:sfshelllin) and phi_2
(thm:sfperp) in exact rational arithmetic.  The proofs in the text derive
every value checked here by hand; this script only repeats them."""
import itertools
from fractions import Fraction as Fr

n = 16                       # fine side of the test torus (coarse side 8)
y = (3, 3, 3, 3)

def add(z, *vs):
    z = list(z)
    for v in vs:
        z[v] = (z[v] + 1) % n
    return tuple(z)

def block(z):
    return tuple(c // 2 for c in z)

def is_tree(z, mu):
    o = [c % 2 for c in z]
    return o[mu] == 0 and all(o[nu] == 0 for nu in range(mu + 1, 4))

def compensate(main):
    """Put -v/7 on the other seven private links of the pair of each main link."""
    phi = dict(main)
    pairs = {}
    for (z, mu), v in main.items():
        assert z[mu] % 2 == 1, "main link not private"
        k = (block(z), mu)
        assert k not in pairs, "two main links in one pair"
        pairs[k] = v
    for (b, mu), v in pairs.items():
        for x in itertools.product(*[(2 * c, 2 * c + 1) for c in b]):
            if x[mu] % 2 == 1 and (x, mu) not in main:
                phi[(x, mu)] = phi.get((x, mu), 0) - v / 7
    return phi, pairs

def V(phi, b, mu):
    s = Fr(0)
    for x in itertools.product(*[(2 * c, 2 * c + 1) for c in b]):
        if x[mu] % 2 == 1:
            s += Fr(1, 8) * phi.get((x, mu), 0)
        else:
            s += Fr(1, 16) * phi.get((x, mu), 0)
    bb = list(b); bb[mu] = (bb[mu] + 1) % (n // 2); bb = tuple(bb)
    for x in itertools.product(*[(2 * c, 2 * c + 1) for c in bb]):
        if x[mu] % 2 == 0:
            s += Fr(1, 16) * phi.get((x, mu), 0)
    return s

def F(phi, z, mu, nu):
    return (phi.get((z, mu), 0) + phi.get((add(z, mu), nu), 0)
            - phi.get((add(z, nu), mu), 0) - phi.get((z, nu), 0))

def plaquettes_touching(phi):
    pl = set()
    for (z, mu) in phi:
        for nu in range(4):
            if nu == mu:
                continue
            a, b = min(mu, nu), max(mu, nu)
            pl.add((z, a, b))
            w = list(z); w[nu] = (w[nu] - 1) % n
            pl.add((tuple(w), a, b))
    return pl

def links_of(z, mu, nu):
    return [(z, mu), (add(z, mu), nu), (add(z, nu), mu), (z, nu)]

def common_checks(phi):
    assert all(not is_tree(z, mu) for (z, mu) in phi), "tree link used"
    blocks = {block(z) for (z, mu) in phi}
    bad = [(b, mu) for b0 in blocks for mu in range(4)
           for b in [b0, tuple((b0[i] - (i == mu)) % (n // 2) for i in range(4))]
           if V(phi, b, mu) != 0]
    assert not bad, bad
    return max(abs(v) for v in phi.values())

def VX_links():
    S = set()
    for x in itertools.product(*[(2 * c, 2 * c + 1) for c in y]):
        for mu in range(4):
            if x[mu] % 2 == 0 and not is_tree(x, mu):
                S.add((x, mu))
    for mu in range(4):
        for b in [y, tuple((y[i] - (i == mu)) for i in range(4))]:
            for x in itertools.product(*[(2 * c, 2 * c + 1) for c in b]):
                if x[mu] % 2 == 1:
                    S.add((x, mu))
    return S

# ---------------------------------------------------------------- phi_1
print("== phi_1 (shell) ==")
x = tuple(2 * c + o for c, o in zip(y, (1, 1, 1, 0)))
h = Fr(1, 2)
main1 = {(add(x, 1), 0): -h, (add(x, 2), 0): h, (add(x, 0), 1): h,
         (add(x, 2), 1): -h, (add(x, 0), 2): -h, (add(x, 1), 2): h}
phi1, pairs1 = compensate(main1)
print("pairs", sorted(pairs1))
print("max |phi_1|", common_checks(phi1), "  V = 0 on every coarse link, no tree link")
faces = [(x, 0, 1), (x, 0, 2), (x, 1, 2)]
outer = [(add(x, 2), 0, 1), (add(x, 1), 0, 2), (add(x, 0), 1, 2)]
print("F1,F2,F3", [F(phi1, *p) for p in faces], " outer", [F(phi1, *p) for p in outer])
others = [p for p in plaquettes_touching(phi1) if p not in faces + outer]
print("max |F| off the six faces", max(abs(F(phi1, *p)) for p in others), "(5/7 =", Fr(5, 7), ")")
VX = VX_links()
print("phi_1 vanishes on V_X:", all(l not in VX for l in phi1))
for p in outer:
    assert all(l not in VX for l in links_of(*p)), "outer face meets V_X"
print("outer faces contain no link of V_X; owners", [block(p[0]) for p in outer])
# direction w: (x,1) -1, path link (2y+e1,1) +1 ; second: (x,2) -1, (2y+e2,2) +1
for mu in (1, 2):
    path = (tuple(2 * c + (i == mu) for i, c in enumerate(y)), mu)
    d = {(x, mu): Fr(-1), path: Fr(1)}
    assert (x, mu) in VX and path in VX
    dF = [F(d, *p) for p in faces]
    mx = max(abs(F(d, *p)) for p in plaquettes_touching(d))
    print(" direction", mu, ": change of (F1,F2,F3)", dF, " largest change", mx,
          " V unchanged", V(d, y, mu) == 0)

# ---------------------------------------------------------------- phi_2
print("== phi_2 (Z_perp) ==")
xp = tuple(2 * c + 1 for c in y)
s, t = Fr(11, 20), Fr(1, 4)
main2 = {}
def put(z, mu, v):
    main2[(z, mu)] = main2.get((z, mu), 0) + v
put(add(xp, 2), 0, s); put(add(xp, 2), 1, -s); put(add(xp, 0), 2, -s); put(add(xp, 1), 2, s)
put(add(xp, 3), 0, -s); put(add(xp, 3), 1, s); put(add(xp, 0), 3, s); put(add(xp, 1), 3, -s)
put(add(xp, 0, 2), 1, -2 * s / 3); put(add(xp, 0, 1), 2, -s / 3)
put(add(xp, 0, 3), 1, 2 * s / 3); put(add(xp, 0, 1), 3, s / 3)
put(add(xp, 1, 3), 2, t); put(add(xp, 1, 2), 3, -t); put(add(xp, 0, 2), 3, t); put(add(xp, 0, 3), 2, -t)
phi2, pairs2 = compensate(main2)
print("16 pairs distinct:", len(pairs2) == 16)
print("max |phi_2|", common_checks(phi2), "  V = 0 on every coarse link, no tree link")
print("phi_2 vanishes on V_X:", all(l not in VX for l in phi2))
oA = [(add(xp, 2), 0, 1), (add(xp, 1), 0, 2), (add(xp, 0), 1, 2)]
oB = [(add(xp, 3), 0, 1), (add(xp, 1), 0, 3), (add(xp, 0), 1, 3)]
for p in oA + oB:
    assert all(l not in VX for l in links_of(*p))
vA = [F(phi2, *p) for p in oA]; vB = [F(phi2, *p) for p in oB]
RA = vA[0] - vA[1] + vA[2]; RB = vB[0] - vB[1] + vB[2]
print("outer A", vA, "outer B", vB, "R_A", RA, "R_B", RB, "R_A-2", RA - 2, "2+R_B", 2 + RB)
nony = [p for p in plaquettes_touching(phi2) if block(p[0]) != y]
print("max |F| on plaquettes not owned by y", max(abs(F(phi2, *p)) for p in nony))
names = {}
for (z, mu) in main2:
    d = tuple((z[i] - xp[i]) % n for i in range(4))
    names[(z, mu)] = "x'" + "".join("+e%d" % i for i in range(4) if d[i]) + ",%d" % mu
print("plaquettes not owned by y with two or more main links:")
cnt = 0
for p in sorted(nony):
    ls = links_of(*p)
    m = [(i, l) for i, l in enumerate(ls) if l in main2]
    if len(m) < 2:
        continue
    cnt += 1
    d = tuple((p[0][i] - xp[i]) % n for i in range(4))
    nm = "x'" + "".join("%+de%d" % ((-1 if d[i] == n - 1 else 1), i) for i in range(4) if d[i])
    sg = [1, 1, -1, -1]
    comp = [(i, phi2[l]) for i, l in enumerate(ls) if l in phi2 and l not in main2]
    print("  (%s;%d%d) main" % (nm, p[1], p[2]),
          [(names[l], sg[i] * main2[l]) for i, l in m],
          "comp", [(sg[i] * v) for i, v in comp], "F =", F(phi2, *p), "owner", block(p[0]))
print("count", cnt)
one = [p for p in nony if sum(l in main2 for l in links_of(*p)) <= 1]
print("max |F| with at most one main link", max(abs(F(phi2, *p)) for p in one), "<= 10s/7 =", 10 * s / 7)
print("site bound 3 max|phi_2| =", 3 * max(abs(v) for v in phi2.values()))
