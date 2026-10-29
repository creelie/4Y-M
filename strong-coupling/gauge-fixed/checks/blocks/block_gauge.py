"""Block-rooted axial gauge with full-block conditioning: the counts behind the
uniform gap and the localisation of the block-constrained cube responses.

Blocks are the 2-cubes 2y + {0,1}^4. A link (x, mu) with offset o = x mod 2 is a
tree link when o_mu = 0 and o_nu = 0 for every nu > mu. For a coarse link (y, mu)
the full-block average V_{y,mu} is the mean, over the 16 straight paths of length
two in direction mu starting in the block, of the sum of the field along the path.

Parts:
 (A) tree: 15 links per block, a spanning tree of the block; 17 active intra-block
     links, 32 inter-block links, 49 classes.
 (B) chain: every active intra-block link is a signed sum of k <= 3 plaquette
     curvatures of its block (all other links of the plaquettes being tree links
     or earlier links of the chain), and every plaquette carries load <= 8.
 (C) faces: the 12 face plaquettes of a coarse link are the 12 edges of the
     transverse 3-cube; V = (1/8) sum(face links) + (S_y + S_{y+e_mu})/16; every
     intra link gets weight 9/2 + 1 = 11/2 when the face bounds are summed.
 (D) constants of the saddle-point system: lambda_min(B B^T) >= 1/16, row sums
     of |H| at most 24, of the scaled system at most 192 (s = 96), and the
     resulting localisation constants.
 (E) exact linear algebra in floating point (comparison only): the bulk gap of
     the constrained symbol and the gap and responses on the cubes n = 2, 4.
 (F) the invisible line class of the line-field gauge (face averages, n = 4).
Run: python3 block_gauge.py. Needs numpy and scipy."""
import itertools, math, sys
from fractions import Fraction as Fr
import numpy as np

E = np.eye(4, dtype=int)
OFFS = list(itertools.product(range(2), repeat=4))
ok = []

def check(name, cond):
    ok.append(bool(cond)); print(('ok  ' if cond else 'FAIL'), name)

def tree(o, mu):
    o = [c % 2 for c in o]
    return o[mu] == 0 and all(o[nu] == 0 for nu in range(mu + 1, 4))

# (A)
tl = [(o, mu) for o in OFFS for mu in range(4) if tree(o, mu)]
parent = {o: o for o in OFFS}
def find(a):
    while parent[a] != a: a = parent[a]
    return a
cyc = False
for (o, mu) in tl:
    a, b = find(o), find(tuple(np.array(o) + E[mu]))
    if a == b: cyc = True
    parent[a] = b
check('(A) 15 tree links forming a spanning tree of the block',
      len(tl) == 15 and not cyc and len({find(o) for o in OFFS}) == 1)
intra = [(o, mu) for o in OFFS for mu in range(4) if o[mu] == 0 and not tree(o, mu)]
inter = [(o, mu) for o in OFFS for mu in range(4) if o[mu] == 1]
check('(A) 17 active intra-block links, 32 inter-block links, 49 classes',
      len(intra) == 17 and len(inter) == 32 and len(intra) + len(inter) == 49)

# (B)
load, kmax = {}, 0
for (z, mu) in intra:
    ones = [r for r in range(mu + 1, 4) if z[r] == 1]
    kmax = max(kmax, len(ones)); x = list(z)
    for nu in sorted(ones, reverse=True):
        x[nu] = 0
        xa = np.array(x)
        # plaquette at x in the plane (mu, nu): (x,mu) earlier, (x,nu), (x+e_mu,nu) tree
        assert tree(tuple(x), nu) and tree(tuple(xa + E[mu]), nu)
        key = (tuple(x), mu, nu); load[key] = load.get(key, 0) + len(ones)
    assert tree(tuple(x), mu)
check('(B) chains of length <= 3 ending in tree links', kmax == 3)
check('(B) plaquette load at most 8 (so I <= 8 Q_in)', max(load.values()) == 8)

# (C)
good = True
for mu in range(4):
    face = [o for o in OFFS if o[mu] == 1]
    pl = [(o, nu) for o in face for nu in range(4) if nu != mu and o[nu] == 0]
    edges = {frozenset([o, tuple(np.array(o) + E[nu])]) for (o, nu) in pl}
    good &= len(pl) == 12 and len(edges) == 12
    coef = {}
    for o in OFFS:
        for k in range(2):
            y = np.array(o) + k * E[mu]
            key = ('F', tuple(y)) if y[mu] == 1 else ('I', tuple(y % 2), int(y[mu] // 2))
            coef[key] = coef.get(key, 0) + Fr(1, 16)
    good &= all(v == (Fr(1, 8) if k[0] == 'F' else Fr(1, 16)) for k, v in coef.items())
    good &= sum(1 for k in coef if k[0] == 'F') == 8 and sum(1 for k in coef if k[0] == 'I') == 16
check('(C) 12 face plaquettes = 12 edges of the 3-cube; V = (1/8)sum F + (S+S\')/16', good)
w = {}
for mu in range(4):
    for o in OFFS:
        for nu in range(4):
            if nu == mu or o[nu] != 0: continue
            if o[mu] == 1:   # (x,nu) in block y, used by face (y,mu)
                w[(o, nu)] = w.get((o, nu), 0) + Fr(3, 2)
            else:            # (x+e_mu,nu) in block y+e_mu, used by face (y,mu)
                w[(o, nu)] = w.get((o, nu), 0) + Fr(3, 2)
    for o in OFFS:
        if o[mu] == 0:
            w[(o, mu)] = w.get((o, mu), 0) + 2 * Fr(1, 2)
check('(C) every intra link has total weight 11/2', set(w.values()) == {Fr(11, 2)})
a = Fr(1, 1) / (Fr(13, 2) * 8)
check('(C) gap constant 1/52', a == Fr(1, 52))

# (D)
s, h, b = 96, 24, Fr(1, 16)
x1 = 1 / (s * math.sqrt(float(b)))   # ||x_1|| per unit ||g'||
xn = math.hypot(52, 52 * h * x1 + x1)  # ||x|| for ||f||^2 + ||g'||^2 = 1
yn = (1 + h * xn) / (s * math.sqrt(float(b)))
M0 = math.hypot(xn, yn)
S = max(h + s * Fr(1, 8), 2 * s)
R = 1  # midpoints of links sharing a plaquette, and of a row and its links, are within distance 1
lam = math.log(1 + 1 / (2 * float(S) * M0)) / R
q = math.exp(-2 * lam)
sum3 = q * (1 + 4 * q + q * q) / (1 - q) ** 4; sum1 = q / (1 - q) ** 2
W = 4 * math.exp(2 * lam) * (1 + 64 * sum3 + 16 * sum1)
rhs = 6 * math.sqrt(3) + s / 16
M = math.sqrt(W) * 2 * M0 * math.exp(R * lam) * rhs
print(f'    M0 <= {M0:.2f}, S = {S}, lambda = {lam:.4e}, sqrt(W) = {math.sqrt(W):.3e}, M <= {M:.3e}')
check('(D) M0 <= 105, S = 192, lambda >= 2.4e-5, ||u_j||_1 <= 10^14',
      M0 <= 105 and S == 192 and lam >= 2.4e-5 and M <= 1e14)

mid = lambda x, mu: np.array(x, float) + 0.5 * E[mu]
dmax = 0
for mu in range(4):
    c = 0.5 * np.ones(4) + E[mu]
    for o in OFFS:
        for k in range(2):
            y = np.array(o) + k * E[mu]
            dmax = max(dmax, np.abs(mid(y, mu) - c).max())
for o in OFFS:
    x = np.array(o)
    for mu in range(4):
        for nu in range(mu + 1, 4):
            ms = [mid(x, mu), mid(x + E[mu], nu), mid(x + E[nu], mu), mid(x, nu)]
            dmax = max(dmax, max(np.abs(p - q).max() for p in ms for q in ms))
check('(D) R = 1: links of a plaquette, and a row and its links, within distance 1', dmax <= 1)

# (E) floating-point comparison
import scipy.sparse as sp
def block_rows(L=2):
    links = [(o, mu) for o in OFFS for mu in range(4)]
    act = [l for l in links if not tree(*l)]; idx = {l: i for i, l in enumerate(act)}
    def red(x, mu):
        x = np.array(x); return (tuple(x % 2), mu), tuple(x // 2)
    P = []
    for o in OFFS:
        x = np.array(o)
        for mu in range(4):
            for nu in range(mu + 1, 4):
                r = []
                for (y, m, sg) in [(x, mu, 1), (x + E[mu], nu, 1), (x + E[nu], mu, -1), (x, nu, -1)]:
                    l, c = red(y, m)
                    if l in idx: r.append((idx[l], c, sg))
                P.append(r)
    C = []
    for mu in range(4):
        r = {}
        for o in OFFS:
            for k in range(2):
                l, c = red(np.array(o) + k * E[mu], mu)
                if l in idx: r[(idx[l], c)] = r.get((idx[l], c), 0) + 1 / 16
        C.append([(i, c, v) for (i, c), v in r.items()])
    return P, C, len(act)
P, C, nA = block_rows()
def sym(rows, q):
    Mx = np.zeros((len(rows), nA), complex)
    for i, r in enumerate(rows):
        for (j, c, v) in r: Mx[i, j] += v * np.exp(1j * np.dot(q, c))
    return Mx
best = 9
for qq in itertools.product(range(4), repeat=4):
    q = 2 * np.pi * np.array(qq) / 4
    D = sym(P, q); H = D.conj().T @ D; B = sym(C, q)
    u, sv, vh = np.linalg.svd(B); N = vh[(sv > 1e-10).sum():].conj().T
    best = min(best, np.linalg.eigvalsh(N.conj().T @ H @ N)[0])
print(f'    bulk min lambda_min(H_c) on the 4^4 grid = {best:.12f}')
check('(E) bulk constrained gap >= 1/52', best >= 1 / 52)

def cube(n):
    inV = lambda x: all(0 <= c < n for c in x)
    V = [(x, mu) for x in itertools.product(range(n), repeat=4) for mu in range(4) if not tree(x, mu)]
    vi = {l: i for i, l in enumerate(V)}; ext = {}
    def var(x, mu):
        x = tuple(int(c) for c in x)
        if tree(x, mu): return None
        if (x, mu) in vi: return ('V', vi[(x, mu)])
        ext.setdefault((x, mu), len(ext)); return ('E', ext[(x, mu)])
    rows = []
    for x in itertools.product(range(-1, n), repeat=4):
        xa = np.array(x)
        for mu in range(4):
            for nu in range(mu + 1, 4):
                ls = [(xa, mu, 1), (xa + E[mu], nu, 1), (xa + E[nu], mu, -1), (xa, nu, -1)]
                if any(inV(tuple(y)) and not tree(tuple(y), m) for (y, m, sg) in ls):
                    rows.append([(var(y, m), sg) for (y, m, sg) in ls])
    cons = []
    for y in itertools.product(range(n // 2), repeat=4):
        for mu in range(4):
            r = {}
            for o in OFFS:
                for k in range(2):
                    v = var(2 * np.array(y) + np.array(o) + k * E[mu], mu)
                    if v: r[v] = r.get(v, 0) + 1 / 16
            cons.append(list(r.items()))
    def mat(rs, ncol, kind):
        I, J, Dd = [], [], []
        for i, r in enumerate(rs):
            for (v, sg) in r:
                if v and v[0] == kind: I.append(i); J.append(v[1]); Dd.append(sg)
        return sp.csr_matrix((Dd, (I, J)), shape=(len(rs), ncol))
    DV, DE = mat(rows, len(V), 'V'), mat(rows, len(ext), 'E')
    BV, BE = mat(cons, len(V), 'V'), mat(cons, len(ext), 'E')
    H = (DV.T @ DV).toarray(); B = BV.toarray()
    u, sv, vh = np.linalg.svd(B); N = vh[(sv > 1e-12).sum():].T
    gap = np.linalg.eigvalsh(N.T @ H @ N)[0]
    K = np.block([[H, B.T], [B, np.zeros((B.shape[0],) * 2)]])
    X = np.linalg.solve(K, np.vstack([-(DV.T @ DE).toarray(), -BE.toarray()]))
    return gap, np.abs(X[:len(V)]).sum(axis=0).max(), np.linalg.eigvalsh(B @ B.T)[0], np.linalg.eigvalsh(H)[-1]
for n in (2, 4):
    g, Cm, bb, hh = cube(n)
    print(f'    n = {n}: gap {g:.6f}, max ||u_j||_1 {Cm:.4f}, lambda_min(BB^T) {bb:.4f}, ||H|| {hh:.4f}')
    check(f'(E) cube n = {n}: gap >= 1/52, BB^T >= 1/16, ||H|| <= 24', g >= 1 / 52 and bb >= 1 / 16 and hh <= 24)

# (F) invisible class of the line field (face averages, L = 2, corners even)
SIGMA = {(0,0,0,0):0,(1,0,0,0):0,(0,1,0,0):2,(1,1,0,0):3,(0,0,1,0):3,(1,0,1,0):1,(0,1,1,0):2,(1,1,1,0):1,
         (0,0,0,1):1,(1,0,0,1):2,(0,1,0,1):1,(1,1,0,1):3,(0,0,1,1):3,(1,0,1,1):2,(0,1,1,1):0,(1,1,1,1):0}
sig = lambda x: SIGMA[tuple(int(c) % 2 for c in x)]
n = 4
def h(x):  # indicator of the line through (*,1,1,1)+..., direction 0, transverse coords (1,1,1)
    return 1 if tuple(x[1:]) == (1, 1, 1) else 0
ok_inv, cls = True, set()
for y in itertools.product(range(-1, n // 2 + 1), repeat=4):
    for mu in range(4):
        tot = 0
        for o in OFFS:
            if o[mu] != 0: continue
            z = 2 * np.array(y) + np.array(o)
            tot += h(z) - h(z + 2 * E[mu])
        ok_inv &= tot == 0
for x in itertools.product(range(2), repeat=4):
    for mu in range(4):
        if sig(x) == mu: continue
        xa = np.array(x)
        if abs(h(xa) - h(xa + E[mu])) == 1: cls.add((x, mu))
check('(F) the invisible line changes no face average', ok_inv)
check('(F) exactly twelve classes meet an invisible line once', len(cls) == 12)

print(f'{sum(ok)} checks passed, {len(ok) - sum(ok)} failed')
sys.exit(0 if all(ok) else 1)
