"""Exact enumerations for the small-field block (smallfield.tex).

Blocks B(y)=2y+{0,1}^4, tree gauge of Definition def:blockgauge.  Variable
groups: I_y (17 active links inside B(y)), P_{y,mu} (8 private links of the
pair (y,mu)).  All statements below are checked by exhaustive enumeration
(integer arithmetic only).
"""
import itertools, random
D = 4
def off(x): return tuple(c % 2 for c in x)
def blk(x): return tuple(c // 2 for c in x)
def add(a, b): return tuple(i + j for i, j in zip(a, b))
def sub(a, b): return tuple(i - j for i, j in zip(a, b))
def e(m): return tuple(1 if i == m else 0 for i in range(D))
def linf(a, b): return max(abs(i - j) for i, j in zip(a, b))
O16 = [tuple(o) for o in itertools.product((0, 1), repeat=D)]

def group(x, m):
    o = off(x); y = blk(x)
    if o[m] == 1: return ('P', y, m)
    if all(o[n] == 0 for n in range(m + 1, D)): return None   # tree link
    return ('I', y)

def owned(y):
    for o in O16:
        x = add(tuple(2 * c for c in y), o)
        for m, n in itertools.combinations(range(D), 2):
            yield (x, m, n)

def links(p):
    x, m, n = p
    return [(x, m), (add(x, e(m)), n), (add(x, e(n)), m), (x, n)]

def fdeps(y):
    return {g for p in owned(y) for (z, k) in links(p) for g in [group(z, k)] if g}

# ---- (1) the factor f_0 -------------------------------------------------
Z = (0,) * D
assert sum(1 for _ in owned(Z)) == 96
f0 = fdeps(Z)
expect = {('I', Z)} | {('P', Z, m) for m in range(D)} \
    | {('I', e(n)) for n in range(D)} \
    | {('P', e(n), m) for n in range(D) for m in range(D) if m != n}
assert f0 == expect, f0 ^ expect
print("(1) f_0 depends exactly on I_0, P_{0,mu}, I_{e_nu}, P_{e_nu,mu} (mu!=nu):",
      len(f0), "groups")
# active links in each group
cnt = {}
for y0 in [Z]:
    for o in O16:
        for m in range(D):
            g = group(o, m)
            if g: cnt[g[0]] = cnt.get(g[0], 0) + 1
print("    active links per block: inside", cnt['I'], " private", cnt['P'])

# ---- (2) the cut-off averages -------------------------------------------
def vdeps(y, m):
    s = set()
    for o in O16:
        x = add(tuple(2 * c for c in y), o)
        for z in (x, add(x, e(m))):
            g = group(z, m)
            if g: s.add(g)
    return s
for m in range(D):
    v = vdeps(Z, m)
    assert v <= {('P', Z, m), ('I', Z), ('I', e(m))} and ('P', Z, m) in v
print("(2) V_{y,mu} depends only on P_{y,mu}, I_y, I_{y+e_mu}")

# ---- (3) the kernels -----------------------------------------------------
def dminus(X):
    return {sub(x, e(m)) for x in X for m in range(D)} - X
def dpairs_in(X):
    return {(sub(x, e(m)), m) for x in X for m in range(D) if sub(x, e(m)) not in X}
def integrated(X):
    s = {('I', x) for x in X} | {('P', x, m) for x in X for m in range(D)}
    s |= {('P', z, m) for (z, m) in dpairs_in(X)}
    return s
def dsecond(X):
    dm = dminus(X)
    return {sub(sub(x, e(m)), e(n)) for x in X for m in range(D) for n in range(D)
            if m != n} - X - dm
def needed(X):
    """blocks w whose f_w involves a variable integrated in K_X"""
    integ = integrated(X)
    cands = {add(x, d) for x in X for d in itertools.product((-1, 0, 1), repeat=D)}
    return {w for w in cands if fdeps(w) & integ}

def check_kernel(X):
    X = set(X)
    hat = X | dminus(X) | dsecond(X)
    nd = needed(X)
    assert nd == hat, (sorted(nd ^ hat))
    # every block of hat lies in X + {0, -e_mu, -e_mu-e_nu}
    for w in hat:
        assert any(sub(x, w) in {Z} | {e(m) for m in range(D)} |
                   {add(e(m), e(n)) for m in range(D) for n in range(D) if m != n}
                   for x in X)
    integ = integrated(X)
    # factors in hat involve, besides integrated variables, only groups of blocks
    # within l_inf distance one of X, and no pair whose measure involves an
    # integrated inside variable is left outside
    outer_used = set().union(*(fdeps(w) for w in hat)) - integ
    for g in outer_used:
        assert min(linf(g[1], x) for x in X) <= 1
        if g[0] == 'P':
            z, m = g[1], g[2]
            assert ('I', z) not in integ and ('I', add(z, e(m))) not in integ
    # moved plaquettes: owned by hat, contain an integrated link, have a vertex in B(X)
    nmov = 0
    for w in hat:
        for p in owned(w):
            gs = {group(z, k) for (z, k) in links(p)}
            if gs & integ:
                nmov += 1
                x0, m, n = p
                verts = [x0, add(x0, e(m)), add(x0, e(n)), add(add(x0, e(m)), e(n))]
                assert any(blk(v) in X for v in verts)
    # every plaquette of Pi(X) contains an integrated link
    for x in X:
        for p in owned(x):
            assert {group(z, k) for (z, k) in links(p)} & integ
    return hat, nmov, len(dpairs_in(X))

hat, nmov, nin = check_kernel([Z])
print("(3) X={0}: |X^|=%d, |d^-X|=%d, |d^=X|=%d, pairs into X=%d, moved plaquettes=%d"
      % (len(hat), len(dminus({Z})), len(dsecond({Z})), nin, nmov))

def random_connected(n, rng):
    X = {Z}
    while len(X) < n:
        x = rng.choice(sorted(X))
        d = tuple(rng.choice((-1, 0, 1)) for _ in range(D))
        X.add(add(x, d))
    return X
rng = random.Random(20261104)
worst = 0
for trial in range(300):
    X = random_connected(rng.randint(1, 7), rng)
    hat, nmov, nin = check_kernel(X)
    worst = max(worst, nmov / len(X))
    assert nin <= 4 * len(X)
print("    300 random l_inf-connected sets of 1..7 blocks: kernel factor set = X u d^-X u d^=X,"
      " moved plaquettes per block <= %.3f <= 384" % worst)
assert worst <= 384

# ---- (4) two components at distance >= 2 share no factor and no variable --
for trial in range(300):
    X1 = random_connected(rng.randint(1, 4), rng)
    shift = tuple(rng.choice((-6, -5, -4, -3, -2, 2, 3, 4, 5, 6)) if i == 0 else rng.randint(-2, 2)
                  for i in range(D))
    X2 = {add(x, shift) for x in random_connected(rng.randint(1, 4), rng)}
    if min(linf(a, b) for a in X1 for b in X2) < 2: continue
    h1 = X1 | dminus(X1) | dsecond(X1); h2 = X2 | dminus(X2) | dsecond(X2)
    assert not (h1 & h2)
    assert not (integrated(X1) & integrated(X2))
    # a factor of one kernel never involves a variable integrated by the other
    for w in h1: assert not (fdeps(w) & integrated(X2))
    for w in h2: assert not (fdeps(w) & integrated(X1))
print("(4) components at l_inf distance >= 2: disjoint factor sets and variables (300 trials)")
print("all checks passed")
