"""Checks for the loop remainders of the one-step bound.

(A) s(g_1...g_A) <= A sum s(g_j) on random unitary matrices.
(B) The ordered filling by adjacent transpositions: on random closed lattice
    paths in Z^4 with random SU(3) links, V_C equals the ordered product of the
    transported plaquettes, the number of plaquettes is the number of
    inversions <= l(l-1)/2, and every plaquette lies in the bounding box.
(C) At a flat background, the second derivative of s(V_C(t)) along V_0 e^{tw}
    is at most (A_C/N) sum_p m_Cp ||D_p w||^2 (finite differences).
(D) R_m = eps (1 - cos t) cos(m t): |R_m| <= eps s and R_m'' = -eps m^2 at pi/2.
(E) The constants of the corollary for SU(2), lambda = 1e100.
"""
import itertools, math, random
import numpy as np
from mpmath import mp, mpf, sqrt, sin, pi

rng = np.random.default_rng(1)
random.seed(1)
fails = 0
def check(ok, msg):
    global fails
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails += 1

def haar_unitary(n, special=False):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / math.sqrt(2)
    q, r = np.linalg.qr(z)
    q = q @ np.diag(np.diag(r) / abs(np.diag(r)))
    if special:
        q = q / np.linalg.det(q) ** (1.0 / n)
    return q

def s(g):
    return 1 - np.trace(g).real / g.shape[0]

# (A)
worst = 0.0
for _ in range(2000):
    A = random.randint(1, 8); N = random.randint(1, 4)
    r = 10.0 ** random.uniform(-4, 0.5)
    gs = []
    for _ in range(A):
        z = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
        X = r * (z - z.conj().T) / 2
        w, U = np.linalg.eigh(-1j * X)
        gs.append(U @ np.diag(np.exp(1j * w)) @ U.conj().T)
    lhs = s(np.linalg.multi_dot(gs + [np.eye(N)]))
    rhs = A * sum(s(g) for g in gs)
    worst = max(worst, (lhs - rhs) / max(rhs, 1e-300))
check(worst <= 1e-9, f"(A) product inequality on 2000 random products near and far from I, max (lhs-rhs)/rhs = {worst:.2e}")

# (B) lattice paths
D = 4
def steps_of_random_closed_path(l):
    while True:
        st = [(random.randrange(D), random.choice((1, -1))) for _ in range(l - 1)]
        pos = [0] * D
        for mu, e in st:
            pos[mu] += e
        nz = [mu for mu in range(D) if pos[mu] != 0]
        if sum(abs(p) for p in pos) == 1:
            mu = nz[0]
            st.append((mu, -pos[mu]))
            return st

class Field:
    def __init__(self, N=3):
        self.N = N; self.links = {}
    def link(self, x, mu):
        key = (tuple(x), mu)
        if key not in self.links:
            self.links[key] = haar_unitary(self.N, special=True)
        return self.links[key]
    def step(self, x, mu, e):
        if e == 1:
            return self.link(x, mu)
        y = list(x); y[mu] -= 1
        return self.link(y, mu).conj().T
    def hol(self, x, st):
        x = list(x); g = np.eye(self.N, dtype=complex)
        for mu, e in st:
            g = g @ self.step(x, mu, e); x[mu] += e
        return g, x

def fill(x0, st):
    """Bubble sort by direction; return list of (h-path prefix, x, plaquette word)."""
    st = list(st); factors = []
    changed = True
    while changed:
        changed = False
        for i in range(len(st) - 1):
            if st[i][0] > st[i + 1][0]:
                factors.append((list(st[:i]), st[i], st[i + 1]))
                st[i], st[i + 1] = st[i + 1], st[i]
                changed = True
    return factors, st

def inversions(st):
    return sum(1 for i in range(len(st)) for j in range(i + 1, len(st)) if st[i][0] > st[j][0])

maxerr = 0.0; okcount = True; okbox = True; npaths = 0
for l in range(4, 15, 2):
    for _ in range(40):
        st = steps_of_random_closed_path(l)
        F = Field(3)
        x0 = [0] * D
        VC, end = F.hol(x0, st)
        assert end == x0
        factors, sorted_st = fill(x0, st)
        okcount &= len(factors) == inversions(st) <= l * (l - 1) // 2
        # bounding box of the path
        pts = [list(x0)]; x = list(x0)
        for mu, e in st:
            x[mu] += e; pts.append(list(x))
        lo = [min(p[m] for p in pts) for m in range(D)]; hi = [max(p[m] for p in pts) for m in range(D)]
        prod = np.eye(3, dtype=complex)
        for pre, u, v in factors:
            h, xs = F.hol(x0, pre)
            # loop u v u'^{-1} v'^{-1} from xs
            P, _ = F.hol(xs, [u, v, (u[0], -u[1]), (v[0], -v[1])])
            # plaquette corners
            corners = [xs, [xs[m] + (u[1] if m == u[0] else 0) for m in range(D)],
                       [xs[m] + (v[1] if m == v[0] else 0) for m in range(D)],
                       [xs[m] + (u[1] if m == u[0] else 0) + (v[1] if m == v[0] else 0) for m in range(D)]]
            okbox &= all(lo[m] <= c[m] <= hi[m] for c in corners for m in range(D))
            prod = prod @ (h @ P @ h.conj().T)
        Vsorted, _ = F.hol(x0, sorted_st)
        maxerr = max(maxerr, np.abs(Vsorted - np.eye(3)).max(), np.abs(prod - VC).max())
        npaths += 1
check(maxerr < 1e-10, f"(B) ordered filling identity on {npaths} random closed paths, lengths 4-14, max error {maxerr:.1e}")
check(okcount, "(B) number of plaquettes = number of inversions <= l(l-1)/2")
check(okbox, "(B) every plaquette of the filling lies in the bounding box of the path")

# (C) flat Hessian
def su_alg(N):
    z = rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    X = (z - z.conj().T) / 2
    return X - np.trace(X) / N * np.eye(N)

def expm(X):
    w, U = np.linalg.eigh(-1j * X)
    return U @ np.diag(np.exp(1j * w)) @ U.conj().T

worst = -1e9
for trial in range(30):
    N = 2 + trial % 2
    l = random.choice((4, 6, 8, 10))
    st = steps_of_random_closed_path(l)
    # flat background: pure gauge V_{x,mu} = g_x g_{x+mu}^{-1}
    gs = {}
    def gsite(x):
        k = tuple(x)
        if k not in gs: gs[k] = haar_unitary(N, special=True)
        return gs[k]
    x0 = [0] * D
    factors, _ = fill(x0, st)
    # collect links of the path and of the plaquettes
    def path_links(x, word):
        x = list(x); out = []
        for mu, e in word:
            y = list(x)
            if e == -1: y[mu] -= 1
            out.append((tuple(y), mu)); x[mu] += e
        return out
    W = {}
    allk = set(path_links(x0, st))
    plist = []
    for pre, u, v in factors:
        x = list(x0)
        for mu, e in pre: x[mu] += e
        word = [u, v, (u[0], -u[1]), (v[0], -v[1])]
        plist.append((x, word)); allk |= set(path_links(x, word))
    for k in allk:
        W[k] = su_alg(N)
    def field_at(t):
        F = Field(N)
        for (y, mu) in allk:
            yp = list(y); yp[mu] += 1
            F.links[(y, mu)] = gsite(y) @ gsite(yp).conj().T @ expm(t * W[(y, mu)])
        return F
    h = 1e-4
    def sC(t):
        F = field_at(t); g = F.hol(x0, st)[0]
        return np.sum(np.abs(np.eye(N) - g) ** 2) / (2 * N)   # = s(g) for unitary g, without cancellation
    second = (sC(h) - 2 * sC(0) + sC(-h)) / h ** 2
    # D_p w via derivative of plaquette holonomy
    bound = 0.0
    for x, word in plist:
        Fp, Fm = field_at(h), field_at(-h)
        dP = (Fp.hol(x, word)[0] - Fm.hol(x, word)[0]) / (2 * h)
        bound += np.sum(np.abs(dP) ** 2)
    bound *= len(plist) / N
    worst = max(worst, (second - bound) / (bound + 1e-3))
    assert abs(sC(0)) < 1e-12
check(worst <= 1e-6, f"(C) flat pure-gauge backgrounds: d2 s(V_C) <= (A_C/N) sum ||D_p w||^2, max excess / (bound + 1e-3) {worst:.1e}")

# (D)
mp.dps = 40
ok = True
eps = mpf(1)
for m in (4, 8, 12, 40):
    f = lambda t: eps * (1 - mp.cos(t)) * mp.cos(m * t)
    d2 = mp.diff(f, mp.pi / 2, 2)
    ok &= abs(d2 + eps * m * m) < mpf(10) ** -20
    ok &= all(abs(f(t)) <= eps * (1 - mp.cos(t)) + mpf(10) ** -30 for t in [mpf(k) / 37 for k in range(0, 233)])
check(ok, "(D) R_m'' = -eps m^2 at pi/2 for m = 4, 8, 12, 40, and |R_m| <= eps s")

# (E)
cstar = 1 / (1 + sqrt(mpf(21) / 8)) ** 2
lam = mpf(10) ** 100; dG = mpf(1) / 3; e = mpf('1.63e-9')
def epsdd(n):
    Cn = 1 + 1 / (4 * sin(pi / (2 * n)) ** 2)
    return (e + dG * Cn / lam) / (2 * cstar)
print("2 c_* =", mp.nstr(2 * cstar, 8))
for n in (2, 1024, 10 ** 6, 10 ** 40):
    print(f"   n = {n}: eps'' = {mp.nstr(epsdd(n), 6)}")
check(epsdd(10 ** 40) <= mpf('5.60e-9'), "(E) eps'' <= 5.60e-9 for n <= 1e40")
lo_, hi_ = mpf(10) ** 40, mpf(10) ** 60
for _ in range(200):
    mid = sqrt(lo_ * hi_)
    if epsdd(mid) <= mpf('0.1'): lo_ = mid
    else: hi_ = mid
print("   largest n with eps'' <= 0.1:", mp.nstr(lo_, 6))
check(epsdd(mpf('9.2e49')) <= mpf('0.1'), "(E) eps'' <= 0.1 for n <= 9.2e49")
print(f"{fails} failures")
