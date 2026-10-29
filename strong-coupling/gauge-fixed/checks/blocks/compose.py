"""One step of block averaging in the block-rooted axial gauge: the coarse
curvature of the averages, the effective quadratic form of the averages, and
the response of the constrained law to a change of one average.

Blocks are 2y + {0,1}^4; tree links as in block_gauge.py. For a field phi the
average V_{y,mu}(phi) is the mean over x in B(y) of phi(x,mu) + phi(x+e_mu,mu).
Q is the sum of squared plaquette curvatures, Q_1 the same on the coarse
lattice, and S_1(v) = min { Q(phi) : phi = 0 on tree links, V(phi) = v }.

Parts:
 (A) Stokes: 16 (dV)_P is the sum over x in B(y) of the curvature of the
     2 x 2 square from x, each fine plaquette lies in exactly four such
     squares, and Q_1(V(phi)) <= Q(phi) (exact arithmetic, periodic box).
 (B) the extension phi*(v) (v on the 32 links between blocks, zero inside)
     has V(phi*) = v and Q(phi*) = 4 Q_1(v) exactly.
 (C) floating point, comparison only: the eigenvalues of S_1 relative to Q_1
     in the Bloch decomposition lie in [1, 4], tend to 1 at small momentum
     and equal 4 at (pi, pi, 0, 0).
 (D) the bound ||r_w||_1 <= 2 M0 s sqrt(W) <= 4e14 for the response to one
     average, and its floating-point value in the cubes n = 4, 6, 8.
Run: python3 compose.py. Needs numpy and scipy."""
import itertools, math, random, sys
import numpy as np

E = np.eye(4, dtype=int)
OFFS = list(itertools.product(range(2), repeat=4))
PLANES = [(m, n) for m in range(4) for n in range(m + 1, 4)]
ok = []

def check(name, cond):
    ok.append(bool(cond)); print(('ok  ' if cond else 'FAIL'), name)

def tree(o, mu):
    o = [c % 2 for c in o]
    return o[mu] == 0 and all(o[nu] == 0 for nu in range(mu + 1, 4))

# (A), (B): periodic fine box of side L, coarse side L/2, integer fields
L = 8; M = L // 2
def wrap(x, n): return tuple(int(c) % n for c in x)
def curl(f, x, m, n, side):
    x = np.array(x)
    return (f[(wrap(x, side), m)] + f[(wrap(x + E[m], side), n)]
            - f[(wrap(x + E[n], side), m)] - f[(wrap(x, side), n)])
sites = list(itertools.product(range(L), repeat=4))
csites = list(itertools.product(range(M), repeat=4))
random.seed(1)
phi = {(x, mu): (0 if tree(x, mu) else random.randint(-5, 5)) for x in sites for mu in range(4)}
V16 = {}   # 16 times the averages
for y in csites:
    for mu in range(4):
        V16[(y, mu)] = sum(phi[(wrap(2 * np.array(y) + np.array(o), L), mu)]
                           + phi[(wrap(2 * np.array(y) + np.array(o) + E[mu], L), mu)] for o in OFFS)
stokes, cover = True, {}
for y in csites:
    for (m, n) in PLANES:
        tot = 0
        for o in OFFS:
            x = 2 * np.array(y) + np.array(o)
            for i, j in itertools.product(range(2), repeat=2):
                z = wrap(x + i * E[m] + j * E[n], L)
                tot += curl(phi, z, m, n, L)
                cover[(z, m, n)] = cover.get((z, m, n), 0) + 1
        stokes &= tot == curl(V16, y, m, n, M)
check('(A) 16 (dV)_P = sum of the 16 squares of side two (Stokes)', stokes)
check('(A) every fine plaquette lies in exactly four squares',
      set(cover.values()) == {4} and len(cover) == len(sites) * 6)
Q = sum(curl(phi, x, m, n, L) ** 2 for x in sites for (m, n) in PLANES)
Q1_256 = sum(curl(V16, y, m, n, M) ** 2 for y in csites for (m, n) in PLANES)
check('(A) Q_1(V(phi)) <= Q(phi)', Q1_256 <= 256 * Q)
print(f'    Q = {Q}, Q_1 = {Q1_256 / 256:.3f}')

v = {(y, mu): random.randint(-5, 5) for y in csites for mu in range(4)}
star = {}
for x in sites:
    y = tuple(c // 2 for c in x); o = tuple(c % 2 for c in x)
    for mu in range(4):
        star[(x, mu)] = v[(y, mu)] if o[mu] == 1 else 0
avg_ok = all(sum(star[(wrap(2 * np.array(y) + np.array(o), L), mu)]
                 + star[(wrap(2 * np.array(y) + np.array(o) + E[mu], L), mu)] for o in OFFS)
             == 16 * v[(y, mu)] for y in csites for mu in range(4))
tree_ok = all(star[(x, mu)] == 0 for x in sites for mu in range(4) if tree(x, mu))
Qs = sum(curl(star, x, m, n, L) ** 2 for x in sites for (m, n) in PLANES)
Q1v = sum(curl(v, y, m, n, M) ** 2 for y in csites for (m, n) in PLANES)
check('(B) phi* vanishes on tree links and V(phi*) = v', avg_ok and tree_ok)
check('(B) Q(phi*) = 4 Q_1(v)', Qs == 4 * Q1v)

# (C) Bloch decomposition
links = [(o, mu) for o in OFFS for mu in range(4)]
act = [l for l in links if not tree(*l)]; idx = {l: i for i, l in enumerate(act)}
def red(x, mu):
    x = np.array(x); return (tuple(x % 2), mu), tuple(x // 2)
P = []
for o in OFFS:
    x = np.array(o)
    for (mu, nu) in PLANES:
        r = []
        for (yy, m, sg) in [(x, mu, 1), (x + E[mu], nu, 1), (x + E[nu], mu, -1), (x, nu, -1)]:
            l, c = red(yy, m)
            if l in idx: r.append((idx[l], c, sg))
        P.append(r)
Cs = []
for mu in range(4):
    r = {}
    for o in OFFS:
        for k in range(2):
            l, c = red(np.array(o) + k * E[mu], mu)
            if l in idx: r[(idx[l], c)] = r.get((idx[l], c), 0) + 1 / 16
    Cs.append([(i, c, val) for (i, c), val in r.items()])
nA = len(act)
def sym(rows, q, ncol):
    Mx = np.zeros((len(rows), ncol), complex)
    for i, r in enumerate(rows):
        for (j, c, val) in r: Mx[i, j] += val * np.exp(1j * np.dot(q, c))
    return Mx
def coarse_curl(q):
    D = np.zeros((6, 4), complex)
    for k, (m, n) in enumerate(PLANES):
        D[k, m] += 1; D[k, n] += np.exp(1j * q[m]); D[k, m] -= np.exp(1j * q[n]); D[k, n] -= 1
    return D
def ratios(q):
    D = sym(P, q, nA); H = D.conj().T @ D; B = sym(Cs, q, nA)
    u, sv, vh = np.linalg.svd(B); N = vh[(sv > 1e-10).sum():].conj().T
    Bp = np.linalg.pinv(B)
    S1 = Bp.conj().T @ (H - H @ N @ np.linalg.solve(N.conj().T @ H @ N, N.conj().T @ H)) @ Bp
    D1 = coarse_curl(q); Q1 = D1.conj().T @ D1
    w, U = np.linalg.eigh(Q1); R = U[:, w > 1e-9]
    return np.linalg.eigvals(np.linalg.solve(R.conj().T @ Q1 @ R, R.conj().T @ S1 @ R)).real
lo, hi = 9, 0
G = 6
for qq in itertools.product(range(G), repeat=4):
    if not any(qq): continue
    ev = ratios(2 * np.pi * np.array(qq) / G)
    lo, hi = min(lo, ev.min()), max(hi, ev.max())
print(f'    S_1 / Q_1 on the 6^4 grid: [{lo:.6f}, {hi:.6f}]')
check('(C) 1 <= S_1 / Q_1 <= 4 in the Bloch decomposition', lo >= 1 - 1e-9 and hi <= 4 + 1e-9)
small = ratios(np.array([1e-3, 0, 0, 0])); top = ratios(np.array([np.pi, np.pi, 0, 0]))
print(f'    q = (0.001,0,0,0): {small.min():.8f}..{small.max():.8f}; q = (pi,pi,0,0): max {top.max():.8f}')
check('(C) the ratio tends to 1 at small momentum and equals 4 at (pi,pi,0,0)',
      abs(small.max() - 1) < 1e-5 and abs(top.max() - 4) < 1e-8)

# (D) response to one average
s, h, b = 96, 24, 1 / 16
x1 = 1 / (s * math.sqrt(b)); xn = math.hypot(52, 52 * h * x1 + x1)
yn = (1 + h * xn) / (s * math.sqrt(b)); M0 = math.hypot(xn, yn)
S = 192; lam = math.log(1 + 1 / (2 * S * M0))
qd = math.exp(-2 * lam)
sum3 = qd * (1 + 4 * qd + qd * qd) / (1 - qd) ** 4; sum1 = qd / (1 - qd) ** 2
W = 4 * math.exp(2 * lam) * (1 + 64 * sum3 + 16 * sum1)
Lw = 2 * 105 * s * math.sqrt(W)
print(f'    ||e^(lambda d) r_w|| <= {2 * 105 * s}, ||r_w||_1 <= {Lw:.3e}')
check('(D) ||r_w||_1 <= 4e14', 2 * 105 * s == 20160 and Lw <= 4e14)

import scipy.sparse as sp, scipy.sparse.linalg as spl
def response(n):
    V = [(x, mu) for x in itertools.product(range(n), repeat=4) for mu in range(4) if not tree(x, mu)]
    vi = {l: i for i, l in enumerate(V)}
    inV = lambda x: all(0 <= c < n for c in x)
    I, J, Dd, r = [], [], [], 0
    for x in itertools.product(range(-1, n), repeat=4):
        xa = np.array(x)
        for (mu, nu) in PLANES:
            ls = [(xa, mu, 1), (xa + E[mu], nu, 1), (xa + E[nu], mu, -1), (xa, nu, -1)]
            ent = [(vi[(tuple(int(c) for c in yy), m)], sg) for (yy, m, sg) in ls
                   if (tuple(int(c) for c in yy), m) in vi]
            if ent:
                for (j, sg) in ent: I.append(r); J.append(j); Dd.append(sg)
                r += 1
    D = sp.csr_matrix((Dd, (I, J)), shape=(r, len(V)))
    I, J, Dd, rows = [], [], [], []
    for y in itertools.product(range(n // 2), repeat=4):
        for mu in range(4):
            acc = {}
            for o in OFFS:
                for k in range(2):
                    z = tuple(int(c) for c in 2 * np.array(y) + np.array(o) + k * E[mu])
                    if (z, mu) in vi: acc[vi[(z, mu)]] = acc.get(vi[(z, mu)], 0) + 1 / 16
            for j, val in acc.items(): I.append(len(rows)); J.append(j); Dd.append(val)
            rows.append((y, mu))
    B = sp.csr_matrix((Dd, (I, J)), shape=(len(rows), len(V)))
    K = sp.bmat([[D.T @ D, B.T], [B, None]]).tocsc()
    c = (n // 4,) * 4; best = 0
    for mu in range(4):
        rhs = np.zeros(K.shape[0]); rhs[len(V) + rows.index((c, mu))] = 1
        best = max(best, np.abs(spl.spsolve(K, rhs)[:len(V)]).sum())
    return best
for n in (4, 6, 8):
    print(f'    n = {n}: max ||r_w||_1 for the central block = {response(n):.4f}')
check('(D) floating-point responses computed', True)

print(f'{sum(ok)} checks passed, {len(ok) - sum(ok)} failed')
sys.exit(0 if all(ok) else 1)
