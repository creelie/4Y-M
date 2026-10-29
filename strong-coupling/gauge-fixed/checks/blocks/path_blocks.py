"""Checks for the path blocks: the anchored cube, the coercivity bound 224 on
group-valued fields, the coarse curvature and reconstruction bounds, and the
thin-path ratio.

Parts:
 (A) anchored cube: the Laplacian of {0,1}^3 with the vertex 0 deleted has
     characteristic polynomial (x-3)(x^2-6x+2)(x-2)^2(x-4)^2 (exact), so its
     eigenvalues are 3, 3 +- sqrt 7, and 2, 4 twice; the least, 3 - sqrt 7,
     exceeds 1/3.
 (B) SU(2) fields in the tree gauge on the periodic lattice of side 4:
     D_v^2 <= 224 Q, Q_1(V) <= 4 Q, Q(R(v)) = 4 Q_1(v), Q <= 8 Q_1 + 48 D_v^2
     (floating point, random fields at five scales).
 (C) thin paths: Q / Q_k = s sin^2(pi/(2s)) for s = 2, 4, ..., 64, from the
     s/2 changed links with six plaquettes each and the six coarse plaquettes
     of the one nontrivial coarse link.
 (D) slice windings: three covariant continuous maps of eight slice links
     (u_1 and alternating products u_a u_b^-1 u_c ...) satisfy F(g u g^-1) =
     g F(u) g^-1 and F(g u) = g F(u), and the degrees of their coordinate loops
     add up to one; the normalised quaternion mean of the eight links, which
     treats them alike, takes values 2 apart at fields arbitrarily close to
     the field with u_t = -I on four links and I on the others.
 (E) the same for SU(3) on its maximal torus: for covariant continuous maps of
     six slice links the integer matrices W_t on H_1 of the torus add up to the
     identity.
 (F) the cut-off average V# = V^p exp(chi/16 sum_x log(V^p^-1 P_x)) on SU(2),
     rho = 1, c = pi/3: covariance, differential at I equal to the block
     average, the bounds 896, 56 c^2, 448 + 896 c^2, 8 + 2688 c^2 on random
     fields at five scales, slice windings 1 at the path and 0 elsewhere,
     and the coarse link e^{X/8} of a single small link e^X.
Run: python3 path_blocks.py. Needs numpy."""
import itertools, math
import numpy as np

ok = []
def check(name, cond):
    ok.append(bool(cond)); print(('ok  ' if cond else 'FAIL'), name)

# (A)
V = [a for a in itertools.product(range(2), repeat=3) if any(a)]
L0 = np.zeros((7, 7), dtype=int)
for i, a in enumerate(V):
    L0[i, i] = 3
    for j, b in enumerate(V):
        if sum(x != y for x, y in zip(a, b)) == 1: L0[i, j] = -1
ev = np.sort(np.linalg.eigvalsh(L0))
want = np.sort([3, 3 + math.sqrt(7), 3 - math.sqrt(7), 2, 2, 4, 4])
check('(A) spectrum of the grounded cube Laplacian is 3, 3 +- sqrt7, 2, 2, 4, 4',
      np.allclose(ev, want))
# exact: det(L0 - 3) in integers equals the value of (t^3-7t)(t^2-6t+8)^2 at t=0... use
# the integer characteristic polynomial via Faddeev-LeVerrier
def charpoly(M):
    n = len(M); M = [[int(x) for x in r] for r in M]
    from fractions import Fraction as F
    c = [F(1)]; Mk = [[F(0)] * n for _ in range(n)]
    I = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    A = [[F(x) for x in r] for r in M]
    for k in range(1, n + 1):
        Mk = [[sum(A[i][l] * Mk[l][j] for l in range(n)) + c[-1] * I[i][j] for j in range(n)] for i in range(n)]
        AM = [[sum(A[i][l] * Mk[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        c.append(-sum(AM[i][i] for i in range(n)) / k)
    return [int(x) for x in c]
# expected (x-3)(x^2-6x+2)(x-2)^2(x-4)^2
def mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return r
exp = [1]
for f in ([1, -3], [1, -6, 2], [1, -2], [1, -2], [1, -4], [1, -4]): exp = mul(exp, f)
check('(A) exact characteristic polynomial (x-3)(x^2-6x+2)(x-2)^2(x-4)^2', charpoly(L0) == exp)
check('(A) 3 - sqrt7 > 1/3, i.e. 64 > 63 after squaring 8/3 > sqrt7', 64 > 9 * 7)

# (B)
rng = np.random.default_rng(7)
S = 4; d = 4; E = np.eye(d, dtype=int)
sig = np.array([[[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]])
def su2(a):
    n = np.linalg.norm(a)
    if n == 0: return np.eye(2, dtype=complex)
    return math.cos(n) * np.eye(2) + math.sin(n) / n * 1j * np.tensordot(a, sig, 1)
def tree(o, mu): return o[mu] == 0 and all(o[n] == 0 for n in range(mu + 1, d))
sites = list(itertools.product(range(S), repeat=d))
sh = lambda x, mu, k=1: tuple((np.array(x) + k * E[mu]) % S)
def hs2(A): return float(np.sum(np.abs(A) ** 2))
def energies(U, Sz=S, step=1):
    Q = 0.0
    pts = [x for x in itertools.product(range(0, Sz, step), repeat=d)]
    for x in pts:
        for mu in range(d):
            for nu in range(mu + 1, d):
                P = U[(x, mu)] @ U[(tuple((np.array(x) + step * E[mu]) % Sz), nu)] @ \
                    U[(tuple((np.array(x) + step * E[nu]) % Sz), mu)].conj().T @ U[(x, nu)].conj().T
                Q += hs2(P - np.eye(2))
    return Q
worst = [0, 0, 0, 0]
for sc in [0.05, 0.3, 1.0, 2.0, 6.0]:
    for rep in range(2):
        U = {}
        for x in sites:
            o = [c % 2 for c in x]
            for mu in range(d):
                U[(x, mu)] = np.eye(2, dtype=complex) if tree(o, mu) else su2(sc * rng.normal(size=3))
        Q = energies(U)
        v = {(y, mu): U[(sh(y, mu), mu)] for y in itertools.product(range(0, S, 2), repeat=d) for mu in range(d)}
        R, D = {}, 0.0
        for x in sites:
            o = [c % 2 for c in x]; y = tuple(c // 2 * 2 for c in x)
            for mu in range(d):
                R[(x, mu)] = v[(y, mu)] if o[mu] == 1 else np.eye(2, dtype=complex)
                D += hs2(U[(x, mu)] - R[(x, mu)])
        # coarse plaquettes: holonomy of 2x2 squares at even sites (coarse torus side 2 counts each twice per plane)
        Vc = {}
        for y in itertools.product(range(0, S, 2), repeat=d):
            for mu in range(d): Vc[(y, mu)] = U[(y, mu)] @ U[(sh(y, mu), mu)]
        Q1 = energies(Vc, S, 2)
        QR = energies(R)
        worst[0] = max(worst[0], D / Q); worst[1] = max(worst[1], Q1 / Q)
        worst[2] = max(worst[2], abs(QR - 4 * Q1) / max(1, Q1)); worst[3] = max(worst[3], Q / (8 * Q1 + 48 * D))
print('    max D^2/Q = %.4f, max Q_1/Q = %.4f, max |Q(R)-4Q_1| = %.2e, max Q/(8Q_1+48D^2) = %.4f' % tuple(worst))
check('(B) D_v^2 <= 224 Q on random SU(2) fields', worst[0] <= 224)
check('(B) Q_1(V) <= 4 Q', worst[1] <= 4)
check('(B) Q(R(v)) = 4 Q_1(v)', worst[2] < 1e-9)
check('(B) Q <= 8 Q_1 + 48 D_v^2', worst[3] <= 1)

# (C)
T = np.diag([1j * math.pi / 2, -1j * math.pi / 2])
def expm(Hd): return np.diag(np.exp(np.diag(Hd)))
good = True
for k in range(1, 7):
    s = 2 ** k; n = 3 * s
    ch = expm(2 * T / s)
    # plaquettes containing a changed link (j,0,0,0), direction 0, j odd < s: six each, distinct
    Qf = 0.0
    for j in range(1, s, 2):
        for nu in range(1, 4):
            # plaquette at (j,0,..) and at (j,0,..)-e_nu in the plane {0,nu}
            Qf += 2 * hs2(ch - np.eye(2))
    Qk = 6 * hs2(expm(T) - np.eye(2))
    r = Qf / Qk; want = s * math.sin(math.pi / (2 * s)) ** 2
    good &= abs(r - want) < 1e-12 and r <= math.pi ** 2 / (4 * s) + 1e-15
    print('    s = %3d: Q/Q_k = %.6f, s sin^2(pi/2s) = %.6f' % (s, r, want))
check('(C) thin-path ratio s sin^2(pi/(2s)) <= pi^2/(4s) for s = 2..64', good)

# (D)
X = np.array([[1j, 0], [0, -1j]])
def ex(th): return np.diag([np.exp(1j * th), np.exp(-1j * th)])
def angle(A): return math.atan2(A[0, 0].imag, A[0, 0].real)
def degree(F, t, n=8, N=400):
    tot, prev = 0.0, None
    for k in range(N + 1):
        u = [np.eye(2, dtype=complex)] * n; u = list(u); u[t] = ex(2 * math.pi * k / N)
        a = angle(F(u))
        if prev is not None:
            d = a - prev; d -= 2 * math.pi * round(d / (2 * math.pi)); tot += d
        prev = a
    return round(tot / (2 * math.pi))
inv = lambda A: A.conj().T
maps = {'path map u_1': lambda u: u[0],
        'u_1 u_2^-1 u_3': lambda u: u[0] @ inv(u[1]) @ u[2],
        'u_4 u_7^-1 u_2 u_5^-1 u_8': lambda u: u[3] @ inv(u[6]) @ u[1] @ inv(u[4]) @ u[7]}
rng = np.random.default_rng(7)
def haar():
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    return np.array([[q[0] + 1j * q[1], q[2] + 1j * q[3]], [-q[2] + 1j * q[3], q[0] - 1j * q[1]]])
good = True
for F in maps.values():
    for _ in range(20):
        u = [haar() for _ in range(8)]; g = haar()
        good &= np.allclose(F([g @ a @ inv(g) for a in u]), g @ F(u) @ inv(g))
        good &= np.allclose(F([g @ a for a in u]), g @ F(u))
check('(D) the three maps are covariant on slice fields', good)
good = True
for name, F in maps.items():
    w = [degree(F, t) for t in range(8)]
    print('    %-30s w = %s' % (name, w)); good &= sum(w) == 1
check('(D) slice windings add up to one for covariant continuous maps', good)
def qmean(u):
    q = sum(u) / len(u); r = math.sqrt(abs(np.linalg.det(q)))
    return q / r
eps = 1e-6
Y = np.array([[0, 1], [-1, 0]], dtype=complex)
def near(Xd):
    return [ex(math.pi) @ (np.eye(2) * math.cos(eps) + Xd * math.sin(eps))] * 4 + [np.eye(2, dtype=complex)] * 4
jump = np.linalg.norm(qmean(near(X)) - qmean(near(Y)))
print('    quaternion mean at fields within %g of (-I)^4 I^4 in two directions differs by %.4f' % (eps, jump))
check('(D) the symmetric quaternion mean jumps by 2 near (-I,...,-I,I,...,I)', abs(jump - 2) < 1e-4)
# (E)
def tor3(a, b): return np.diag([np.exp(1j * a), np.exp(1j * b), np.exp(-1j * (a + b))])
def wind(vals):
    tot = 0.0
    for v0, v1 in zip(vals, vals[1:]):
        d = math.atan2((v1 / v0).imag, (v1 / v0).real); tot += d
    return round(tot / (2 * math.pi))
def wmatrix(F, t, n=6, N=400):
    W = np.zeros((2, 2), dtype=int)
    for j in range(2):
        d0, d1 = [], []
        for k in range(N + 1):
            u = [np.eye(3, dtype=complex) for _ in range(n)]
            th = 2 * math.pi * k / N
            u[t] = tor3(th, 0) if j == 0 else tor3(0, th)
            R = F(u); d0.append(R[0, 0]); d1.append(R[1, 1])
        W[0, j], W[1, j] = wind(d0), wind(d1)
    return W
def haar3():
    z = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    q, r = np.linalg.qr(z); q = q @ np.diag(np.diag(r) / abs(np.diag(r)))
    return q / np.linalg.det(q) ** (1 / 3)
maps3 = {'u_1': lambda u: u[0], 'u_1 u_2^-1 u_3': lambda u: u[0] @ inv(u[1]) @ u[2],
         'u_6 u_4^-1 u_5 u_2^-1 u_3': lambda u: u[5] @ inv(u[3]) @ u[4] @ inv(u[1]) @ u[2]}
good = True
for F in maps3.values():
    for _ in range(10):
        u = [haar3() for _ in range(6)]; g = haar3()
        good &= np.allclose(F([g @ a @ inv(g) for a in u]), g @ F(u) @ inv(g))
        good &= np.allclose(F([g @ a for a in u]), g @ F(u))
check('(E) the three SU(3) maps are covariant on slice fields', good)
good = True
for name, F in maps3.items():
    S = sum(wmatrix(F, t) for t in range(6))
    print('    %-26s sum of W_t = %s' % (name, S.tolist())); good &= (S == np.eye(2, dtype=int)).all()
check('(E) for SU(3) the torus matrices W_t add up to the identity', good)
# (F)
Sd = 4; I2 = np.eye(2, dtype=complex); c = math.pi / 3; rho = 1.0
shf = lambda x, mu, k=1: tuple((np.array(x) + k * E[mu]) % Sd)
def slog(g):  # principal logarithm of SU(2) near I
    n = math.acos(max(-1.0, min(1.0, g[0, 0].real / 2 + g[1, 1].real / 2)))
    A = (g - g.conj().T) / 2
    return A if n < 1e-12 else n / math.sin(n) * A
def sexp(Z):
    n = math.sqrt(hs2(Z) / 2)
    return I2 * math.cos(n) + (Z if n < 1e-12 else math.sin(n) / n * Z)
def phi(t): return 1.0 if t <= 0.5 else (0.0 if t >= 1 else 2 * (1 - t))
blk = list(itertools.product(range(2), repeat=d))
def sharp(U, y, mu):
    P = [U.get((tuple((np.array(y) + o) % Sd), mu), I2) @
         U.get((shf(tuple((np.array(y) + o) % Sd), mu), mu), I2) for o in blk]
    v = [inv(P[0]) @ p for p in P]
    m = max(math.sqrt(hs2(w - I2)) for w in v); chi = phi(m / rho)
    M = I2 if chi == 0 else sexp(chi * sum(slog(w) for w in v) / 16)
    return P[0] @ M, P[0], M, P
evens = list(itertools.product(range(0, Sd, 2), repeat=d))
def field(sc):
    U = {}
    for x in sites:
        o = [q % 2 for q in x]
        for mu in range(d):
            U[(x, mu)] = I2 if tree(o, mu) else su2(sc * rng.normal(size=3))
    return U
# covariance under the gauge transformations constant on blocks
good = True
for sc in [0.05, 0.3, 2.0]:
    U = field(sc); h = {y: haar() for y in evens}
    hb = lambda x: h[tuple(q // 2 * 2 for q in x)]
    Uh = {(x, mu): hb(x) @ U[(x, mu)] @ inv(hb(shf(x, mu))) for (x, mu) in U}
    for y in evens:
        for mu in range(d):
            good &= np.allclose(sharp(Uh, y, mu)[0], h[y] @ sharp(U, y, mu)[0] @ inv(h[shf(y, mu, 2)]))
check('(F) the cut-off average is covariant under block gauge transformations', good)
# linearisation: the block average
t = 1e-6; err = 0.0
for rep in range(3):
    a = {}
    for x in sites:
        o = [q % 2 for q in x]
        for mu in range(d):
            a[(x, mu)] = np.zeros(3) if tree(o, mu) else rng.normal(size=3)
    U = {k: su2(t * a[k]) for k in a}
    for y in evens:
        for mu in range(d):
            avg = sum(a[(tuple((np.array(y) + o) % Sd), mu)] + a[(shf(tuple((np.array(y) + o) % Sd), mu), mu)]
                      for o in blk) / 16
            lin = (sharp(U, y, mu)[0] - I2) / t
            err = max(err, np.abs(lin - 1j * np.tensordot(avg, sig, 1)).max())
print('    max |(V#(e^{t phi}) - I)/t - V(phi)| at t = 1e-6: %.2e' % err)
check('(F) the differential of the cut-off average at I is the block average', err < 1e-4)
# bounds for all fields
worst = [0, 0, 0, 0]
for sc in [0.05, 0.3, 1.0, 2.0, 6.0]:
    for rep in range(2):
        U = field(sc); Q = energies(U, Sd)
        Vs, sP, sM = {}, 0.0, 0.0
        for y in evens:
            for mu in range(d):
                V, P0, M, P = sharp(U, y, mu); Vs[(y, mu)] = V
                sP += sum(hs2(p - P0) for p in P); sM += hs2(M - I2)
        D = 0.0
        for x in sites:
            o = [q % 2 for q in x]; y = tuple(q // 2 * 2 for q in x)
            for mu in range(d):
                D += hs2(U[(x, mu)] - (Vs[(y, mu)] if o[mu] == 1 else I2))
        Q1 = energies(Vs, Sd, 2)
        worst = [max(worst[0], sP / Q), max(worst[1], sM / Q), max(worst[2], D / Q), max(worst[3], Q1 / Q)]
print('    max sum d(P_x,P_2y)^2/Q = %.3f, sum d(M,I)^2/Q = %.4f, D^2/Q = %.3f, Q_1/Q = %.3f' % tuple(worst))
check('(F) sum of d(P_x, V^p)^2 <= 896 Q', worst[0] <= 896)
check('(F) sum of d(V#, V^p)^2 <= 56 c^2 Q, c = pi/3', worst[1] <= 56 * c * c)
check('(F) D_{V#}^2 <= (448 + 896 c^2) Q', worst[2] <= 448 + 896 * c * c)
check('(F) Q_1(V#) <= (8 + 2688 c^2) Q', worst[3] <= 8 + 2688 * c * c)
# slice windings: the identity at the path, zero elsewhere
trans = [a for a in itertools.product(range(2), repeat=3)]
def Fsharp(u):
    U = {((1,) + a, 0): u[i] for i, a in enumerate(trans)}
    return sharp(U, (0, 0, 0, 0), 0)[0]
w = [degree(Fsharp, t) for t in range(8)]
print('    slice windings of the cut-off average: %s' % w)
check('(F) slice windings of the cut-off average are 1 at the path and 0 elsewhere', w == [1] + [0] * 7)
# one link e^X on the straight path or at another transverse position gives e^{X/8}
good = True
for a in trans:
    for th in [0.05, 0.15]:
        X = th * np.array([[1j, 0], [0, -1j]])
        V = sharp({((1,) + a, 0): sexp(X)}, (0, 0, 0, 0), 0)[0]
        good &= np.allclose(V, sexp(X / 8))
check('(F) one small link e^X at any transverse position gives the coarse link e^{X/8}', good)
print('%d checks passed, %d failed' % (sum(ok), len(ok) - sum(ok)))
