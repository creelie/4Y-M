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
 (G) torus-valued fields (random and smooth, under random gauge
     transformations) on the lattice of side 8, two steps with the tree gauge
     restored between them: the cut-off is inactive, the coarse plaquette
     angles are the means of the sixteen 2x2 fluxes, the L1 norm drops by at
     least 4, the L2 norm does not increase, the maximum grows by at most 4,
     and Q_2 <= 1.001 Q when eta_0 A <= (rho/12)^2.
 (H) SU(2) fields of small noncommuting curvature (random and smooth, under
     random SU(2) gauge transformations) on the lattice of side 8, two steps:
     the cut-off is inactive, the coarse plaquette logarithm is the mean of the
     sixteen 2x2 square logarithms up to 3.2 r_P^2 with r_P <= 28.26 eta, and
     the growth bounds of sqrt E, A and eta of the nonabelian step hold.
 (I) large fields under the SU(2) Wilson measure: on the torus of side 8 a
     plaquette lies in exactly four unit cubes, a cube has 24 plaquettes and
     80 neighbours, Z_{y,mu} has 24 cubes and a cube lies in Z_{y,mu} for at
     most eight coarse links; the constants C_0, 4 8 4 24 144 <= 4.5e5 and
     8 b_0 = 11/(3 pi^2) < 4/3; the link-angle bound 1/(24 pi^3 beta^{3/2}) and
     the one-link tail p = v_* e^{-3 beta} by quadrature, with the cap
     probability decreasing in the distance of the centre and increasing in
     kappa; the staple sum -4I of the frustrated exterior on a lattice and the
     conditional probability of p_0 being bad, even for rho = 2.5; and on random SU(2) fields with
     localised large links, chi = 1 at every coarse link whose 24 cubes are
     good for rho = 1/12.
 (J) defects across scales: SU(2) fields of small curvature with a localised
     large link on the lattice of side 8, one step: every plaquette outside
     the propagated defect set B_j (coarse plaquettes whose box 2y+[0,4]^4
     meets B_{j-1}) has |Theta| <= eta_j with eta_{j+1} = 4(1+647 eta_j)eta_j,
     every plaquette of B_j lies within l-infinity distance 3(2^j-1) of B_0,
     and the energy outside the defects obeys the bound of the theorem.
 (K) effective forms uniformly in the scale: in exact integer arithmetic on
     periodic lattices (d = 4 with L = 2, 4; d = 2 with L = 8; d = 3 with
     L = 4), the Whitney extension with nodes at the block centres has k-fold
     averages T u with T = tau^(x d), tau = (1/8, 3/4, 1/8), plaquette energy
     at most Q_k(u), and zero curvature on gradients; the symbol of tau is at
     least 1/2, against 1/L with the coarse sites as nodes.
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
# (G)
def tor(t): return np.diag([np.exp(1j * t), np.exp(-1j * t)])
def ang(g): return math.atan2(g[0, 0].imag, g[0, 0].real)
def shz(x, mu, Sz, k=1): return tuple((np.array(x) + k * E[mu]) % Sz)
def gfix(U, Sz):
    h = {}
    for y in itertools.product(range(0, Sz, 2), repeat=d):
        for o in sorted(blk, key=sum):
            x = tuple((np.array(y) + o) % Sz)
            if not any(o): h[x] = I2; continue
            nu = max(i for i in range(d) if o[i]); p = list(o); p[nu] = 0
            xp = tuple((np.array(y) + p) % Sz); h[x] = h[xp] @ U[(xp, nu)]
    return {(x, mu): h[x] @ U[(x, mu)] @ inv(h[shz(x, mu, Sz)]) for (x, mu) in U}
def sharpz(U, y, mu, Sz):
    P = []
    for o in blk:
        x = tuple((np.array(y) + o) % Sz); P.append(U[(x, mu)] @ U[(shz(x, mu, Sz), mu)])
    v = [inv(P[0]) @ p for p in P]
    m = max(math.sqrt(hs2(w - I2)) for w in v); chi = phi(m / rho)
    M = I2 if chi == 0 else sexp(chi * sum(slog(w) for w in v) / 16)
    return P[0] @ M, chi
def plaq(U, x, mu, nu, Sz):
    return U[(x, mu)] @ U[(shz(x, mu, Sz), nu)] @ inv(U[(shz(x, nu, Sz), mu)]) @ inv(U[(x, nu)])
def angles(U, Sz):
    return {(x, mu, nu): ang(plaq(U, x, mu, nu, Sz)) for x in itertools.product(range(Sz), repeat=d)
            for mu in range(d) for nu in range(mu + 1, d)}
def step(U, Sz):
    U = gfix(U, Sz); W, cmin = {}, 1.0
    for y in itertools.product(range(0, Sz, 2), repeat=d):
        for mu in range(d):
            W[(tuple(q // 2 for q in y), mu)], c = sharpz(U, y, mu, Sz); cmin = min(cmin, c)
    return W, cmin
Sz = 8; eta_G = rho / 12
good_formula, good_tree, rows = True, True, []
cph = rng.uniform(0, 2 * math.pi, size=(d, d))
for eps in [1e-5, 3e-3, 0.004]:
    U = {}
    g = {x: tor(rng.uniform(0, 2 * math.pi)) for x in itertools.product(range(Sz), repeat=d)}
    for x in itertools.product(range(Sz), repeat=d):
        for mu in range(d):
            a = eps * rng.normal() if eps != 0.004 else \
                eps * sum(math.sin(2 * math.pi * x[n] / Sz + cph[mu, n]) for n in range(d) if n != mu)
            U[(x, mu)] = g[x] @ tor(a) @ inv(g[shz(x, mu, Sz)])
    Th = angles(U, Sz); Q0 = sum(8 * math.sin(t / 2) ** 2 for t in Th.values())
    A = sum(math.sqrt(2) * abs(t) for t in Th.values()); eta0 = max(math.sqrt(2) * abs(t) for t in Th.values())
    cur, cs = U, Sz
    for j in range(2):
        Ug = gfix(cur, cs)
        good_tree &= all(np.allclose(Ug[(tuple((np.array(y) + o) % cs), mu)], I2)
                         for y in itertools.product(range(0, cs, 2), repeat=d) for o in blk for mu in range(d)
                         if tree(list(o), mu))
        Tc = angles(cur, cs); mx = max(math.sqrt(2) * abs(t) for t in Tc.values())
        nxt, cmin = step(cur, cs); ns = cs // 2; Tn = angles(nxt, ns)
        err = 0.0
        if mx <= eta_G:
            for (y, mu, nu), t in Tn.items():
                f = 0.0
                for o in blk:
                    x = tuple((2 * np.array(y) + o) % cs)
                    for a, b in [(0, 0), (1, 0), (0, 1), (1, 1)]:
                        z = tuple((np.array(x) + a * E[mu] + b * E[nu]) % cs); f += Tc[(z, mu, nu)]
                err = max(err, abs(t - f / 16))
            good_formula &= err < 1e-9 and cmin == 1.0
        L1c = sum(abs(t) for t in Tc.values()); L1n = sum(abs(t) for t in Tn.values())
        L2c = sum(t * t for t in Tc.values()); L2n = sum(t * t for t in Tn.values())
        mxn = max(abs(t) for t in Tn.values())
        rows.append((eps, j, mx <= eta_G, cmin, err, L1n / L1c, L2n / L2c, mxn / max(abs(t) for t in Tc.values())))
        cur, cs = nxt, ns
    Qk = sum(8 * math.sin(t / 2) ** 2 for t in angles(cur, cs).values())
    rows.append(('Q', eps, eta0 * A <= eta_G ** 2, Qk / Q0))
for r in rows:
    if r[0] == 'Q': print('    eps = %g: eta0*A <= eta_G^2: %s, Q_2/Q = %.4f' % r[1:])
    else: print('    eps = %g, step %d: hypothesis %s, min chi %.3f, formula err %.1e, L1 ratio %.4f, L2 ratio %.4f, sup ratio %.3f' % r)
check('(G) the tree gauge fixing puts every tree link at I', good_tree)
check('(G) torus fields: coarse angles are the means of the 2x2 fluxes, cut-off inactive', good_formula)
check('(G) L1 ratio <= 1/4, L2 ratio <= 1, sup ratio <= 4 at every step',
      all(r[5] <= 0.25 + 1e-12 and r[6] <= 1 + 1e-12 and r[7] <= 4 + 1e-12 for r in rows if r[0] != 'Q' and r[2]))
check('(G) Q_2(V#^2) <= 1.001 Q when eta0 * A <= eta_G^2',
      all(r[3] <= 1.001 for r in rows if r[0] == 'Q' and r[2]) and any(r[2] for r in rows if r[0] == 'Q'))
# (H)
def nrm(Z): return math.sqrt(hs2(Z))
def logs(U, Sz):
    return {(x, mu, nu): slog(plaq(U, x, mu, nu, Sz)) for x in itertools.product(range(Sz), repeat=d)
            for mu in range(d) for nu in range(mu + 1, d)}
def stats(L):
    v = [nrm(Z) for Z in L.values()]
    return sum(t * t for t in v), sum(v), max(v)
def pathp(U, x, mu, Sz): return U[(x, mu)] @ U[(shz(x, mu, Sz), mu)]
def nonab_step(cur, cs):
    Ug = gfix(cur, cs); nxt, cmin = step(cur, cs); ns = cs // 2
    Ln = logs(nxt, ns); worst, rmax = 0.0, 0.0
    for (yc, mu, nu), Th in Ln.items():
        y = 2 * np.array(yc)
        corner = [(y, mu), (y + 2 * E[mu], nu), (y + 2 * E[nu], mu), (y, nu)]
        K, Dmax = [], [0.0] * 4
        for o in blk:
            Pl = [pathp(Ug, tuple((c + o) % cs), m, cs) for c, m in corner]
            K.append(slog(Pl[0] @ Pl[1] @ inv(Pl[2]) @ inv(Pl[3])))
            for i, (c, m) in enumerate(corner):
                P0 = pathp(Ug, tuple(c % cs), m, cs)
                Dmax[i] = max(Dmax[i], nrm(slog(inv(P0) @ Pl[i])))
        Fbar = sum(K) / 16; r = sum(Dmax) + nrm(K[0])
        rmax = max(rmax, r); worst = max(worst, nrm(Th - Fbar) / (3.2 * r * r))
    return nxt, ns, cmin, worst, rmax
good_lem, good_steps, rows = True, True, []
for kind, eps in [('random', 4e-5), ('random', 7e-5), ('smooth', 1e-4)]:
    U = {}
    g = {x: haar() for x in itertools.product(range(Sz), repeat=d)}
    for x in itertools.product(range(Sz), repeat=d):
        for mu in range(d):
            a = eps * rng.normal(size=3) if kind == 'random' else \
                eps * np.array([sum(math.sin(2 * math.pi * x[n] / Sz + cph[mu, n] + k) for n in range(d) if n != mu)
                                for k in range(3)])
            U[(x, mu)] = g[x] @ su2(a) @ inv(g[shz(x, mu, Sz)])
    cur, cs = U, Sz
    for j in range(2):
        E0, A0, e0 = stats(logs(cur, cs))
        nxt, ns, cmin, worst, rmax = nonab_step(cur, cs)
        E1, A1, e1 = stats(logs(nxt, ns))
        ok_h = e0 <= rho / 1000
        good_lem &= (not ok_h) or (cmin == 1.0 and worst <= 1 and rmax <= 28.26 * e0)
        good_steps &= (not ok_h) or (math.sqrt(E1) <= (1 + 6030 * e0) * math.sqrt(E0) and
                                      A1 <= A0 / 4 * (1 + 9147 * e0) and e1 <= 4 * e0 * (1 + 647 * e0))
        rows.append((kind, eps, j, e0, ok_h, cmin, worst, math.sqrt(E1 / E0), 4 * A1 / A0, e1 / e0))
        cur, cs = nxt, ns
for r in rows:
    print('    %s eps = %g, step %d: eta = %.2e (<= rho/1000: %s), min chi %.0f, max err/(3.2 r^2) = %.3f, '
          'sqrt(E1/E0) = %.4f, 4 A1/A0 = %.4f, eta1/eta0 = %.3f' % r)
check('(H) nonabelian small fields: cut-off inactive and |log V#_P - Fbar_P| <= 3.2 r_P^2, r_P <= 28.26 eta', good_lem)
check('(H) sqrt E grows by at most 1 + 6030 eta, A by (1/4)(1 + 9147 eta), eta by 4(1 + 647 eta)', good_steps)
check('(H) the nonabelian error term is present (max err/(3.2 r^2) > 0)', max(r[6] for r in rows) > 1e-4)
# (I)
def quad(f, a, b, n=20000):
    t = np.linspace(a, b, n + 1); w = np.full(n + 1, 2.0); w[0] = w[-1] = 1; w[1:-1:2] = 4
    return float(np.sum(w * f(t)) * (b - a) / (3 * n))
Lt = 8
cubes = list(itertools.product(range(Lt), repeat=d))
def cube_plaqs(z):
    out = set()
    for w in itertools.product(range(2), repeat=d):
        for mu in range(d):
            for nu in range(mu + 1, d):
                if w[mu] == 0 and w[nu] == 0:
                    out.add((tuple((np.array(z) + w) % Lt), mu, nu))
    return out
cp = {z: cube_plaqs(z) for z in cubes}
inc = {}
for z, P in cp.items():
    for p in P: inc[p] = inc.get(p, 0) + 1
check('(I) a unit cube has 24 plaquettes and a plaquette lies in exactly 4 unit cubes',
      all(len(P) == 24 for P in cp.values()) and set(inc.values()) == {4} and len(inc) == 6 * Lt ** 4)
nb = len([w for w in itertools.product((-1, 0, 1), repeat=d) if any(w)])
check('(I) a unit cube has 80 neighbours at l-infinity distance one', nb == 80 and 80 ** 2 <= 6400)
def Zset(y, mu):
    return {tuple((2 * np.array(y) + t * E[mu] + np.array(w)) % Lt) for t in range(3)
            for w in itertools.product(range(2), repeat=d) if w[mu] == 0}
Zs = {(y, mu): Zset(y, mu) for y in itertools.product(range(Lt // 2), repeat=d) for mu in range(d)}
cnt = {}
for k, Z in Zs.items():
    for z in Z: cnt[z] = cnt.get(z, 0) + 1
check('(I) |Z_{y,mu}| = 24 and a cube lies in Z_{y,mu} for at most 8 coarse links',
      all(len(Z) == 24 for Z in Zs.values()) and max(cnt.values()) <= 8)
C0 = 0.5 + 2 / 3 * math.log(24 * math.pi ** 3)
b0 = 11 * 2 / (48 * math.pi ** 2)
print('    C_0 = %.4f, 4*8*4*24*144 = %d, 8 b_0 = %.4f' % (C0, 4 * 8 * 4 * 24 * 144, 8 * b0))
check('(I) 4*8*4*24*144 <= 4.5e5 and the one-loop slope 8 b_0 = 11/(3 pi^2) is below 4/3',
      4 * 8 * 4 * 24 * 144 <= 4.5e5 and abs(8 * b0 - 11 / (3 * math.pi ** 2)) < 1e-12 and 8 * b0 < 4 / 3)
good = True
for beta in [1, 3, 10, 100, 1e4]:
    r = 1 / (4 * math.sqrt(beta))
    pr = 2 / math.pi * quad(lambda t: np.sin(t) ** 2, 0, r)
    good &= pr >= 1 / (24 * math.pi ** 3 * beta ** 1.5)
check('(I) a Haar link has angle <= (4 sqrt beta)^-1 with probability >= 1/(24 pi^3 beta^{3/2})', good)
vstar = 2 / math.pi * quad(lambda t: np.sin(t) ** 2, math.acos(0.75), math.pi / 3)
good = True; rows = []
for beta in [0.5, 1, 2, 5, 10]:
    for rho_ in [1.0, 0.5, 1 / 12]:
        t0 = math.acos(1 - rho_ ** 2 / 4)
        for kap in np.linspace(0, 6 * beta, 7):
            m = kap  # stabilise the exponentials
            num = quad(lambda t: np.exp(kap * np.cos(t) - m) * np.sin(t) ** 2, t0, math.pi)
            den = quad(lambda t: np.exp(kap * np.cos(t) - m) * np.sin(t) ** 2, 0, math.pi)
            good &= num / den >= vstar * math.exp(-3 * beta) * (1 - 1e-9)
    rows.append((beta, vstar * math.exp(-3 * beta)))
print('    v_* = %.5f; p = %s' % (vstar, ', '.join('%.3e at beta %g' % (q, b) for b, q in rows)))
check('(I) the angle exceeds arccos(1 - rho^2/4) with probability >= v_* e^{-3 beta} for kappa <= 6 beta', good)
# cap about I under the density e^{kappa h.u}: decreasing in the angle phi of h, increasing in kappa
def capprob(phi_, kap, th0, n=241):
    th = np.linspace(0, math.pi, n); al = np.linspace(0, math.pi, n)
    T, A = np.meshgrid(th, al, indexing='ij')
    f = np.exp(kap * (math.cos(phi_) * np.cos(T) + math.sin(phi_) * np.sin(T) * np.cos(A) - 1)) * np.sin(T) ** 2 * np.sin(A)
    W = np.outer(np.where(th <= th0, 1.0, 0.0), np.ones(n))
    return float(np.trapezoid(np.trapezoid(f * W, al, axis=1), th) / np.trapezoid(np.trapezoid(f, al, axis=1), th))
good = True
for kap in [1.0, 4.0, 12.0]:
    vals = [capprob(ph, kap, 0.6) for ph in np.linspace(0, math.pi, 9)]
    good &= all(vals[i] >= vals[i + 1] - 1e-6 for i in range(len(vals) - 1))
kv = [capprob(0.0, k, 0.6) for k in [0.0, 1.0, 4.0, 12.0]]
good &= all(kv[i] <= kv[i + 1] + 1e-6 for i in range(3))
check('(I) the probability of a cap about I is largest at h = I and increases with kappa', good)
# frustrated exterior on a lattice of side 4
S4 = 4; x0 = (1, 1, 1, 1)
sh4 = lambda x, mu, k=1: tuple((np.array(x) + k * E[mu]) % S4)
Uf = {(x, mu): I2.copy() for x in itertools.product(range(S4), repeat=d) for mu in range(d)}
for nu, sgn in [(2, 1), (3, 1), (1, -1), (2, -1), (3, -1)]:
    Uf[(sh4(x0, nu, sgn), 0)] = -I2
def pl4(U, x, mu, nu):
    return U[(x, mu)] @ U[(sh4(x, mu), nu)] @ inv(U[(sh4(x, nu), mu)]) @ inv(U[(x, nu)])
good = True
for rep in range(5):
    Ue = haar(); Uf[(x0, 0)] = Ue
    tot = 0.0
    for nu in range(1, d):
        tot += pl4(Uf, x0, 0, nu).trace().real / 2 + pl4(Uf, sh4(x0, nu, -1), 0, nu).trace().real / 2
    good &= abs(tot - (-2 * Ue.trace().real)) < 1e-10 and np.allclose(pl4(Uf, x0, 0, 1), Ue)
check('(I) frustrated exterior: the six plaquettes of e add to -2 Tr U_e and p_0 has holonomy U_e', good)
fr = []
for beta in [0.25, 0.5, 1, 2, 10]:
    t0 = math.acos(1 - 2.5 ** 2 / 4)
    num = quad(lambda t: np.exp(-4 * beta * (np.cos(t) + 1)) * np.sin(t) ** 2, t0, math.pi)
    den = quad(lambda t: np.exp(-4 * beta * (np.cos(t) + 1)) * np.sin(t) ** 2, 0, math.pi)
    fr.append(num / den)
print('    frustrated exterior, rho = 2.5: P(p_0 bad) = %s at beta = 0.25, 0.5, 1, 2, 10' % ', '.join('%.6f' % v for v in fr))
check('(I) under the frustrated exterior p_0 is bad with probability increasing to 1',
      all(fr[i] <= fr[i + 1] + 1e-12 for i in range(4)) and fr[-1] > 0.999)
# chi = 1 off the bad cubes
rho_b = 1 / 12; good, nact, nfree, nnear = True, 0, 0, 0
for rep in range(4):
    g = {x: haar() for x in itertools.product(range(Sz), repeat=d)}
    centres = [tuple(rng.integers(0, Sz, size=d)) for _ in range(2)]
    U = {}
    for x in itertools.product(range(Sz), repeat=d):
        for mu in range(d):
            big = any(max(min((x[n] - c[n]) % Sz, (c[n] - x[n]) % Sz) for n in range(d)) == 0 for c in centres)
            a = rng.normal(size=3) * (0.8 if big else 1.5e-3)
            U[(x, mu)] = g[x] @ su2(a) @ inv(g[shz(x, mu, Sz)])
    bad = set()
    for x in itertools.product(range(Sz), repeat=d):
        for mu in range(d):
            for nu in range(mu + 1, d):
                if nrm(plaq(U, x, mu, nu, Sz) - I2) > rho_b:
                    for w in itertools.product(range(2), repeat=d):
                        if w[mu] == 0 and w[nu] == 0:
                            bad.add(tuple((np.array(x) - np.array(w)) % Sz))
    Ug = gfix(U, Sz)
    for y in itertools.product(range(0, Sz, 2), repeat=d):
        for mu in range(d):
            Z = {tuple((np.array(y) + t * E[mu] + np.array(w)) % Sz) for t in range(3)
                 for w in itertools.product(range(2), repeat=d) if w[mu] == 0}
            chi = sharpz(Ug, y, mu, Sz)[1]
            if Z & bad: nnear += 1; nact += chi < 1
            else: nfree += 1; good &= chi == 1.0
print('    coarse links with a bad cube in Z: %d (chi < 1 at %d); without: %d' % (nnear, nact, nfree))
check('(I) chi = 1 at every coarse link whose 24 cubes are good (rho = 1/12)', good and nfree > 0 and nnear > 0)
# (J)
def plaq_pts(x, mu, nu, Sz):
    return [tuple((np.array(x) + a * E[mu] + b * E[nu]) % Sz) for a in (0, 1) for b in (0, 1)]
def linf(p, q, Sz):
    return max(min((p[n] - q[n]) % Sz, (q[n] - p[n]) % Sz) for n in range(d))
def propagate(Bset, cs):
    out = set()
    for y in itertools.product(range(cs // 2), repeat=d):
        box = {tuple((2 * np.array(y) + np.array(o)) % cs) for o in itertools.product(range(4), repeat=d)}
        for mu in range(d):
            for nu in range(mu + 1, d):
                if any(x in box for (x, m, n) in Bset): out.add((y, mu, nu))
    return out
good_def, good_dist, good_E, rows = True, True, True, []
nonempty = True
for rep in range(3):
    g = {x: haar() for x in itertools.product(range(Sz), repeat=d)}
    c0 = tuple(rng.integers(0, Sz, size=d))
    U = {}
    for x in itertools.product(range(Sz), repeat=d):
        for mu in range(d):
            big = x == c0 and mu == 0
            a = rng.normal(size=3) * (0.9 if big else 1e-5)
            U[(x, mu)] = g[x] @ su2(a) @ inv(g[shz(x, mu, Sz)])
    eta0 = 4e-4
    L0 = logs(U, Sz)
    B = {k for k, Z in L0.items() if nrm(plaq(U, k[0], k[1], k[2], Sz) - I2) > 1 or nrm(Z) > eta0}
    B0pts = [pt for (x, m, n) in B for pt in plaq_pts(x, m, n, Sz)]
    Eo0 = sum(nrm(Z) ** 2 for k, Z in L0.items() if k not in B)
    cur, cs, eta, sgm, Bj = U, Sz, eta0, 0.0, B
    for j in range(1):
        nxt, cmin = step(cur, cs); ns = cs // 2
        Bn = propagate(Bj, cs); sgm += eta; eta = 4 * (1 + 647 * eta) * eta
        Ln = {}
        for (y, mu, nu) in itertools.product(itertools.product(range(ns), repeat=d), range(d), range(d)):
            if mu < nu and (y, mu, nu) not in Bn:
                P = plaq(nxt, y, mu, nu, ns)
                good_def &= nrm(P - I2) <= 1 and nrm(slog(P)) <= eta
                Ln[(y, mu, nu)] = slog(P)
        scale = 2 ** (j + 1)
        for (y, mu, nu) in Bn:
            for pt in plaq_pts(y, mu, nu, ns):
                fine = tuple(scale * q for q in pt)
                good_dist &= min(linf(fine, b, Sz) for b in B0pts) <= 3 * (scale - 1)
        Eo = sum(nrm(Z) ** 2 for Z in Ln.values())
        good_E &= Eo <= math.exp(12060 * sgm) * Eo0
        nonempty &= len(Ln) > 0 and len(B) > 0
        rows.append((rep, j + 1, len(B), len(Bn), max(nrm(Z) for Z in Ln.values()) / eta if Ln else 0.0, Eo / Eo0))
        cur, cs, Bj = nxt, ns, Bn
for r in rows:
    print('    field %d, step %d: |B_0| = %d, |B_j| = %d, max |Theta| off B_j / eta_j = %.3f, E_j off B_j / E_0 off B_0 = %.4f' % r)
check('(J) off the propagated defects |Theta^(j)| <= eta_j and delta <= 1', good_def and nonempty)
check('(J) the defects of the j-th iterate lie within 3(2^j - 1) fine spacings of B_0', good_dist)
check('(J) the energy off the defects obeys E_j <= e^{12060 sigma} E_0', good_E)
# (K) uniform upper bound for the effective forms: Whitney extension with nodes at the block centres
def whitney_int(u, L, dim):
    # returns Pi u times S = 16 L^4 (dim = 4) as an integer array; u integer, periodic coarse lattice
    n = u.shape[0]; F = n * L
    X = np.indices((F,) * dim)
    S = 2 * (2 * L) ** (dim - 1) * L
    phi = np.zeros((F,) * dim + (dim,), dtype=np.int64)
    for mu in range(dim):
        acc = np.zeros((F,) * dim, dtype=np.int64)
        for z in itertools.product(range(n), repeat=dim):
            w = np.ones((F,) * dim, dtype=np.int64)
            for nu in range(dim):
                c2 = 2 * L * z[nu] + L - 1          # twice the centre
                x2 = 2 * X[nu]
                tot = np.zeros((F,) * dim, dtype=np.int64)
                for sh in (-2 * F, 0, 2 * F):
                    if nu == mu:   # twice the overlap of [x, x+1) with [c, c+L)
                        tot += np.maximum(0, np.minimum(x2 + 2, c2 + sh + 2 * L) - np.maximum(x2, c2 + sh))
                    else:          # 2L times the hat function
                        tot += np.maximum(0, 2 * L - np.abs(x2 - c2 - sh))
                w = w * tot
            acc += u[z + (mu,)] * w
        phi[..., mu] = acc
    return phi, S              # phi = S (Pi u): overlap doubled, hats times 2L, and the 1/L of omega
def curv_int(f, dim):
    out = []
    for mu in range(dim):
        for nu in range(mu + 1, dim):
            out.append(f[..., mu] + np.roll(f[..., nu], -1, mu) - np.roll(f[..., mu], -1, nu) - f[..., nu])
    return np.array(out)
rngK = np.random.default_rng(7)
good_T, good_Q, good_g = True, True, True
for (dim, n, L) in ((4, 3, 2), (4, 3, 4), (2, 4, 8), (3, 3, 4)):
    for trial in range(2):
        u = rngK.integers(-3, 4, size=(n,) * dim + (dim,))
        phi, S = whitney_int(u, L, dim)
        F = n * L
        # L^4 S times the k-fold averages, from straight paths of length L
        lhs = np.zeros((n,) * dim + (dim,), dtype=np.int64)
        for mu in range(dim):
            ps = sum(np.roll(phi[..., mu], -t, mu) for t in range(L))
            for z in itertools.product(range(n), repeat=dim):
                sl = tuple(slice(L * zi, L * zi + L) for zi in z)
                lhs[z + (mu,)] = ps[sl].sum()
        Tu = u.astype(np.int64) * 1
        for nu in range(dim):
            Tu = 6 * Tu + np.roll(Tu, 1, nu) + np.roll(Tu, -1, nu)     # 8 tau
        good_T &= np.array_equal(lhs * 8 ** dim, Tu * S * L ** dim)
        good_Q &= int((curv_int(phi, dim) ** 2).sum()) <= int((curv_int(u, dim) ** 2).sum()) * S * S
        b = rngK.integers(-3, 4, size=(n,) * dim)
        db = np.stack([np.roll(b, -1, m) - b for m in range(dim)], -1)
        good_g &= not np.any(curv_int(whitney_int(db, L, dim)[0], dim))
check('(K) the k-fold averages of the Whitney extension are T u, T = tau^(x d), tau = (1/8, 3/4, 1/8) (exact)', good_T)
check('(K) Q(Pi u) <= Q_k(u) and Pi(db) has zero curvature (exact)', good_Q and good_g)
th = np.linspace(0, np.pi, 2001)
check('(K) the symbol 3/4 + cos(theta)/4 of tau is at least 1/2, and 1/L with coarse sites as nodes',
      np.min(0.75 + 0.25 * np.cos(th)) >= 0.5 - 1e-15 and all(abs((L + 1) / (2 * L) - (L - 1) / (2 * L) - 1 / L) < 1e-15 for L in (2, 4, 8, 16)))
print('%d checks passed, %d failed' % (sum(ok), len(ok) - sum(ok)))
