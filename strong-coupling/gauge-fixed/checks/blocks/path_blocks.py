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
print('%d checks passed, %d failed' % (sum(ok), len(ok) - sum(ok)))
