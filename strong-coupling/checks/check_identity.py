"""Numerical check of the exact identity (2.15) and inequality (2.8) of
OpenAI, 'Exponential decay in two-dimensional classical O(n) models' (2026-09-23),
on small graphs with prescribed spins, by quadrature over (S^2)^k.

Identity checked (with T' = sum_i u_i conj(u_i)^T, the Hermitian contraction):
   R(Q, conj Q) = sum_{i,j} R(S_ij, conj S_ij) - sum_{i,j} D^4W(u_i a, conj u_i b, conj u_j a, u_j b)
where R(h,k) = -D^2 W(hc,kc), S_ij = S(f_i,f_j), Q = sum_i S(f_i, conj f_i).
Numerical evidence only (quadrature), used to confirm signs and conjugations.
"""
import numpy as np, itertools
rng = np.random.default_rng(1)
E1, E2, E3 = np.eye(3)
def crossmat(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
A, B, C = crossmat(E1), crossmat(E2), crossmat(E3)
assert np.allclose(A @ B - B @ A, C)
GEN = {'a': A, 'b': B, 'c': C}

def sphere_nodes(nu, nphi):
    u, w = np.polynomial.legendre.leggauss(nu)
    phi = 2 * np.pi * (np.arange(nphi) + 0.5) / nphi
    U, P = np.meshgrid(u, phi, indexing='ij')
    s = np.sqrt(1 - U ** 2)
    pts = np.stack([s * np.cos(P), s * np.sin(P), U], -1).reshape(-1, 3)
    wts = (np.repeat(w, nphi) / 2 / nphi)
    return pts, wts

def sym_prod(gs):
    k = len(gs)
    M = np.zeros((3, 3))
    for perm in itertools.permutations(range(k)):
        P = np.eye(3)
        for p in perm:
            P = P @ GEN[gs[p]]
        M += P
    return M / __import__("math").factorial(k)

def run(nu=24, nphi=32, nint=2, extra_edges=True, beta=0.9, nfam=3, seed=0):
    rng = np.random.default_rng(seed)
    # vertices: 0..nint-1 integrated, then two pinned
    npin = 2
    V = nint + npin
    pinned = {nint + k: (lambda v: v / np.linalg.norm(v))(rng.normal(size=3)) for k in range(npin)}
    edges = [(nint, 0)] + [(i, i + 1) for i in range(nint - 1)] + [(nint - 1, nint + 1)]
    if extra_edges:
        edges += [(nint, nint - 1), (0, nint + 1)]
    b = rng.uniform(0.2, beta, size=len(edges))
    pts, wts = sphere_nodes(nu, nphi)
    # configurations: product grid over integrated vertices
    grids = np.meshgrid(*[np.arange(len(wts))] * nint, indexing='ij')
    idx = [g.ravel() for g in grids]
    W = np.prod([wts[i] for i in idx], axis=0)
    S = {}
    for v in range(V):
        S[v] = pts[idx[v]] if v < nint else np.broadcast_to(pinned[v], (len(W), 3))
    energy = sum(be * np.einsum('ni,ni->n', S[x], S[y]) for be, (x, y) in zip(b, edges))
    gibbs = W * np.exp(energy)
    gibbs /= gibbs.sum()
    def Je(gs):  # array (nconf, nedges)
        M = sym_prod(gs)
        return np.stack([be * np.einsum('ni,ij,nj->n', S[x], M, S[y]) for be, (x, y) in zip(b, edges)], -1)
    cache = {}
    def J(gs):
        key = tuple(gs)
        if key not in cache: cache[key] = Je(gs)
        return cache[key]
    def Hk(args):  # args list of (v, g) with v edge form (complex array over edges)
        gs = [g for _, g in args]
        prod = np.ones(len(edges), complex)
        for v, _ in args: prod = prod * v
        return J(gs) @ prod
    def X(v, g): return Hk([(v, g)])
    def E(z): return np.sum(gibbs * z)
    def D2W(h, k, g='c'): return E(Hk([(h, g), (k, g)]) + X(h, g) * X(k, g))
    def R(h, k): return -D2W(h, k)
    def D4W(args):
        idxs = range(4); tot = 0
        # all set partitions of {0,1,2,3}
        def partitions(s):
            if not s: yield []; return
            first, rest = s[0], s[1:]
            for p in partitions(rest):
                for i in range(len(p)):
                    yield p[:i] + [[first] + p[i]] + p[i + 1:]
                yield [[first]] + p
        for part in partitions(list(idxs)):
            term = np.ones(len(gibbs), complex)
            for block in part:
                term = term * Hk([args[r] for r in block])
            tot = tot + term
        return E(tot)
    # family of complex vertex functions vanishing on pinned vertices
    fam = []
    for _ in range(nfam):
        f = np.zeros(V, complex); f[:nint] = rng.normal(size=nint) + 1j * rng.normal(size=nint)
        fam.append(f)
    def d(f): return np.array([f[y] - f[x] for (x, y) in edges])
    def Sform(f, g): return np.array([0.5 * (f[x] * g[y] - g[x] * f[y]) for (x, y) in edges])
    us = [d(f) for f in fam]
    Q = sum(Sform(f, f.conj()) for f in fam)
    lhs = R(Q, Q.conj())
    cross = sum(R(Sform(fi, fj), Sform(fi, fj).conj()) for fi in fam for fj in fam)
    d4 = sum(D4W([(ui, 'a'), (ui.conj(), 'b'), (uj.conj(), 'a'), (uj, 'b')]) for ui in us for uj in us)
    rhs = cross - d4
    # inequality (2.8) with Hermitian T'
    Tp = sum(np.outer(u, u.conj()) for u in us)
    T = sum(np.outer(u, u) for u in us)
    Krot = 6 * beta ** 2 + 4 * beta ** 1.5 + beta
    def mn(fi, fj):
        return min(np.max(abs(fi)) ** 2 * np.sum(abs(d(fj)) ** 2), np.max(abs(fj)) ** 2 * np.sum(abs(d(fi)) ** 2))
    bound = Krot * np.sum(abs(Tp) ** 2) + beta * sum(mn(fi, fj) for fi in fam for fj in fam)
    # gradient invariance check R(h, d f) = 0
    h = rng.normal(size=len(edges)); gi = R(h, d(fam[0]))
    return dict(lhs=lhs, rhs=rhs, err=abs(lhs - rhs), bound=bound, d4=d4, cross=cross,
                gradinv=gi, HS_T=np.sum(abs(T) ** 2), HS_Tp=np.sum(abs(Tp) ** 2))

if __name__ == '__main__':
    for seed in range(4):
        for nint in (2, 3):
            nu, nphi = (20, 24) if nint == 2 else (8, 10)
            r = run(nu=nu, nphi=nphi, nint=nint, seed=seed)
            print(f"seed={seed} nint={nint}  R(Q,Qbar)={r['lhs'].real:+.6e}{r['lhs'].imag:+.1e}i  "
                  f"rhs={r['rhs'].real:+.6e}{r['rhs'].imag:+.1e}i  |diff|={r['err']:.1e}  "
                  f"bound(2.8)={r['bound']:.4e}  R(h,df)={abs(r['gradinv']):.1e}")
