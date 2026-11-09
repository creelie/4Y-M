"""Gaussian (abelian, beta-free) exterior constant rho_1 of Thm 114 (X), first chaos.
Model: fine torus N=bL, field A on links, energy |dA|^2/2. Gauge-invariant linear functionals = divergence-free currents u,
Var(u) = u^T K^+ u = u^T (-Lap)^+ u (Hodge identity on div-free currents; k=0 harmonic mode dropped).
H_c: div-free combinations of straight path currents P g (coarse div-free g). Cube V' = n'^4 sites at offset o;
'in' links have both ends in V'. F = gauge-invariant data outside: Var(u|F) = u_in^T K_in^+ u_in.
T = coarse links whose path has an in-link; H_c,out = div-free g supported off T (F-measurable).
C' = H_c ⊖ H_c,out (in the Var metric). rho_1^2 = 1 - min_{u in C'} Var(u|F)/Var(u)."""
import sys, json, itertools, numpy as np, scipy.linalg as sl

def sites(L): return np.array(list(itertools.product(range(L), repeat=4)))
def sidx(x, L): x = np.mod(x, L); return ((x[..., 0]*L + x[..., 1])*L + x[..., 2])*L + x[..., 3]

def d0T(L):  # divergence: (L^4) x (4 L^4), column (y,mu): +1 at y+mu, -1 at y
    n = L**4; D = np.zeros((n, 4*n)); X = sites(L)
    for mu in range(4):
        e = np.eye(4, dtype=int)[mu]; cols = sidx(X, L)*4 + mu
        D[sidx(X+e, L), cols] += 1; D[sidx(X, L), cols] -= 1
    return D

def scalar_green(N):
    k = 2*np.pi*np.fft.fftfreq(N); lam = sum(np.meshgrid(*(2-2*np.cos(k),)*4, indexing='ij'))
    with np.errstate(divide='ignore'): ig = np.where(lam > 1e-12, 1/lam, 0.0)
    return np.fft.ifftn(ig).real

def coarse_M(b, L):  # M[(y,mu),(y',nu)] = Var-inner product of path currents
    N = b*L; g = scalar_green(N); n = L**4; X = sites(L); M = np.zeros((4*n, 4*n))
    for mu in range(4):
        e = np.eye(4, dtype=int)[mu]
        # kernel over coarse displacement z: sum_{i,j} g(b z + (i-j) e)
        Z = sites(L)
        ker = sum(g[tuple(np.mod(b*Z + (i-j)*e, N).T)] for i in range(b) for j in range(b))
        yi = sidx(X, L)
        for a in range(n):
            M[a*4+mu, yi*4+mu] = ker[sidx(X - X[a], L)]
    return M

def cube_K(b, L, npr, off):
    N = b*L; lo = np.array(off); X = lo + sites(npr)
    inl = {}
    for x in X:
        for mu in range(4):
            y = x.copy(); y[mu] += 1
            if np.all(y < lo + npr): inl[(tuple(np.mod(x, N)), mu)] = len(inl)
    m = len(inl); K = np.zeros((m, m))
    seen = set()
    for (x, mu) in list(inl):
        x = np.array(x)
        for nu in range(4):
            if nu == mu: continue
            for s in (0, 1):  # plaquettes (x,mu,nu) and (x-nu,mu,nu)
                base = x.copy();
                if s: base[nu] -= 1
                key = (tuple(np.mod(base, N)), min(mu, nu), max(mu, nu))
                if key in seen: continue
                seen.add(key); a, c = min(mu, nu), max(mu, nu)
                ea, ec = np.eye(4, dtype=int)[a], np.eye(4, dtype=int)[c]
                links = [(base, a, 1), (base+ea, c, 1), (base+ec, a, -1), (base, c, -1)]
                ids = [(inl.get((tuple(np.mod(p, N)), d)), sg) for p, d, sg in links]
                for i1, s1 in ids:
                    for i2, s2 in ids:
                        if i1 is not None and i2 is not None: K[i1, i2] += s1*s2
    return inl, K

def rho1(b, L, npr, off):
    n = L**4; N = b*L
    inl, K = cube_K(b, L, npr, off)
    Pin = np.zeros((len(inl), 4*n)); X = sites(L)
    for a, y in enumerate(X):
        for mu in range(4):
            for i in range(b):
                f = b*y; f[mu] += i; r = inl.get((tuple(np.mod(f, N)), mu))
                if r is not None: Pin[r, a*4+mu] = 1
    T = np.where(np.abs(Pin).sum(0) > 0)[0]; O = np.setdiff1d(np.arange(4*n), T)
    W = Pin.T @ np.linalg.pinv(K, hermitian=True) @ Pin
    M = coarse_M(b, L); D = d0T(L)
    Z = sl.null_space(D); ZO = np.zeros((4*n, 0))
    zo = sl.null_space(D[:, O]); ZO = np.zeros((4*n, zo.shape[1])); ZO[O] = zo
    Bc = sl.null_space(ZO.T @ M @ Z); U = Z @ Bc
    Mc = U.T @ M @ U; Wc = U.T @ W @ U
    ev = sl.eigh(Wc, Mc, eigvals_only=True)
    return dict(b=b, L=L, nprime=npr, off=list(off), nT=int(len(T)), dimC=int(U.shape[1]), n_in=len(inl),
                min_ratio=float(ev.min()), max_ratio=float(ev.max()), rho1=float(np.sqrt(max(0.0, 1-ev.min()))))

def selftest():
    # (1) Hodge identity: u^T K^+ u == u^T (-Lap)^+ u for div-free u on fine torus (dense, N=4)
    N = 4; n = N**4; D = d0T(N); Kfull = np.zeros((4*n, 4*n))
    # K = d1^T d1 built from cube_K on the whole torus (cube = whole lattice wraps -> build directly)
    X = sites(N); rows = []
    for x in X:
        for a in range(4):
            for c in range(a+1, 4):
                r = np.zeros(4*n); ea, ec = np.eye(4, dtype=int)[a], np.eye(4, dtype=int)[c]
                r[sidx(x, N)*4+a] += 1; r[sidx(x+ea, N)*4+c] += 1; r[sidx(x+ec, N)*4+a] -= 1; r[sidx(x, N)*4+c] -= 1
                rows.append(r)
    d1 = np.array(rows); Kfull = d1.T @ d1
    u = sl.null_space(D) @ np.random.default_rng(0).standard_normal(sl.null_space(D).shape[1])
    u -= np.repeat(u.reshape(n, 4).mean(0)[None], n, 0).ravel()  # drop harmonic part
    g = scalar_green(N); G = np.zeros((4*n, 4*n))
    for a in range(n):
        for mu in range(4): G[a*4+mu, sidx(X, N)*4+mu] = g[tuple(np.mod(X - X[a], N).T)]
    v1 = u @ np.linalg.pinv(Kfull, hermitian=True) @ u; v2 = u @ G @ u
    print(f"hodge: K+ {v1:.10f}  G {v2:.10f}  rel {abs(v1-v2)/v1:.1e}")
    assert abs(v1-v2)/v1 < 1e-8
    # (2) cube_K on a cube must agree with Kfull block
    inl, K = cube_K(1, N, 3, (0, 0, 0, 0)); idx = [sidx(np.array(x), N)*4+mu for (x, mu) in inl]
    err = np.abs(K - Kfull[np.ix_(idx, idx)]).max(); print("cubeK block err", err); assert err < 1e-12
    # (3) small rho run: ratio in [0,1]
    r = rho1(3, 4, 3, (0, 0, 0, 0)); print(r); assert -1e-9 <= r['min_ratio'] <= 1+1e-9
    print("SELFTEST PASS")

if __name__ == '__main__':
    if sys.argv[1] == 'selftest': selftest()
    else:
        b, L, npr = map(int, sys.argv[1:4]); off = tuple(map(int, sys.argv[4].split(',')))
        print(json.dumps(rho1(b, L, npr, off)), flush=True)
