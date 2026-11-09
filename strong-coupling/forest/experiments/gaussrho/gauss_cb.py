"""First-chaos lower bound on the optimal split constant C_B (A_GAUSSIAN_FIBRE_4D §8.3), b = 2, n' = 2n.
Test functions f = w.c, linear in the coarse variables c = Pth^T A (so Var(f|c) = 0, F_f = f), with w
divergence-free and zero-mean per direction. For these, S(2n,n;C_A,C_B) forces C_B >= Q_c(w)/Q_f(w) for every C_A:
  Q_c(w) = n^-4 sum_y Var(w.c | coarse links outside cube y)  = n^-4 sum_y w_in^T ((M^+)_in,in)^+ w_in
  Q_f(w) = (2n)^-4 sum_x Var(f | fine links outside cube x)  = (2n)^-4 sum_x (Pth w)_in^T (K_in)^+ (Pth w)_in
C_B^lin = top generalized eigenvalue of (Q_c, Q_f) on that subspace."""
import sys, json, itertools, numpy as np, scipy.linalg as sl
from gauss_rho import sites, sidx, d0T, coarse_M, cube_K
from gauss_angle import nullsp, whiten

def scatter(Q, block, Y, mus, L, scale):
    for z in sites(L):
        idx = sidx(Y + z, L)*4 + mus; Q[np.ix_(idx, idx)] += scale*block

def coarse_cond_block(Mp, L, n):  # precision-based conditional covariance of in-links, cube at 0
    inl, _ = cube_K(1, L, n, (0,)*4)
    Y = np.array([x for (x, mu) in inl]); mus = np.array([mu for (x, mu) in inl])
    gin = sidx(Y, L)*4 + mus
    return Y, mus, gin, np.linalg.pinv(Mp[np.ix_(gin, gin)], hermitian=True)

def build(L, n):
    b, npr = 2, 2*n; N = b*L; nc = 4*L**4
    M = coarse_M(b, L); Mp = np.linalg.pinv(M, hermitian=True)
    Qc = np.zeros((nc, nc)); Y, mus, gin, Bc = coarse_cond_block(Mp, L, n)
    scatter(Qc, Bc, Y, mus, L, n**-4.0)
    # fine: K_in identical for every fine cube; the coarse columns it touches depend on x mod b
    inl, Kin = cube_K(b, L, npr, (0,)*4); Kp = np.linalg.pinv(Kin, hermitian=True)
    Qf = np.zeros((nc, nc))
    for r in itertools.product(range(b), repeat=4):
        inl, _ = cube_K(b, L, npr, r)
        Bx = np.zeros((len(inl), nc))
        for (x, mu), i in inl.items():
            x = np.array(x)  # straight path of coarse (y,mu): fine links (b y + i e_mu, mu), i < b
            if all(x[nu] % b == 0 for nu in range(4) if nu != mu): Bx[i, sidx(x // b, L)*4 + mu] = 1.0
        cols = np.where(np.abs(Bx).sum(0) > 0)[0]; blk = Bx[:, cols].T @ Kp @ Bx[:, cols]
        Yr = np.array([[c//4 // L**(3-k) % L for k in range(4)] for c in cols]); mr = cols % 4
        # sidx convention check
        assert np.all(sidx(Yr, L)*4 + mr == cols)
        scatter(Qf, blk, Yr, mr, L, (2.0*n)**-4.0)
    D = d0T(L); Z = nullsp(D); Cst = np.kron(np.ones((L**4, 1)), np.eye(4))
    W = Z @ nullsp(Cst.T @ Z)
    return M, Mp, Qc, Qf, W

def run(L, n):
    M, Mp, Qc, Qf, W = build(L, n)
    A = W.T @ Qc @ W; B = W.T @ Qf @ W
    wB = np.linalg.eigvalsh((B+B.T)/2); Wh = whiten(B)
    ev = np.linalg.eigvalsh(Wh.T @ A @ Wh)
    print(json.dumps(dict(L=L, n=n, dimW=W.shape[1], Bmin_over_max=float(wB.min()/wB.max()),
                          CB_lin=float(ev.max()), top=ev[-5:][::-1].tolist(), median=float(np.median(ev)))), flush=True)

def selftest():
    L, n = 4, 2; b = 2; M = coarse_M(b, L); Mp = np.linalg.pinv(M, hermitian=True)
    Y, mus, gin, Bc = coarse_cond_block(Mp, L, n)
    D = d0T(L); Z = nullsp(D); rng = np.random.default_rng(0); w = Z @ rng.standard_normal((Z.shape[1], 3))
    out = np.setdiff1d(np.arange(4*L**4), gin)
    MoW = M[out] @ w
    direct = np.einsum('ij,ij->j', w, M @ w) - np.einsum('ij,ij->j', MoW, np.linalg.pinv(M[np.ix_(out, out)], hermitian=True) @ MoW)
    viaP = np.einsum('ij,ij->j', w[gin], Bc @ w[gin])
    e = np.abs(direct - viaP).max()/np.abs(direct).max(); print('coarse cond-var rel err', e, direct.tolist()); assert e < 1e-6
    from gauss_rho import sites as _s
    X = _s(L); assert np.all(sidx(X, L) == np.arange(L**4)), 'sidx ordering'
    print('SELFTEST PASS')

if __name__ == '__main__':
    if sys.argv[1] == 'selftest': selftest()
    else: run(int(sys.argv[1]), int(sys.argv[2]))
