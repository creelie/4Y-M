"""k -> 0 limit of the transverse massive (U) test (A_GAUSSIAN_FIBRE_4D §10). At k -> 0, S -> mu^2 I on the transverse block and
Qhat(0) = q I (cubic symmetry), q = e^T (K_in + mu^2)^{-1} e, e = indicator of interior 0-direction links. Hence
varpi_lim(n) = mu^2 n^-4 q, mu = m ell / n. Sparse K_in, CG solve, so n can go well beyond 8."""
import sys, json, itertools, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla

def sparse_K(n):
    # links (x,mu) with x in [0,n)^4, x_mu <= n-2 are interior; plaquettes with base in [-1,n-1]^4
    def lid(x, d):
        x = np.asarray(x)
        ok = np.all((x >= 0) & (x < n), axis=-1) & (x[..., d] <= n-2)
        idx = d*n**4 + np.ravel_multi_index(np.clip(x, 0, n-1).T, (n,)*4)
        return np.where(ok, idx, -1)
    B = np.array(list(itertools.product(range(-1, n), repeat=4)))
    rows, cols, vals = [], [], []; r0 = 0
    for a in range(4):
        for c in range(a+1, 4):
            ea, ec = np.eye(4, dtype=int)[a], np.eye(4, dtype=int)[c]
            for p, d, s in ((B, a, 1), (B+ea, c, 1), (B+ec, a, -1), (B, c, -1)):
                l = lid(p, d); m = l >= 0
                rows.append(r0+np.where(m)[0]); cols.append(l[m]); vals.append(np.full(m.sum(), s))
            r0 += len(B)
    d1 = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(r0, 4*n**4))
    X = np.array(list(itertools.product(range(n), repeat=4)))
    live = np.concatenate([lid(X, d) for d in range(4)]); live = np.sort(live[live >= 0])
    d1 = d1[:, live]; K = (d1.T @ d1).tocsr()
    e = (live < n**4).astype(float)   # direction-0 links
    return K, e, live

def run(n, mells):
    K, e, _ = sparse_K(n)
    for ml in mells:
        mu = ml/n; A = K + mu**2*sp.identity(K.shape[0], format='csr')
        x, info = sla.cg(A, e, rtol=1e-10, maxiter=200000)
        q = float(e @ x); print(json.dumps(dict(n=n, m_ell=ml, q_over_n6=q/n**6, varpi_lim=mu**2*q/n**4, cg_info=int(info))), flush=True)

def selftest():
    from gauss_rho import cube_K
    for n in (2, 3):
        K, e, live = sparse_K(n); inl, Kd = cube_K(1, n+4, n, (0,)*4)
        # map dense ordering to sparse ordering
        idx = {d*n**4 + np.ravel_multi_index(x, (n,)*4): i for (x, d), i in inl.items()}
        perm = np.array([idx[l] for l in live]); err = np.abs(K.toarray() - Kd[np.ix_(perm, perm)]).max()
        print('n', n, 'K err', err); assert err < 1e-12
    print('SELFTEST PASS')

if __name__ == '__main__':
    if sys.argv[1] == 'selftest': selftest()
    else: run(int(sys.argv[1]), [float(v) for v in sys.argv[2:]])
