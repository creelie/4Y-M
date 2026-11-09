"""Debiased per-shell path-score covariance (run 2).
Split the chain into K contiguous blocks. For every pair (e,j) and entry, the cross-block product
  Q = (S^2 - sum_k C_k^2) / (K(K-1)),  S = sum_k C_k
is unbiased for Cov^2 when blocks are independent (blocks >> autocorrelation time).
Per shell r: Q_r = sum over pairs in shell of ||Cov||_F^2 (unbiased), jackknife error over blocks.
Reported: Q_r, z_r = Q_r/err_r, and the Cauchy-Schwarz bound sum_r sqrt(n_r Q_r) >= sum ||Cov||_F."""
import json, sys, glob, numpy as np

def consts(bW, b):
    BW = 6*bW; qW = 4*bW/(4*bW+1); m = b-2
    eta = 48*(1+BW)**2*qW**(m-2)
    h = (BW**2+BW)*eta/(1-eta)+BW**2*eta**2/(1-eta)**2 if eta < 1 else float('nan')
    K = BW*eta/(1-eta) if eta < 1 else float('nan')
    return dict(eta=eta, h=h, K=K, thr=(1-h)/3 if eta < 1 else float('nan'))

def covfield(A, B=None):  # A: (T, L,L,L,L, 4, 3) -> C[dx(4), mu, nu, a, b]
    A = A - A.mean(axis=(0, 1, 2, 3, 4), keepdims=True)
    B = A if B is None else B - B.mean(axis=(0, 1, 2, 3, 4), keepdims=True)
    F = np.fft.fftn(B, axes=(1, 2, 3, 4)); F0 = np.fft.fftn(A, axes=(1, 2, 3, 4))
    # C(D) = mean_y <a_y a_{y+D}>
    G = np.einsum('t...ma,t...nb->...mnab', F0.conj(), F) / A.shape[0]
    C = np.fft.ifftn(G, axes=(0, 1, 2, 3)).real / np.prod(A.shape[1:5])
    return C


def shellidx(L):
    idx = np.indices((L,)*4).reshape(4, -1).T
    return np.minimum(idx, L-idx).max(axis=1)

def blockcovs(A, K):
    T = A.shape[0] // K
    A = A - A.mean(axis=(0, 1, 2, 3, 4), keepdims=True)   # global mean; O(1/T) bias, negligible vs block noise
    return [covfield(A[k*T:(k+1)*T]) for k in range(K)]   # each: (L,L,L,L,4,4,3,3)

def shellQ(Cs, L):
    K = len(Cs); d = shellidx(L)
    Cs = np.stack([C.reshape(L**4, 4, 4, 3, 3) for C in Cs])
    def Q(Cs):
        k = len(Cs); S = Cs.sum(0)
        q = ((S**2).sum((-1, -2)) - (Cs**2).sum(0).sum((-1, -2))) / (k*(k-1))   # (L^4, mu, nu)
        for nu in range(4): q[0, nu, nu] = 0     # e != j
        return np.array([q[d == r].sum() / 4 for r in range(L//2+1)])  # per reference link j (avg over nu)
    full = Q(Cs)
    jk = np.array([Q(np.delete(Cs, i, 0)) for i in range(K)])
    err = np.sqrt((K-1)/K*((jk - jk.mean(0))**2).sum(0))
    npair = np.array([(d == r).sum()*4 for r in range(L//2+1)], float); npair[0] -= 1
    var = float(np.mean([np.linalg.norm(Cs.mean(0)[0, m, m], 2) for m in range(4)]))
    return full, err, npair, var

def report(A, L, K, meta):
    Q, e, n, var = shellQ(blockcovs(A, K), L)
    meta.update(K=K, var=var, shellQ=Q.tolist(), shellQerr=e.tolist(), shellZ=(Q/np.maximum(e, 1e-300)).tolist(),
                CSbound=float(np.sqrt(np.clip(Q, 0, None)*n).sum()),
                CSbound_sig=float(np.sqrt(np.clip(Q, 0, None)*n)[Q > 3*e].sum()))
    return meta

if __name__ == '__main__':
    K = 20; rows = []
    if sys.argv[1] == 'control':     # synthetic bphidden-type control: hidden global sign s mixes product laws
        rng = np.random.default_rng(1)
        for L in (4, 6, 8):
            for amp in (0.0, 0.1):
                T = 2000; s = rng.choice([-1., 1.], T)
                A = rng.standard_normal((T, L, L, L, L, 4, 3)) + amp*s[:, None, None, None, None, None, None]
                m = report(A, L, K, dict(control=True, L=L, amp=amp, T=T,
                                         ))
                m['truth_per_shell'] = (9*amp**4*np.array([(shellidx(L) == r).sum()*4 - (r == 0) for r in range(L//2+1)])).tolist()
                rows.append(m); print(json.dumps(m), flush=True)
    else:
        for lg in sorted(glob.glob(sys.argv[1] + '/*.log')):
            meta = json.loads(open(lg).readline()); L = meta['L']
            A = np.fromfile(lg[:-4] + '.bin', dtype=np.float32)
            T = A.size // (L**4*12); A = A[:T*L**4*12].reshape(T, L, L, L, L, 4, 3).astype(np.float64)
            m = report(A, L, K, dict(meta, T=T, **consts(meta['beta'], meta['b'])))
            rows.append(m); print(json.dumps(m), flush=True)
    json.dump(rows, open(sys.argv[2], 'w'), indent=1)
