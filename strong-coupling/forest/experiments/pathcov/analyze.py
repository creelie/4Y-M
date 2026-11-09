"""Covariance sums of path scores from pathcov .bin output (V = identity, translation invariant)."""
import json, sys, glob, numpy as np
from scipy.special import iv

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

def summarize(C, L):
    idx = np.indices((L,)*4).reshape(4, -1).T
    d = np.minimum(idx, L-idx).max(axis=1)
    nrm = np.linalg.norm(C.reshape(-1, 4, 4, 3, 3), ord=2, axis=(-2, -1))  # (L^4, mu, nu)
    var = nrm[0].diagonal().mean()
    out = {'var': float(var)}
    for nu in range(4):
        col = nrm[:, :, nu].copy(); col[0, nu] = 0  # e != j
        sh = [float(col[d == r].sum()) for r in range(L//2+1)]
        out.setdefault('shells', []).append(sh)
    sh = np.array(out['shells']).mean(axis=0)
    out['shell_mean'] = sh.tolist(); out['total'] = float(sh.sum())
    out['same_base'] = float(np.mean([nrm[0, mu, nu] for nu in range(4) for mu in range(4) if mu != nu]) * 3)
    return out

rows = []
for lg in sorted(glob.glob(sys.argv[1] + '/*.log')):
    meta = json.loads(open(lg).readline()); L = meta['L']
    A = np.fromfile(lg[:-4] + '.bin', dtype=np.float32).astype(np.float64)
    T = A.size // (L**4*4*3); A = A[:T*L**4*12].reshape(T, L, L, L, L, 4, 3)
    # coarse layout: path e = y*4+mu, y = c0 + L c1 + ... -> reshape gives (c3,c2,c1,c0); fine for translation sums
    full = summarize(covfield(A), L)
    h1 = summarize(covfield(A[:T//2]), L); h2 = summarize(covfield(A[T//2:]), L)
    # noise floor: same estimator on time-permuted pairing (true covariance zero)
    null = summarize(covfield(A, A[np.random.default_rng(0).permutation(T)]), L)
    meta.update(T=T, **{k: full[k] for k in ('var', 'total', 'same_base', 'shell_mean')},
                total_halves=[h1['total'], h2['total']], null_total=null['total'], null_shell=null['shell_mean'], ratio=full['total']/full['var'], **consts(meta['beta'], meta['b']))
    rows.append(meta); print(json.dumps(meta), flush=True)
json.dump(rows, open(sys.argv[1] + '/summary.json', 'w'), indent=1)
