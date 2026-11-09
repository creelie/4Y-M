"""Step 3: two-scale comparison in the free massive field. Same lattice (mu = m a fixed), cubes of side n and n0 = n/2.
For each k, the transverse BP symbols H_n(k) (side n) and H_n0(k) (side n0) are 2x2/3x3 Hermitian; report
C(n) = inf_k lambda_min(H_n0^{-1/2} H_n H_n0^{-1/2}), i.e. the best constant in E_n >= C E_n0. Uniform C in n at fixed
n/n0 = 2 (fixed ell/ell0) is the free-field version of the §8.4 two-scale comparison."""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from gauss_k0 import sparse_K
from gauss_massU import Khat

class Sym:
    def __init__(s, n, mu):
        s.n, s.mu = n, mu; K, e, live = sparse_K(n); s.A = (K + mu**2*sp.identity(K.shape[0], format='csr')).tocsr()
        s.d = live // n**4; x = np.array(np.unravel_index(live % n**4, (n,)*4)).T; s.R = x + 0.5*np.eye(4)[s.d]
    def H(s, k, B):
        ph = np.exp(1j*s.R @ k); E = np.zeros((len(ph), 4), complex); E[np.arange(len(ph)), s.d] = ph
        X = np.column_stack([sla.cg(s.A, E[:, j].real, rtol=1e-10)[0] + 1j*sla.cg(s.A, E[:, j].imag, rtol=1e-10)[0] for j in range(4)])
        Q = E.conj().T @ X / s.n**4; Q = (Q + Q.conj().T)/2
        S = Khat(k) + s.mu**2*np.eye(4); w, V = np.linalg.eigh(S); Sh = V @ np.diag(np.sqrt(w)) @ V.T
        return B.conj().T @ Sh @ Q @ Sh @ B

def basis(k):
    c = 2*np.sin(k/2)
    if c @ c < 1e-14: return np.eye(4)
    P = np.eye(4) - np.outer(c, c)/(c @ c); w, V = np.linalg.eigh(P); return V[:, w > 0.5]

def run(n, mell, nt):
    mu = mell/n; big, small = Sym(n, mu), Sym(n//2, mu)
    dirs = {'axis': [1, 0, 0, 0], 'diag2': [1, 1, 0, 0], 'diag4': [1, 1, 1, 1], 'gen': [1, 2, 3, 5]}
    best = (np.inf, None)
    for name, v in dirs.items():
        v = np.array(v, float); v /= np.linalg.norm(v)
        for t in [0.0] + list(np.geomspace(0.05, np.pi*n, nt)):
            k = t/n*v; B = basis(k); Hb, Hs = big.H(k, B), small.H(k, B)
            w, V = np.linalg.eigh(Hs); Hi = V @ np.diag(w**-0.5) @ V.conj().T
            r = float(np.linalg.eigvalsh(Hi @ Hb @ Hi).min())
            print(json.dumps(dict(n=n, n0=n//2, m_ell=mell, dir=name, t=t, ratio=r, varpi_n=float(np.linalg.eigvalsh(Hb).min()), varpi_n0=float(np.linalg.eigvalsh(Hs).min()))), flush=True)
            if r < best[0]: best = (r, [name, t])
    print(json.dumps(dict(summary=True, n=n, n0=n//2, m_ell=mell, C=best[0], at=best[1])), flush=True)

if __name__ == '__main__':
    run(int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]))
