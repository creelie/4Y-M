"""Fine radial k-scan of the transverse massive (U) symbol near k = 0 (open step 1: is inf_k attained as k -> 0?).
Hm(k) = S^{1/2} n^-4 Qhat(k) S^{1/2} on the transverse block; Qhat(k) = E_k^* (K_in+mu^2)^{-1} E_k via 4 sparse CG solves."""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from gauss_k0 import sparse_K
from gauss_massU import Khat

def run(n, mell, nt):
    K, e, live = sparse_K(n); mu = mell/n; A = (K + mu**2*sp.identity(K.shape[0], format='csr')).tocsr()
    d = live // n**4; x = np.array(np.unravel_index(live % n**4, (n,)*4)).T
    R = x + 0.5*np.eye(4)[d]
    dirs = {'axis': [1, 0, 0, 0], 'diag2': [1, 1, 0, 0], 'diag3': [1, 1, 1, 0], 'diag4': [1, 1, 1, 1], 'gen': [1, 2, 3, 5]}
    for name, v in dirs.items():
        v = np.array(v, float); v /= np.linalg.norm(v)
        for t in [0.0] + list(np.geomspace(0.02, np.pi*n, nt)):
            k = t/n*v; ph = np.exp(1j*R @ k); E = np.zeros((len(live), 4), complex); E[np.arange(len(live)), d] = ph
            X = np.column_stack([sla.cg(A, E[:, j].real, rtol=1e-10)[0] + 1j*sla.cg(A, E[:, j].imag, rtol=1e-10)[0] if np.abs(E[:, j].imag).max() > 0 else sla.cg(A, E[:, j].real, rtol=1e-10)[0] for j in range(4)])
            Q = E.conj().T @ X / n**4; Q = (Q + Q.conj().T)/2
            S = Khat(k) + mu**2*np.eye(4); w, V = np.linalg.eigh(S); Sh = V @ np.diag(np.sqrt(w)) @ V.T
            Hm = Sh @ Q @ Sh
            c = 2*np.sin(k/2)
            if c @ c > 1e-14:
                P = np.eye(4) - np.outer(c, c)/(c @ c); w, V = np.linalg.eigh(P); B = V[:, w > 0.5]; Hm = B.T @ Hm @ B
            print(json.dumps(dict(n=n, m_ell=mell, dir=name, t=t, varpi_k=float(np.linalg.eigvalsh(Hm).min()))), flush=True)

if __name__ == '__main__':
    if sys.argv[1] == 'selftest':   # must reproduce gauss_massU transverse values at small n
        run(3, 1.0, 4)
    else: run(int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]))
