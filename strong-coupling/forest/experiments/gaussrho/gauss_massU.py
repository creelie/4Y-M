"""Does a gap imply (U)? (A_GAUSSIAN_FIBRE_4D §8.4) Free massive abelian field (Proca, lattice): precision K + mu^2,
K = d1^T d1. First chaos: f = u.A, Var(f) = u^*(K^+mu^2)^{-1}u. Infinite volume, Fourier symbols (4x4):
  E_n(f) symbol  = n^-4 Qhat(k),  Qhat(k)_{mu nu} = sum_{i,j in cube} [(K_in+mu^2)^{-1}]_{ij} e^{ik(r_i-r_j)} (dir mu_i=mu, mu_j=nu)
  varpi_n = inf_k lambda_min( S^{1/2} n^-4 Qhat S^{1/2} ),  S = Khat + mu^2.
Lattice spacing a = ell/n, mu = m a = (m ell)/n. (U) holds iff varpi_n stays bounded below as n grows at fixed m*ell."""
import sys, json, itertools, numpy as np
from gauss_rho import cube_K

def Khat(k):
    c = 2*np.sin(k/2); return np.sum(c**2)*np.eye(4) - np.outer(c, c)  # midpoint convention for link positions

def run(n, mell, npts, transverse=False):
    mu = mell/n; inl, Kin = cube_K(1, n+4, n, (0,)*4)
    Ci = np.linalg.inv(Kin + mu**2*np.eye(len(inl)))
    R = np.array([np.array(x)+0.5*np.eye(4)[d] for (x, d) in inl]); D = np.array([d for (x, d) in inl])
    E = np.zeros((4, len(inl))); E[D, np.arange(len(inl))] = 1
    ks = [0.0] + list(np.geomspace(1e-3, np.pi, npts-1))
    best = (np.inf, None)
    for k in itertools.product(ks, repeat=4):
        if list(k) != sorted(k): continue  # axis-permutation symmetry
        k = np.array(k); ph = np.exp(1j*R @ k)
        Q = (E*ph) @ Ci @ (E*ph).conj().T / n**4
        S = Khat(k) + mu**2*np.eye(4); w, V = np.linalg.eigh(S); Sh = V @ np.diag(np.sqrt(w)) @ V.conj().T
        Hm = Sh @ Q @ Sh
        if transverse:  # gauge-invariant (divergence-free) u only: u orthogonal to c(k); undefined at k = 0
            c = 2*np.sin(k/2)
            if c @ c < 1e-14: continue
            P = np.eye(4) - np.outer(c, c)/(c @ c); w, V = np.linalg.eigh(P); B = V[:, w > 0.5]; Hm = B.T @ Hm @ B
        lm = np.linalg.eigvalsh(Hm).min()
        if lm < best[0]: best = (float(lm), k.tolist())
    vp = best[0]
    print(json.dumps(dict(n=n, m_ell=mell, transverse=transverse, mu=mu, varpi=vp, kmin=best[1], gapbound_over_m=16*vp/(625*np.e**2*mell))), flush=True)

def selftest():  # Khat check: Khat(k) annihilates pure-gauge s-direction and matches finite-difference K on a small torus
    k = np.array([0.3, 0.7, 1.1, 0.2]); g = 2*np.sin(k/2)
    print('gauge null', np.abs(Khat(k) @ g).max()); assert np.abs(Khat(k) @ g).max() < 1e-12
    from gauss_angle import full_K
    N = 4; K = full_K(N).toarray(); kk = 2*np.pi*np.array([1, 0, 2, 3])/N
    from gauss_rho import sites, sidx
    X = sites(N); v = np.zeros((4*N**4, 4), complex)
    for d in range(4): v[sidx(X, N)*4+d, d] = np.exp(1j*(X + 0.5*np.eye(4)[d]) @ kk)
    M = v.conj().T @ K @ v / N**4; err = np.abs(M - Khat(kk)).max(); print('Khat vs K', err); assert err < 1e-10
    print('SELFTEST PASS')

if __name__ == '__main__':
    if sys.argv[1] == 'selftest': selftest()
    else: run(int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), len(sys.argv) > 4 and sys.argv[4] == 't')
