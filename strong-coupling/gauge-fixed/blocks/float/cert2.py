"""Conditioning route with a shared interaction: K <= 1/4 + (1+beta) J_theta.
V(y,z) = (psi - E psi) rho23 on (x3=y, x2=z); V = Py + Pz + R; T2 = theta R + Pz, T3 = (1-theta) R + Py."""
import numpy as np
from tri import qmul, qconj, Z0, Z0p

def parts(b1, b2, b3, beta, h, g):
    x = g.x; w8 = g.w
    W = qmul(qconj(x)[:, None], qconj(x)[None, :])            # [y=x3, z=x2], w = conj(x2 x3)
    C1 = b1 + 4*beta*W; K1 = np.linalg.norm(C1, axis=-1)
    rho = np.exp(x @ b3)[:, None]*np.exp(x @ b2)[None, :]*Z0(K1); Z = w8 @ rho @ w8; rho /= Z
    psi = Z0p(K1)*(C1 @ h)/np.maximum(K1, 1e-14)/Z0(K1)
    m = w8 @ (rho*psi) @ w8
    V = rho*(psi - m)
    Py = V @ w8; Pz = w8 @ V; R = V - Py[:, None] - Pz[None, :]
    return V, Py, Pz, R, rho

def J_theta(b1, b2, b3, h, beta, g, thetas=np.linspace(0, 1, 21), cs=False, K0=None):
    V, Py, Pz, R, rho = parts(b1, b2, b3, beta, h, g)
    w8 = g.w
    cR2, cP2 = R @ g.BW.T, Pz @ g.BW.T            # rows y, harmonics in z (x2)
    cR3, cP3 = R.T @ g.BW.T, Py @ g.BW.T          # rows z, harmonics in y (x3)
    d = g.deg
    a_y = rho @ w8; a_z = w8 @ rho                # marginal densities of x3, x2
    def q(c):
        return np.maximum((c**2/d).sum(-1), 0)
    best = None
    for t in thetas:
        q2 = q(t*cR2 + cP2[None]); q3 = q((1-t)*cR3 + cP3[None])
        if cs:
            v = 0.5*np.sqrt(w8 @ (q2/a_y)) + 0.5*np.sqrt(w8 @ (q3/a_z))
        else:
            v = 0.5*(w8 @ np.sqrt(q2)) + 0.5*(w8 @ np.sqrt(q3))
        if best is None or v < best[0]: best = (v, t)
    return best

def K_cert(b1, b2, b3, h, beta, g, **kw):
    J, t = J_theta(b1, b2, b3, h, beta, g, **kw)
    return 0.25 + (1+beta)*J, J, t
