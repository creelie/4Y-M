"""Nystrom check of Theorem A2 (finite-volume gap and Doeblin bound) for the
smallest nontrivial case: 3D SU(2) Wilson theory, spatial 2-torus of side L=1
(one site, two links, one plaquette U1 U2 U1^-1 U2^-1).

Transfer kernel on L^2(SU(2)^2):
  k(U,V) = exp(beta/2 (Re tr P(U) + Re tr P(V)))
           * int dh exp(beta Re tr(U1 h V1^-1 h^-1) + beta Re tr(U2 h V2^-1 h^-1))
Theorem A2 predicts lambda1/lambda0 <= 1 - exp(-2 beta N (|P|+|E|)), N=2, |P|=1, |E|=2.
Run: python3 -I torus_gap_check.py
"""
import numpy as np

rng = np.random.default_rng(1)

def haar(n):
    q = rng.normal(size=(n, 4))
    return q / np.linalg.norm(q, axis=1, keepdims=True)

def mul(a, b):
    a0, av = a[..., 0], a[..., 1:]
    b0, bv = b[..., 0], b[..., 1:]
    s = a0 * b0 - np.sum(av * bv, axis=-1)
    v = a0[..., None] * bv + b0[..., None] * av + np.cross(av, bv)
    return np.concatenate([s[..., None], v], axis=-1)

def conj(a):
    return np.concatenate([a[..., :1], -a[..., 1:]], axis=-1)

def run(beta, n=700, m=60):
    U1, U2 = haar(n), haar(n)
    h = haar(m)
    h = np.concatenate([h, conj(h)])  # symmetric sample in h <-> h^-1
    plaq = 2 * mul(mul(U1, U2), mul(conj(U1), conj(U2)))[:, 0]  # Re tr
    def W(V):  # h V^-1 h^-1, shape (n, 2m, 4)
        return mul(mul(h[None, :, :], conj(V)[:, None, :]), conj(h)[None, :, :])
    S = np.zeros((n, n, h.shape[0]))
    for U, V in ((U1, U1), (U2, U2)):
        wb = conj(W(V))  # Re tr(U W) = 2 <U, conj W>
        S += 2 * np.einsum('ia,kja->ikj', U, wb)
    # factor out the maximum to avoid overflow; kernel ratios are unchanged
    kin = np.exp(beta * (S - 4)).mean(axis=2)
    K = np.exp(beta / 2 * (plaq[:, None] + plaq[None, :] - 4)) * kin
    K = (K + K.T) / 2 / n
    ev = np.sort(np.linalg.eigvalsh(K))[::-1]
    ratio = ev[1] / ev[0]
    bound = 1 - np.exp(-2 * beta * 2 * (1 + 2))
    return ratio, bound

if __name__ == '__main__':
    for beta in (0.25, 0.5, 1.0, 2.0, 4.0):
        r, b = run(beta)
        print(f"beta={beta:5.2f}  lambda1/lambda0={r:.4f}  Doeblin bound={b:.10f}  ok={r <= b}"
              f"  gap={-np.log(r):.4f}  bound_gap={-np.log(b):.3e}")
