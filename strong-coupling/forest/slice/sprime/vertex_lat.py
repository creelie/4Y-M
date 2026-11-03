"""Quadratic vertex of the SU(2) Wilson-flow drift on a 3D lattice slice (check of Lemma 3 in LATTICE_SPRIME.md).

Drift Z_mu(x) = vec(U_mu(x) * sum of staples) (the su(2) part of the plaquette force, up to sign/normalisation).
For plane waves A1 = cos(k.x) e1 t^a, A2 = cos(l.x) e2 t^b on the links, the bilinear part
Q = d^2 Z / d alpha d beta at 0 is computed by a symmetric finite difference (exact up to O(h^2)), Fourier
transformed, and its amplitude at k + l compared with |khat| + |lhat| + |Phat|.
    python3 vertex_lat.py [L NPAIRS NMAX]   (NMAX: momentum components n_i in [-NMAX, NMAX], default all)
"""
import sys
import numpy as np

L = int(sys.argv[1]) if len(sys.argv) > 1 else 8
NP = int(sys.argv[2]) if len(sys.argv) > 2 else 400
NMAX = int(sys.argv[3]) if len(sys.argv) > 3 else L // 2
rng = np.random.default_rng(1)
x = np.stack(np.meshgrid(*([np.arange(L)] * 3), indexing="ij"))


def qmul(a, b):
    a0, av = a[..., 0], a[..., 1:]
    b0, bv = b[..., 0], b[..., 1:]
    return np.concatenate([(a0 * b0 - np.sum(av * bv, -1))[..., None],
                           a0[..., None] * bv + b0[..., None] * av + np.cross(av, bv)], -1)


def qconj(a):
    return a * np.array([1, -1, -1, -1])


def qexp(v):                                   # v: (..., 3) -> unit quaternion exp(v . i sigma / 2)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    s = np.where(n > 0, np.sin(n / 2) / np.where(n > 0, n, 1), 0.5)
    return np.concatenate([np.cos(n / 2), s * v], -1)


def sh(f, mu, d):                              # f(x + d e_mu)
    return np.roll(f, -d, axis=mu)


def drift(A):                                  # A: (3 dirs, L, L, L, 3 colours)
    U = qexp(A)
    Z = np.zeros_like(A)
    for mu in range(3):
        S = np.zeros(U.shape[1:])
        for nu in range(3):
            if nu == mu:
                continue
            Unu_xmu = sh(U[nu], mu, 1)
            S += qmul(qmul(Unu_xmu, qconj(sh(U[mu], nu, 1))), qconj(U[nu]))
            Unu_xmu_mnu = sh(sh(U[nu], mu, 1), nu, -1)
            S += qmul(qmul(qconj(Unu_xmu_mnu), qconj(sh(U[mu], nu, -1))), sh(U[nu], nu, -1))
        Z[mu] = qmul(U[mu], S)[..., 1:]
    return Z


def wave(k, d, c):
    f = np.zeros((3, L, L, L, 3))
    f[d, ..., c] = np.cos(np.tensordot(k, x, axes=1))
    return f


h = 1e-3
khat = lambda k: np.sqrt(np.sum(4 * np.sin(k / 2) ** 2))
rows = []
for _ in range(NP):
    n1, n2 = rng.integers(-NMAX, NMAX + 1, 3) % L, rng.integers(-NMAX, NMAX + 1, 3) % L
    k, l = 2 * np.pi * n1 / L, 2 * np.pi * n2 / L
    d1, d2 = rng.integers(0, 3, 2)
    c1, c2 = rng.choice(3, 2, replace=False)  # different colours: the abelian part drops out of the bilinear term
    a, b = h * wave(k, d1, c1), h * wave(l, d2, c2)
    Q = (drift(a + b) - drift(a - b) - drift(-a + b) + drift(-a - b)) / (4 * h * h)
    Qf = np.fft.fftn(Q, axes=(1, 2, 3)) / L ** 3
    m = tuple(((n1 + n2) % L).tolist())
    amp = np.linalg.norm(Qf[(slice(None),) + m])
    P = k + l
    rows.append((khat(k), khat(l), khat(P), amp))
r = np.array(rows)
den = r[:, 0] + r[:, 1] + r[:, 2]
ok = den > 0
print(f"L = {L}, NMAX = {NMAX}, {NP} random pairs (k, l), plane waves of unit amplitude in different colours")
print(f"max |Q(k+l)| over pairs with k = l = 0: {r[~ok, 3].max() if (~ok).any() else 'none sampled'}")
ratio = r[ok, 3] / den[ok]
print(f"max |Q(k+l)| / (|khat|+|lhat|+|Phat|) = {ratio.max():.4f}")
for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 11)):
    s = (den[ok] >= lo) & (den[ok] < hi)
    if s.any():
        print(f"  |khat|+|lhat|+|Phat| in [{lo},{hi}): {s.sum():4d} pairs, max ratio {ratio[s].max():.4f}")
# constant fields: the bilinear term must vanish
Q0 = []
for d1 in range(3):
    for d2 in range(3):
        a, b = h * wave(np.zeros(3), d1, 0), h * wave(np.zeros(3), d2, 1)
        Q0.append(np.abs((drift(a + b) - drift(a - b) - drift(-a + b) + drift(-a - b)) / (4 * h * h)).max())
print(f"constant fields, all direction pairs: max |Q| = {max(Q0):.2e}")
