"""Slice ratio for composite, gauge-invariant slice fields of the free lattice gauge field.

Remark freeslice of the forest paper treats the linear field. The slice field of Sec.
"Renormalised flowed observables" is a gauge-invariant function of the spatial links of one
slice, composed with the spatial flow. Its free (abelian, one colour; N^2-1 colours multiply D by
N^2-1) analogue is the magnetic energy  Phi_B(n) = sum_x sum_{i<j spatial} F_ij(n,x)^2.  For
comparison we also treat the full energy Phi_E+B, which contains F_0i and is not a function of
the spatial links of one slice.

F = forward-difference curl of A on Z^4, Feynman gauge, and the lattice heat flow of the spatial
links for a physical time t_f (factor exp(-t_f phat^2) per field, phat the spatial momentum).
For a Gaussian field Cov(F^2, F^2) = 2 G^2, so per unit spatial volume (zero total momentum)
    D(n) = (2/V) sum_{p in BZ} sum_{(mu nu), (rho sigma)} |Gt_{mu nu, rho sigma}(n, p)|^2 ,
Gt(n,p) the mixed (time, spatial momentum) propagator of F. With d_mu = e^{i p_mu} - 1,
    <F_{mu nu} F_{rho sigma}~>(p) = (d_mu db_rho dl_{nu sigma} - d_mu db_sigma dl_{nu rho}
                                    - d_nu db_rho dl_{mu sigma} + d_nu db_sigma dl_{mu rho}) / |d|^2;
the numerator is c_{-1} e^{-i th} + c_0 + c_1 e^{i th} in the time momentum th, and the inverse
Fourier transform of 1/(2 - 2cos th + m^2) is K(n) = e^{-w|n|}/(2 sinh w), cosh w = 1 + m^2/2,
so Gt(n) = sum_j c_j K(n+j) exactly. For the magnetic field the numerator does not depend on th.

The continuum reference is the same sum over the momenta 2 pi n/ell of the torus of side ell,
with Gt(r,p) = N(p0 = i|p|) e^{-|p| r}/(2|p|) for r > 0; the infinite-volume continuum value is
also printed. Without the flow the continuum ratio is (r/r')^5 for both fields.

    python3 composite_free.py
"""
import itertools

import numpy as np

SPACE = [(1, 2), (1, 3), (2, 3)]
ALL = [(m, n) for m in range(4) for n in range(m + 1, 4)]


def numer(d, db, mu, nu, rho, sg):
    dl = lambda a, b: 1.0 if a == b else 0.0
    return (d[:, mu] * db[:, rho] * dl(nu, sg) - d[:, mu] * db[:, sg] * dl(nu, rho)
            - d[:, nu] * db[:, rho] * dl(mu, sg) + d[:, nu] * db[:, sg] * dl(mu, rho))


def D_lattice(L, ns, tf_lat, pairs):
    """D(n) for n in ns on an L^3 spatial torus (lattice units), spatial flow time tf_lat."""
    k = 2 * np.pi * np.fft.fftfreq(L)
    p = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1).reshape(-1, 3)
    p = p[np.any(p != 0, axis=1)]
    ds = np.exp(1j * p) - 1.0
    m2 = np.sum(np.abs(ds) ** 2, axis=1)
    w = np.arccosh(1 + m2 / 2)
    flow = np.exp(-2 * tf_lat * m2)
    K = lambda n: np.exp(-w * np.abs(n)) / (2 * np.sinh(w))
    ths = 2 * np.pi * np.arange(3) / 3
    ds_all = [np.concatenate([np.full((len(p), 1), np.exp(1j * th) - 1), ds], axis=1) for th in ths]
    out = np.zeros(len(ns))
    for (mu, nu), (rho, sg) in itertools.product(pairs, pairs):
        vals = np.array([numer(d, np.conj(d), mu, nu, rho, sg) for d in ds_all])
        c = {j: np.mean(vals * np.exp(-1j * j * ths)[:, None], axis=0) for j in (-1, 0, 1)}
        for i, n in enumerate(ns):
            out[i] += np.sum(np.abs(sum(c[j] * K(n + j) for j in (-1, 0, 1)) * flow) ** 2)
    return 2 * out / L ** 3


def D_continuum(ell, rs, tf, pairs, nmax=None):
    """D(r) for r > 0 in the continuum; ell = None is infinite volume (angular average, then a
    radial quadrature), otherwise the torus of side ell."""
    if ell is None:
        rng = np.random.default_rng(0)
        u = rng.normal(size=(20000, 3))
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        d = np.concatenate([-np.ones((len(u), 1)), 1j * u], axis=1)
        db = np.concatenate([np.ones((len(u), 1)), -1j * u], axis=1)
        S = sum(np.abs(numer(d, db, *a, *b)) ** 2 for a, b in itertools.product(pairs, pairs)).mean()
        q = np.linspace(1e-4, 80, 40001)
        return np.array([np.trapezoid(S * q ** 4 / 4 * np.exp(-2 * q * r - 4 * tf * q * q), q)
                         for r in rs]) * 2 / (2 * np.pi ** 2)
    nmax = nmax or int(ell * 10)
    k = 2 * np.pi * np.arange(-nmax, nmax + 1) / ell
    p = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1).reshape(-1, 3)
    p = p[np.any(p != 0, axis=1)]
    w = np.linalg.norm(p, axis=1)
    d = np.concatenate([(-w)[:, None], 1j * p], axis=1)     # d_0 = i p0 at p0 = i|p|
    db = np.concatenate([w[:, None], -1j * p], axis=1)      # conjugate taken before continuation
    S = sum(np.abs(numer(d, db, *a, *b)) ** 2 for a, b in itertools.product(pairs, pairs))
    return np.array([np.sum(S * np.exp(-2 * w * r - 4 * tf * w * w) / (4 * w * w)) for r in rs]) * 2 / ell ** 3


def main():
    rs = np.array([0.25, 0.5, 1.0])
    ell = 6.0
    for name, pairs, spacings in (("magnetic (the paper's slice field)", SPACE, (8, 12, 16, 24)),
                                  ("electric + magnetic", ALL, (8, 12, 16))):
        print(f"== {name}")
        for tf in (0.05, 0.0):
            Di = D_continuum(None, rs, tf, pairs)
            Df = D_continuum(ell, rs, tf, pairs)
            print(f"t_f = {tf}: continuum D(.5)/D(.25), D(1)/D(.5) = {Di[1]/Di[0]:.5f}, {Di[2]/Di[1]:.5f}"
                  f" (infinite volume); {Df[1]/Df[0]:.5f}, {Df[2]/Df[1]:.5f} (side {ell})", flush=True)
            for inv_a in spacings:
                L = int(ell * inv_a)
                ns = np.array([0, 1] + [int(round(r * inv_a)) for r in rs])
                D = D_lattice(L, ns, tf * inv_a ** 2, pairs)
                print(f"   1/a={inv_a:3d} L={L:3d}: D(.5)/D(.25)={D[3]/D[2]:.5f} D(1)/D(.5)={D[4]/D[3]:.5f}"
                      f"  D(a)/D(0)={D[1]/D[0]:.4f}  D(.25)/D(0)={D[2]/D[0]:.3e}", flush=True)


if __name__ == "__main__":
    main()
