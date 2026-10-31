"""Checks of the lower bound for the derivative of the cut-off average along
its eight private links (lemma privatesubmersion).

For one coarse link the private links are the straight link l_0 and the seven
links l_a, a != 0, that leave B(y) in the direction mu.  With inner links u_a,
w_a the paths are P = l_0, P_x = u_a l_a and P_x' = l_a w_a.  The script
(A) checks the constants: the eigenvalue angle bound 2 arcsin(1/(2 sqrt 2)),
    Lambda = alpha/sin(alpha) < 1.0927, 1 - 7 Lambda/8 > 0.0439 and
    sigma_0 = (1 - 7 Lambda/8)/sqrt 2 > 0.031;
(B) on random fields of SU(2) and SU(3), Haar inner links and every value of
    the cut-off, applies the test variation of the proof (straight link and the
    link of the maximising path) and checks |T zeta| >= (1 - 7 Lambda/8)|zeta|,
    with derivatives by central differences;
(C) checks that the least singular value of the private derivative is at least
    sigma_0 on the same fields.
Double precision; the proof in the paper is by hand.
Usage: python3 private_submersion.py
"""
import itertools
import numpy as np
from scipy.linalg import expm, logm

checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)


def basis(N):
    """Orthonormal basis of su(N) for <X,Y> = Re tr(X^* Y)."""
    B = []
    for j in range(N):
        for k in range(j + 1, N):
            E = np.zeros((N, N), complex); E[j, k] = 1; E[k, j] = -1
            B.append(E / np.sqrt(2))
            E = np.zeros((N, N), complex); E[j, k] = 1j; E[k, j] = 1j
            B.append(E / np.sqrt(2))
    for d in range(1, N):
        D = np.zeros((N, N), complex)
        D[:d, :d] = np.eye(d) * 1j; D[d, d] = -1j * d
        B.append(D / np.linalg.norm(D))
    return B


def coords(B, X):
    return np.array([np.real(np.trace(b.conj().T @ X)) for b in B])


def alg(B, c):
    return sum(ci * b for ci, b in zip(c, B))


def haar(N, rng):
    Z = (rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    Q = Q @ np.diag(np.diag(R) / np.abs(np.diag(R)))
    return Q / np.linalg.det(Q) ** (1 / N)


def lg(g):
    X = logm(g)
    return (X - X.conj().T) / 2


def phi(t):
    return min(1.0, max(0.0, 2 - 2 * t))


def cutoff(l, u, w, rho):
    P = l[0]
    Ps = []
    for a in range(1, 8):
        Ps += [u[a] @ l[a], l[a] @ w[a]]
    d = [np.linalg.norm(Px - P) for Px in Ps]
    m = max(d)
    chi = phi(m / rho)
    if chi == 0:
        return P, chi, m, None
    M = sum(lg(P.conj().T @ Px) for Px in Ps) / 16
    return P @ expm(chi * M), chi, m, int(np.argmax(d))


def right_deriv(B, f, h=1e-6):
    """Right-trivialised derivative of t -> f(t) at 0 in the basis B."""
    Vp, Vm, V0 = f(h), f(-h), f(0.0)
    return coords(B, lg(Vp @ V0.conj().T) - lg(Vm @ V0.conj().T)) / (2 * h)


ALPHA = 2 * np.arcsin(1 / (2 * np.sqrt(2)))
LAM = ALPHA / np.sin(ALPHA)
MARGIN = 1 - 7 * LAM / 8
SIGMA0 = MARGIN / np.sqrt(2)
check(f"alpha = 2 arcsin(1/(2 sqrt 2)) = {ALPHA:.6f} < 1.4455/2", ALPHA < 1.4455 / 2)
check(f"Lambda = {LAM:.6f} < 1.0927", LAM < 1.0927)
check(f"1 - 7 Lambda/8 = {MARGIN:.6f} > 0.0439", MARGIN > 0.0439)
check(f"sigma_0 = {SIGMA0:.6f} > 0.031", SIGMA0 > 0.031)
# two eigenvalue angles a, b with 4 sin^2(a/2) + 4 sin^2(b/2) <= 1 have a + b <= 2 alpha
ts = np.linspace(0, np.pi, 200001)
s = 4 * np.sin(ts / 2) ** 2
bmax = np.where(s <= 1, 2 * np.arcsin(np.sqrt(np.clip(1 - s, 0, 1)) / 2), -1)
check("max a + b under 4sin^2(a/2)+4sin^2(b/2) <= 1 equals 2 alpha",
      abs(np.max(np.where(bmax >= 0, ts + bmax, 0)) - 2 * ALPHA) < 1e-4)


def sample(N, rng, rho):
    B = basis(N)
    while True:
        l0 = haar(N, rng)
        sc = rng.uniform(0.05, 0.75)
        l, u, w = [l0], [None], [None]
        for a in range(1, 8):
            ua = haar(N, rng)
            v1 = expm(alg(B, rng.normal(size=len(B)) * sc / np.sqrt(len(B))))
            v2 = expm(alg(B, rng.normal(size=len(B)) * sc / np.sqrt(len(B))))
            la = ua.conj().T @ l0 @ v1
            l.append(la); u.append(ua); w.append(la.conj().T @ l0 @ v2)
        V, chi, m, xs = cutoff(l, u, w, rho)
        if chi > 0:
            return B, l, u, w, chi, m, xs


worst_T, worst_s, counts = 9.0, 9.0, {"switch": 0, "average": 0}
for N, rho, trials in ((2, 1.0, 150), (3, 1.0, 60)):
    rng = np.random.default_rng(N)
    gT, gs = 9.0, 9.0
    for t in range(trials):
        B, l, u, w, chi, m, xs = sample(N, rng, rho)
        counts["switch" if chi < 1 else "average"] += 1
        a_star = xs // 2 + 1
        first = xs % 2 == 0          # maximiser is u_a l_a
        P = l[0]
        zeta = rng.normal(size=len(B)); zeta /= np.linalg.norm(zeta)
        Z = alg(B, zeta)
        eta = P @ Z @ P.conj().T     # P-frame vector zeta, moved to the identity

        def f(tt):
            ll = list(l)
            ll[0] = expm(tt * eta) @ l[0]
            if chi < 1:
                if first:   # eta'_{x} = zeta at x = (a*, 1)
                    xi = u[a_star].conj().T @ eta @ u[a_star]
                else:       # eta'_{x'} = zeta
                    xi = eta
                ll[a_star] = expm(tt * xi) @ l[a_star]
            return cutoff(ll, u, w, rho)[0]

        out = right_deriv(B, f)
        gT = min(gT, np.linalg.norm(out))
        cols = []
        for a in range(8):
            for k in range(len(B)):
                def g(tt, a=a, k=k):
                    ll = list(l); ll[a] = expm(tt * B[k]) @ l[a]
                    return cutoff(ll, u, w, rho)[0]
                cols.append(right_deriv(B, g))
        sv = np.linalg.svd(np.array(cols).T, compute_uv=False)[-1]
        gs = min(gs, sv)
    worst_T, worst_s = min(worst_T, gT), min(worst_s, gs)
    print(f"SU({N}): fields {trials}, least |T zeta| {gT:.4f}, least sigma {gs:.4f}")
check(f"test variation |T zeta| >= 1 - 7 Lambda/8 on every sampled field "
      f"({counts['switch']} in the switch region, {counts['average']} averaging)",
      worst_T >= MARGIN)
check(f"least singular value of the private derivative >= sigma_0 (least {worst_s:.4f})",
      worst_s >= SIGMA0)
print(f"{sum(checks)} of {len(checks)} checks passed")
