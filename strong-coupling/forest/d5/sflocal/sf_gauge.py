"""Gauge obstruction to e^{-lam q(log c)} times a local polymer sum.
(1) |N_1(S)| for the set S of bases of the coarse links changed by a gauge
    transformation at two adjacent blocks y, z = y + e_0 (exact count);
(2) the defect eta = log(g_y g_z^{-1}) - (log g_y - log g_z) on the link
    (y,0), in an orthonormal basis tau_a = i sigma_a / sqrt 2 of su(2) (HS
    norm), enclosed in interval arithmetic; (3) the resulting thresholds."""
import itertools
import numpy as np
from scipy.linalg import logm, expm
from mpmath import iv, mp
iv.dps = 30
# (1) exact count on a coarse torus of side 9
N = 9
def add(y, mu, k):
    y = list(y); y[mu] = (y[mu] + k) % N; return tuple(y)
def dist(a, b):
    return max(min((a[i] - b[i]) % N, (b[i] - a[i]) % N) for i in range(4))
y = (4, 4, 4, 4); z = add(y, 0, 1)
blocks = list(itertools.product(range(N), repeat=4))
changed = set()
for w in blocks:
    for mu in range(4):
        gw = w in (y, z); gwm = add(w, mu, 1) in (y, z)
        if gw or gwm:
            # link (w,mu) -> g_w c g_{w+mu}^{-1} changes iff exactly the
            # gauge values differ; with c = I it is g_w g_{w+mu}^{-1}
            changed.add((w, mu))
S = sorted(set(w for w, mu in changed))
print('changed coarse links:', len(changed), ' their bases:', len(S))
N1 = [b for b in blocks if min(dist(b, s) for s in S) <= 1]
print('blocks within l^infty distance one of a base: |N_1(S)| =', len(N1))
# (2) the defect at phi = pi/2 exactly: g_y = i s1, g_z = i s2
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.diag([1, -1]).astype(complex)
tau = [1j * s / np.sqrt(2) for s in (s1, s2, s3)]
def coords(X):
    return np.array([np.real(np.trace(t.conj().T @ X)) for t in tau])
gy, gz = 1j * s1, 1j * s2
eta = coords(logm(gy @ np.linalg.inv(gz))) - (coords(logm(gy)) - coords(logm(gz)))
print('phi = pi/2: eta in tau-coordinates', np.round(eta, 12), ' (pi/sqrt2)(-1,1,1) =', np.round(np.pi / np.sqrt(2) * np.array([-1, 1, 1]), 12))
print('|eta|^2 =', eta @ eta, ' 3 pi^2/2 =', 3 * np.pi ** 2 / 2)
# general phi: g_y = exp(sqrt2 phi tau_1), g_z = exp(sqrt2 phi tau_2)
def eta2_iv(phi):
    s, c = iv.sin(phi), iv.cos(phi)
    cpsi = c * c
    spsi = iv.sqrt(1 - cpsi * cpsi)
    psi = iv.atan2(spsi, cpsi)
    k = psi / spsi
    return 2 * (2 * (k * s * c - phi) ** 2 + (k * s * s) ** 2)
def eta2_np(phi):
    gy = expm(np.sqrt(2) * phi * tau[0]); gz = expm(np.sqrt(2) * phi * tau[1])
    e = coords(logm(gy @ np.linalg.inv(gz))) - (coords(logm(gy)) - coords(logm(gz)))
    return e @ e
bound = iv.mpf(2) * len(N1) / 10
print('total change of log of the polymer sum at most 2 |N_1(S)|/10 =', bound)
e2 = 3 * iv.pi ** 2 / 2
lam = bound / (6 * e2)
print('phi = pi/2: lambda <=', lam, '  i.e. 5.4/pi^2 =', iv.mpf('5.4') / iv.pi ** 2)
print(' phi      |eta|^2 (iv)                      check(np)      2 sin phi      lambda <= ')
for p in ['0.5', '0.2', '0.1', '0.05', '0.02']:
    ph = iv.mpf(p)
    e2 = eta2_iv(ph)
    print(p, iv.nstr(e2, 10), '%.10g' % eta2_np(float(p)), iv.nstr(2 * iv.sin(ph), 8), iv.nstr(bound / (6 * e2), 10))
# small-angle: |eta|^2 = 2 phi^4 (1 + O(phi^2))
for p in ['0.05', '0.02']:
    ph = iv.mpf(p)
    print('phi', p, ' |eta|^2/(2 phi^4) in', iv.nstr(eta2_iv(ph) / (2 * ph ** 4), 10))
