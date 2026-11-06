"""Numerical sanity check for SHELL_TRANSPORT (SU(2), rho_G = 1).

Paths: P_2y = l0, P_{x_b} = u_b l_b, P_{x'_b} = l_b w_b (b = 0..7, u_0 = w_0 = I).
v_x = P_2y^{-1} P_x, m = max_x |v_x - I|_HS, chi = phi(m), phi(t) = min(1, max(0, 2 - 2t)),
V# = l0 exp(chi/16 sum_x log v_x).
J  = right-trivialised derivative of V# in the eight private links l_0..l_7 (3 x 24).
For an inner link j in {u_b, w_b}, dc = derivative of V# along j, Y = -J^T (J J^T)^{-1} dc.
Checks: sigma_min(J) >= sigma_0 = (1 - 7 Lambda/8)/sqrt2 and |Y|_2 <= (Lambda/16 + 2 pi/3)/sigma_0.
Floating point, finite differences: heuristic evidence only.
"""
import json, sys, time
import numpy as np

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
N_SAMPLES = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
I2 = np.eye(2, dtype=complex)
TAU = [np.array(a, dtype=complex) for a in (
    [[0, 1j], [1j, 0]], [[0, 1], [-1, 0]], [[1j, 0], [0, -1j]])]
TAU = [t / np.sqrt(2) for t in TAU]   # HS-orthonormal basis of su(2): <A,B> = Re tr(A B^*)

alpha = 2 * np.arcsin(1 / (2 * np.sqrt(2)))
LAM = alpha / np.sin(alpha)
SIG0 = (1 - 7 * LAM / 8) / np.sqrt(2)
YBOUND = (LAM / 16 + 2 * np.pi / 3) / SIG0

def alg(x): return sum(c * t for c, t in zip(x, TAU))
def coords(A): return np.array([np.real(np.trace(A @ t.conj().T)) for t in TAU])
def expm(A):
    x = coords(A); th = np.linalg.norm(x) / np.sqrt(2)   # A^2 = -th^2 I
    return I2 * np.cos(th) + (A * np.sin(th) / th if th > 1e-15 else A)
def logm(U):
    c = np.clip(np.real(np.trace(U)) / 2, -1, 1); th = np.arccos(c)
    B = (U - U.conj().T) / 2
    return B * (th / np.sin(th) if th > 1e-12 else 1.0)
def rand_su2(scale):
    return expm(alg(rng.normal(size=3) * scale))
def haar():
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    return np.array([[q[0] + 1j * q[1], q[2] + 1j * q[3]], [-q[2] + 1j * q[3], q[0] - 1j * q[1]]])

def vsharp(l, u, w):
    l0 = l[0]; inv = l0.conj().T
    vs = []
    for b in range(8):
        vs.append(inv @ u[b] @ l[b]); vs.append(inv @ l[b] @ w[b])
    m = max(np.linalg.norm(v - I2) for v in vs)
    chi = min(1.0, max(0.0, 2 - 2 * m))
    if chi == 0: return l0, m
    N = sum(logm(v) for v in vs) / 16
    return l0 @ expm(chi * N), m

def deriv(f, base, h=1e-6):
    """right-trivialised derivative: X -> d/dt f(e^{tX} g) f^{-1}, central differences."""
    V0 = f(base)[0]; cols = []
    for k in range(3):
        Vp = f(expm(h * TAU[k]) @ base)[0]; Vm = f(expm(-h * TAU[k]) @ base)[0]
        cols.append(coords((Vp - Vm) @ V0.conj().T / (2 * h)))
    return np.array(cols).T

def sample(mode):
    l0 = haar()
    s = {'calm': 0.05, 'shell': 0.25, 'wide': 0.6}[mode]
    l = [l0] + [rand_su2(s) @ l0 for _ in range(7)]
    u = [I2] + [rand_su2(s) for _ in range(7)]
    w = [I2] + [rand_su2(s) for _ in range(7)]
    return l, u, w

stats = {}
t0 = time.time()
for mode in ('calm', 'shell', 'wide'):
    worst_sig, worst_Y, counts = np.inf, 0.0, {'m<1/2': 0, 'shell': 0, 'm>1': 0}
    for _ in range(N_SAMPLES):
        l, u, w = sample(mode)
        V, m = vsharp(l, u, w)
        if abs(m - 0.5) < 1e-4 or abs(m - 1) < 1e-4: continue
        key = 'm<1/2' if m < 0.5 else ('shell' if m < 1 else 'm>1')
        counts[key] += 1
        J = np.hstack([deriv(lambda g, k=k: vsharp(l[:k] + [g] + l[k + 1:], u, w), l[k]) for k in range(8)])
        sig = np.linalg.svd(J, compute_uv=False).min()
        worst_sig = min(worst_sig, sig)
        G = np.linalg.inv(J @ J.T)
        for b in range(1, 8):
            for which in ('u', 'w'):
                if which == 'u':
                    dc = deriv(lambda g: vsharp(l, u[:b] + [g] + u[b + 1:], w), u[b])
                else:
                    dc = deriv(lambda g: vsharp(l, u, w[:b] + [g] + w[b + 1:]), w[b])
                Y = -J.T @ G @ dc            # 24 x 3 : transport for each unit direction
                worst_Y = max(worst_Y, np.linalg.norm(Y, 2))
    stats[mode] = {'counts': counts, 'min_sigma': float(worst_sig), 'max_Y_opnorm': float(worst_Y)}

out = {'sigma0': SIG0, 'Lambda': LAM, 'Y_bound': YBOUND, 'samples_per_mode': N_SAMPLES,
       'stats': stats, 'seconds': time.time() - t0,
       'pass': all(s['min_sigma'] >= SIG0 - 1e-6 and s['max_Y_opnorm'] <= YBOUND + 1e-6
                   for s in stats.values())}
print(json.dumps(out, indent=1))
