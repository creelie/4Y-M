"""Checks of the explicit constants c_3 and c_6 of Theorem fibrechart
(proposition explicitconstants).

(A) evaluates the chain of bounds of the proof at the chosen radii and checks
    that m! S(R)/R^m <= 1.5e4 for m = 1, 2, 3 (operator norm), and the bound
    for the derivatives of log j;
(B) adds up the constant c_6 / (dim G c_3^2) of the proof;
(C) for SU(2), on random backgrounds whose sixteen sites lie within 1e-6 of I,
    computes the derivatives of order 1 to 3 of c(X) = log(V(U)^{-1} V(U e^X))
    by finite differences along random directions Y with N(Y) = 1 and checks
    that they are far below c_3.
Double precision; the proofs in the paper are by hand.
Usage: python3 explicit_constants.py
"""
import itertools
import math
import numpy as np
from scipy.linalg import expm, logm

checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)


# (A) the chain of bounds, with a = 1e-6 and x0 = 1e-3 (rho_G <= 1)
A0, X0 = 1e-6, 1e-3
E0 = A0 + 3 * X0
G0 = X0 + 3 * X0 / (1 - E0)


def S(t):
    phi = math.expm1(t)
    psi = (1 + phi) ** 3 - 1
    if E0 + psi >= 1:
        return None
    d = psi / (1 - E0 - psi)
    eta = (1 + phi) * math.exp(d) - 1
    if G0 + eta >= 1:
        return None
    return eta / (1 - G0 - eta)


radii = {1: 1e-3, 2: 0.05, 3: 0.07}
bounds = {m: math.factorial(m) * S(R) / R ** m for m, R in radii.items()}
print("gamma_0 = %.6f" % G0, {m: round(b, 1) for m, b in bounds.items()})
check("m! S(R)/R^m <= 1.5e4 at R = 1e-3, 0.05, 0.07", max(bounds.values()) <= 1.5e4)
s3 = math.sinh(1.5) / 1.5 - 1
mg = -math.log(1 - s3)
R = 0.8
u = 2 * (X0 + R) / 3
supF = mg * u ** 2 / (1 - u ** 2)
print("sinh(1.5)/1.5 - 1 = %.5f, sup|g| <= %.5f, sup|log j| <= %.4f dim G" % (s3, mg, supF))
check("first two derivatives of log j <= dim G", s3 <= 0.4196 and mg <= 0.5441
      and supF / R <= 1 and 2 * supF / R ** 2 <= 1)

# (B) the constant c_6
kb = 34 * 0.63
a, b = 2.1, 2.1 + 2 * 744          # ||A'||_2^2 and ||A''||_1 over dim G c_3^2 T
pi1, pi2 = kb + kb ** 2, (2 * kb ** 3 + 2 * kb ** 2, kb ** 2 + kb)
T2 = 25 * b + 1275 * a
T3 = 0.5 * ((pi2[0] * a + pi2[1] * b) / 3 + (2 + 1 / 9) * 481 ** 2 * a)
total = 745 + T2 + T3
print("||Pi'|| <= %.1f ||A'||, ||Pi''|| <= %.0f ||A'||^2 + %.1f ||A''||; c_6 <= %.0f dim G c_3^2"
      % (pi1, pi2[0], pi2[1], total))
check("481 and 20576, 481 bound the derivatives of Pi", pi1 <= 481 and pi2[0] <= 20576 and pi2[1] <= 481)
check("c_6 <= 6.81e5 dim G c_3^2", total <= 6.81e5)

# SU(2) thresholds
c3s = 1.5e4 * math.sqrt(2)
c6s = 6.81e5 * 3 * c3s ** 2
kmax = 6.0e-13 / (1 + 1.24e5)
C0 = 0.5 + 2 / 3 * math.log(24 * math.pi ** 3)
bq = 1.0
while bq * kmax ** 2 <= 96 * (math.log(bq) + C0):
    bq *= 1.001
lam_calm = 45 * math.log(3e8 * 384 ** 1.5 / kmax ** 3) / (0.75 * kmax ** 2)
print("SU(2): c_3 = %.0f, c_6 = %.3e, 60 c_6 = %.3e, kappa <= %.3e, q < 1 needs beta > %.2e, "
      "calm fibres need beta > %.2e |W|" % (c3s, c6s, 60 * c6s, kmax, bq, 4 * lam_calm))
check("c_3 <= 2.13e4, c_6 <= 9.3e14, 60 c_6 < 5.6e16, 4 * 5.6e16 <= 2.3e17",
      c3s <= 2.13e4 and c6s <= 9.3e14 and 60 * c6s < 5.6e16)
check("4e5 c_3 (kappa + s) <= 1/195 when kappa + s <= 6.0e-13, and then s <= 1/(6000 c_3), kappa + 6s <= 1e-8",
      4e5 * c3s * 6.0e-13 <= 1 / 195 and 6.0e-13 <= 1 / (6000 * c3s) and 6 * 6.0e-13 <= 1e-8)
check("kappa < 4.9e-18, beta_W > 3.8e38 for q < 1, 1.5e39 |W| for calm fibres",
      kmax < 4.9e-18 and 3.8e38 <= bq and 1.5e39 <= 4 * lam_calm)

# (C) SU(2): numerical derivatives
Sg = [np.array(m, complex) for m in ([[0, 1], [-1, 0]], [[0, 1j], [1j, 0]], [[1j, 0], [0, -1j]])]
B = [s / np.sqrt(2) for s in Sg]
I2 = np.eye(2, dtype=complex)


def alg(c):
    return sum(ci * bb for ci, bb in zip(c, B))


def haar(rng):
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    return np.array([[q[0] + 1j * q[1], q[2] + 1j * q[3]], [-q[2] + 1j * q[3], q[0] - 1j * q[1]]])


def lg(g):
    Xm = logm(g)
    return (Xm - Xm.conj().T) / 2


OFFS = [o for o in itertools.product((0, 1), repeat=3) if o != (0, 0, 0)]


def V(links):
    g = links["g"]
    vs = [I2]
    for o in OFFS:
        vs += [g.conj().T @ links["u", o] @ links["l", o], g.conj().T @ links["l", o] @ links["w", o]]
    vs.append(I2)   # the straight path, through a tree link
    return g @ expm(sum(lg(v) for v in vs) / 16)


rng = np.random.default_rng(5)
worst = [0.0, 0.0, 0.0]
for trial in range(40):
    g = haar(rng)
    L = {"g": g}
    for o in OFFS:
        L["u", o] = haar(rng)
        L["l", o] = L["u", o].conj().T @ g @ expm(alg(rng.normal(size=3)) * 1e-7)
        L["w", o] = L["l", o].conj().T @ g @ expm(alg(rng.normal(size=3)) * 1e-7)
    keys = list(L)
    X = {k: alg(rng.normal(size=3)) * 1e-3 / 3 for k in keys}
    Y = {k: alg(rng.normal(size=3)) for k in keys}
    nY = math.sqrt(sum(np.real(np.trace(y.conj().T @ y)) for y in Y.values()))
    Y = {k: y / nY for k, y in Y.items()}
    V0inv = V(L).conj().T

    def c(t):
        return lg(V0inv @ V({k: L[k] @ expm(X[k] + t * Y[k]) for k in keys}))

    h = 2e-2
    f = {j: c(j * h) for j in (-2, -1, 0, 1, 2)}
    d1 = (f[1] - f[-1]) / (2 * h)
    d2 = (f[1] - 2 * f[0] + f[-1]) / h ** 2
    d3 = (f[2] - 2 * f[1] + 2 * f[-1] - f[-2]) / (2 * h ** 3)
    for i, dd in enumerate((d1, d2, d3)):
        worst[i] = max(worst[i], np.linalg.norm(dd))
print("SU(2): largest |D^m c[Y,..,Y]| for m = 1, 2, 3: %.3f, %.3f, %.3f" % tuple(worst))
c3 = 1.5e4 * math.sqrt(2)
check("the derivatives are below c_3 = 1.5e4 sqrt 2 = %.0f" % c3, max(worst) < c3)
print(f"{sum(checks)} of {len(checks)} checks passed")
