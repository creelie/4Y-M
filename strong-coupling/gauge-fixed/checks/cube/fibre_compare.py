"""Checks of the comparison of fibres through the law of the averages
(proposition fibrecompare) and of its corollary on exteriors with a calm
completion (corollary calmexterior).

(A) checks every constant of the two proofs in exact rational arithmetic:
    n_P >= 74, rho < 1.5e-4, the bound 16 sqrt2 rho + 16 rho^2 < 22.63 rho,
    22.63 n_P rho < t/4, 14 * 90.6^3 < 1.05e7, the arithmetic of the
    corollary, and, in double precision, the sufficient condition on beta_W
    for |W| up to 10^6;
(B) on random SU(2) links checks delta <= d, delta <= 2 sqrt 2, the growth
    of delta(U_p, I)^2 when the four links of a plaquette move by at most
    rho, and the Haar mass of the ball of radius rho against rho^3/14;
Double precision in (B); the proofs in the paper are by hand.
Usage: python3 fibre_compare.py
"""
from fractions import Fraction as Fr
from math import sqrt, sin, cos, pi, log, exp
import numpy as np

checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)


# (A) constants
nP_min = Fr(49 * 6, 4)
check("n_P >= 74 (49 links, six plaquettes each, four per plaquette)",
      nP_min > 73)
rho_max = Fr(1) / (Fr(906, 10) * 74)
check("rho = t/(90.6 n_P) < 1.5e-4", rho_max < Fr(15, 100000))
# 16 sqrt2 < 22.6275, and 16 rho^2 <= 16 * 1.5e-4 rho
check("16 sqrt2 < 22.6275", Fr(226275, 10000) ** 2 > 512)
check("16 sqrt2 + 16 * 1.5e-4 < 22.63",
      Fr(226275, 10000) + 16 * Fr(15, 100000) < Fr(2263, 100))
check("22.63 / 90.6 < 1/4", Fr(2263, 100) / Fr(906, 10) < Fr(1, 4))
check("14 * 90.6^3 < 1.05e7", 14 * Fr(906, 10) ** 3 < Fr(105, 10) * 10**6)
# corollary: Q_* <= kappa^2/4, Q_W <= Q_* + kappa^2/2 = 3/4 kappa^2 < kappa^2
check("kappa^2/4 + kappa^2/2 < kappa^2", Fr(1, 4) + Fr(1, 2) < 1)
check("8.4e7 = 1.05e7 * 8", Fr(105, 10) * 10**6 * 8 == 84 * 10**6)
# lambda kappa^2 >= 196 |W| L  =>  (3/8) lam k^2 - 49|W| L >= lam k^2 / 8
check("3/8 - 49/196 = 1/8", Fr(3, 8) - Fr(49, 196) == Fr(1, 8))
K = 4.8e-18
worst = 0.0
for W in [1, 2, 3, 5, 10, 30, 100, 10**3, 10**4, 10**5, 10**6]:
    nP = 384 * W
    k2 = K**2 / (16 * nP)
    L = log(8.4e7 * nP**3 / k2**3)
    beta = 256 * 49 * nP * W * L / K**2
    worst = max(worst, beta / (W**2 * (1 + log(W))))
check(f"beta_W >= 6.4e43 |W|^2 (1+ln|W|) suffices (worst ratio {worst:.3e})",
      worst < 6.4e43)
check("K = 4.8e-18 < 4.9e-18 and kappa <= K/4 <= 1e-8",
      K < 4.9e-18 and K / 4 <= 1e-8)

# (B) SU(2)
rng = np.random.default_rng(27)


def su2(q):
    a, b, c, d = q / np.linalg.norm(q)
    return np.array([[a + 1j * b, c + 1j * d], [-c + 1j * d, a - 1j * b]])


def haar():
    return su2(rng.normal(size=4))


def expm(X):
    # X in su(2): X = i (x . sigma), exp via closed form
    w, V = np.linalg.eig(X)
    return V @ np.diag(np.exp(w)) @ np.linalg.inv(V)


def delta(A, B):
    return np.linalg.norm(A - B)


def dist(A, B):
    w = np.linalg.eigvals(A.conj().T @ B)
    th = np.abs(np.angle(w))
    return sqrt(float(np.sum(th**2)))


def small(r):
    v = rng.normal(size=3)
    v *= r * rng.random() / (np.linalg.norm(v) * sqrt(2))
    X = 1j * (v[0] * np.array([[0, 1], [1, 0]]) +
              v[1] * np.array([[0, -1j], [1j, 0]]) +
              v[2] * np.array([[1, 0], [0, -1]]))
    return expm(X)


I2 = np.eye(2)
ok1 = ok2 = True
for _ in range(20000):
    A, B = haar(), haar()
    ok1 &= delta(A, B) <= dist(A, B) + 1e-12
    ok2 &= delta(A, I2) <= 2 * sqrt(2) + 1e-12
check("delta <= d on 20000 random pairs", ok1)
check("delta(U, I) <= 2 sqrt2", ok2)
worst = 0.0
for _ in range(20000):
    rho = 1.5e-4 * rng.random()
    U = [haar() for _ in range(4)]
    V = [u @ small(rho) for u in U]
    for u, v in zip(U, V):
        assert dist(u, v) <= rho * (1 + 1e-9)
    P = U[0] @ U[1] @ U[2].conj().T @ U[3].conj().T
    Pp = V[0] @ V[1] @ V[2].conj().T @ V[3].conj().T
    worst = max(worst, (delta(Pp, I2)**2 - delta(P, I2)**2) / rho)
check(f"delta(U'_p,I)^2 - delta(U_p,I)^2 < 22.63 rho (worst {worst:.3f})",
      worst < 22.63)
okb = True
for rho in np.linspace(1e-6, 1.5e-4, 200):
    th = rho / sqrt(2)
    okb &= (2 * th - sin(2 * th)) / (2 * pi) >= rho**3 / 14
check("h(B_rho) >= rho^3/14 for rho <= 1.5e-4", okb)

print(f"{sum(checks)} of {len(checks)} checks passed")
