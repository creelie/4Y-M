"""Checks of the weighted operator norms (lemma weightednorms) and of the
explicit rate and amplitude of the chart decay (proposition explicitdecay).

(A) evaluates the chain of bounds of the proof at theta' = 0.197,
    theta'' = 0.6 and alpha = 2.5e-4 and checks every number stated there;
(B) on random Hermitian block matrices of finite range along a chain, checks
    |Re M^w - M| <= 2 (theta/theta')^2 |M|_(theta') with the weighted norm
    replaced by the largest value found over many Lipschitz weights, which
    is smaller than the true norm, so the test is stronger than the lemma;
(C) checks Re (s^* s)^w >= -sigma^2 sinh^2(theta) on random rows s with
    weights that move by at most one;
(D) checks, for random positive matrices of range two with Re M^w >= c,
    that the blocks of M^{-1} are at most e^{-theta d}/c.
Double precision; the proofs in the paper are by hand.
Usage: python3 explicit_decay.py
"""
from math import cosh, sinh, exp, sqrt
import numpy as np

checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)


# (A) the chain of bounds
kap = 1e-8
c3s = 1 / (4e5 * 195)
c3kp = 6 * c3s
epsA = 7.35 * c3s
RB = 2 + 24 * 260 * kap
CB = 2 * (1 / 8 + 260 * kap)
sig = (sqrt(11) / 8 + 260 * kap) / (1 / 8 - 260 * kap)
LAM60 = 60 * 6.81e5


def norms(t):
    b = exp(t) * sqrt(RB * CB)
    amb = exp(t) * sqrt(48) * epsA
    a = b + amb
    gG = 1 / (1 / 16 - 700 * kap - (cosh(2 * t) - 1) * (1 / 16 + 454 * kap))
    E = amb * b
    h = gG / (1 - gG * E)
    g = gG / (1 - gG * (2 * E + amb ** 2))
    y = b * h
    sb = sqrt(1 + 2 * sig ** 2 * exp(2 * t))
    pb = sb * (1 + y * amb)
    return b, amb, a, gG, h, g, y, sb, pb


t1, t2, alpha = 0.197, 0.6, 2.5e-4
b, amb, a, gG, h, g, y, sb, pb = norms(t1)
m = norms(t2)[6] * norms(t2)[1] * norms(t2)[7]
g2 = 1 / (1 - 2 * sig ** 2 * sinh(t1) ** 2 - 2 * m ** 2 * (t1 / t2) ** 2)
h1 = 17 / (1 - 17 * 2 * RB * epsA)
x1 = 2 * RB * h1 * exp(2 * t1) * pb ** 2
aa = 2 * exp(2 * t1) * pb ** 2
abar = 2 * (exp(2 * t1) * pb ** 2 + x1)
k0 = (pb ** 2 + x1) + (g * a * abar + (g + 2 * g ** 2 * a ** 2) * exp(2 * t1) * aa) \
    + (g2 * pb ** 2 * y * abar + (3 * g2 * pb ** 2 * y ** 2 + 2 * g2 ** 2 * pb ** 4 * y ** 2) * exp(2 * t1) * aa)
k1 = 48 * exp(t1) * pb ** 2 + 12 * c3kp * x1
K = k1 + k0 / LAM60
pa = norms(alpha)[8]
print(f"b {b:.4f}  g0 {gG:.4f}  h0 {h:.4f}  g1 {g:.4f}  y {y:.4f}  s,p {sb:.4f} {pb:.4f}  g2 {g2:.4f}")
print(f"sigma {sig:.4f}  h1 {h1:.5f}  k1 {k1:.1f}  k0 {k0:.4e}  K/lambda {K:.1f}")
print(f"alpha_max {t1 / sqrt(240 * K):.4e}  pbar(alpha) {pa:.4f}  c9 {480 * pa ** 2:.1f}")
check("c_3 s <= 1.29e-8 and c_3 kappa' <= 7.7e-8", c3s <= 1.29e-8 and c3kp <= 7.7e-8)
check("b <= e^t (1/2 + 1e-4)^(1/2) and b <= 0.862", RB * CB <= 0.5 + 1e-4 and b <= 0.862)
check("sigma <= 3.32", sig <= 3.32)
check("g0, h0, g1 <= 17.37, y <= 14.96, s, p <= 5.80, g2 <= 7.41",
      max(gG, h, g) <= 17.37 and y <= 14.96 and max(sb, pb) <= 5.80 and g2 <= 7.41)
check("h1 <= 17.001", h1 <= 17.001)
check("k1 <= 1966, k0 <= 4.16e9, K <= 2068 lambda", k1 <= 1966 and k0 <= 4.16e9 and K <= 2068)
check("2.5e-4 <= 0.197 (240 * 2068)^(-1/2)", alpha <= t1 / sqrt(240 * 2068))
check("pbar(alpha) <= 4.8 and 480 * 4.8^2 <= 1.2e4", pa <= 4.8 and 480 * 4.8 ** 2 <= 1.2e4)
check("60 c_6 = 4.086e7 dim G c_3^2", abs(LAM60 - 4.086e7) < 1e3)

rng = np.random.default_rng(5)
n, dblk = 40, 3


def lipschitz(rng, n):
    steps = rng.uniform(-1, 1, n - 1)
    return np.concatenate([[0], np.cumsum(steps)])


def conj(M, phi, t):
    w = np.repeat(np.exp(t * phi), dblk)
    return (w[:, None] * M) / w[None, :]


def banded(rng, R):
    M = np.zeros((n * dblk, n * dblk))
    for i in range(n):
        for j in range(max(0, i - R), min(n, i + R + 1)):
            M[i * dblk:(i + 1) * dblk, j * dblk:(j + 1) * dblk] = rng.normal(size=(dblk, dblk))
    return M


# (B) the Hermitian part of a conjugate
worst = 0.0
for trial in range(30):
    M = banded(rng, 2); M = (M + M.T) / 2
    tp = rng.uniform(0.1, 0.8)
    phis = [lipschitz(rng, n) for _ in range(60)] + [np.abs(np.arange(n) - k) for k in range(0, n, 5)]
    Kp = max(np.linalg.norm(conj(M, p, tp), 2) for p in phis)
    for p in phis[:20]:
        t = tp * rng.uniform(0, 1)
        Mw = conj(M, p, t)
        lhs = np.linalg.norm((Mw + Mw.T) / 2 - M, 2)
        worst = max(worst, lhs / (2 * (t / tp) ** 2 * Kp))
print(f"(B) largest ratio |Re M^w - M| / (2 (t/t')^2 |M|_(t')) {worst:.3f}")
check("Schwarz bound for the Hermitian part of a conjugate", worst <= 1)

# (C) one row of the section
worst = 0.0
for trial in range(2000):
    k = 24
    s = rng.normal(size=(3, 3 * k))
    sg = np.linalg.norm(s, 2)
    t = rng.uniform(0, 1)
    D = np.repeat(np.exp(t * rng.uniform(-1, 1, k)), 3)
    u, v = s * D, s / D
    R = (v.T @ u + u.T @ v) / 2
    worst = max(worst, -np.linalg.eigvalsh(R).min() / (sg ** 2 * sinh(t) ** 2 + 1e-300))
print(f"(C) largest -min Re(s^*s)^w / (sigma^2 sinh^2 t) {worst:.3f}")
check("Re (s^* s)^w >= -sigma^2 sinh^2 t", worst <= 1 + 1e-9)

# (D) decay of inverses
worst = 0.0
for trial in range(30):
    Bm = banded(rng, 1) * 0.15
    M = Bm @ Bm.T + 0.5 * np.eye(n * dblk)
    t = rng.uniform(0.05, 0.4)
    off = M.copy()
    for i in range(n):
        off[i * dblk:(i + 1) * dblk, i * dblk:(i + 1) * dblk] = 0
    rows = [sum(np.linalg.norm(off[i * dblk:(i + 1) * dblk, j * dblk:(j + 1) * dblk], 2) for j in range(n))
            for i in range(n)]
    c = np.linalg.eigvalsh(M).min() - (cosh(2 * t) - 1) * max(rows)
    if c <= 0:
        continue
    Mi = np.linalg.inv(M)
    for i in range(n):
        for j in range(n):
            blk = np.linalg.norm(Mi[i * dblk:(i + 1) * dblk, j * dblk:(j + 1) * dblk], 2)
            worst = max(worst, blk * exp(t * abs(i - j)) * c)
print(f"(D) largest |(M^-1)_ab| e^(t d) c {worst:.3f}")
check("blocks of the inverse at most e^(-t d)/c", worst <= 1 + 1e-9)
print(f"{sum(checks)} of {len(checks)} checks passed")
