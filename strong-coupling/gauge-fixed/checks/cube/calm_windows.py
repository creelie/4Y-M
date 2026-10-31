"""Checks of the windows with a calm completion (corollary localcalm) and of
the obstruction on cubes of fixed physical size (proposition calmrare).

(A) enumerates, for every unit cube x + [0,1]^4, the blocks y whose
    neighbourhood {y}^+ contains a block holding a vertex of the cube, and
    confirms that there are at most 256; checks 16 * 256 = 4096 and
    rho = K/(3072 m) <= K/(8 n_P) for n_P <= 384 m;
(B) checks the constants of proposition calmrare: the integral of
    e^{-s^2} s^2 over [0, pi sqrt3/2] exceeds 0.44, the ratio 9.72, the
    factor (pi/(2 sqrt2))^3 < 1.371, the bound 7e-32, the logarithm 162.8
    and the threshold 10^24, the last four in exact rational arithmetic
    where they are rational;
(C) computes the probability of the cap {angle <= theta} under the density
    proportional to e^{6 beta cos t} sin^2 t by quadrature for beta from 1
    to 10^6 and rho' from 10^-4 to 1, and compares it with
    14 (beta rho'^2)^{3/2}.
Double precision in (B) and (C) where not stated; the proofs in the paper
are by hand.
Usage: python3 calm_windows.py
"""
from fractions import Fraction as Fr
from itertools import product
from math import pi, sqrt, cos, sin, exp, log, floor
import numpy as np
from scipy.integrate import quad

checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)


# (A)
worst = 0
for x in product(range(4), repeat=4):
    blocks = set()
    for v in product(*[(xi, xi + 1) for xi in x]):
        yp = tuple(vi // 2 for vi in v)
        for d in product((-1, 0, 1), repeat=4):
            blocks.add(tuple(a + b for a, b in zip(yp, d)))
    worst = max(worst, len(blocks))
check(f"at most 256 blocks y per unit cube (found {worst})", worst <= 256)
check("16 * 256 = 4096", 16 * 256 == 4096)
check("K/(3072 m) <= K/(8 * 384 m)", Fr(1, 3072) <= Fr(1, 8 * 384))

# (B)
I044 = quad(lambda s: exp(-s * s) * s * s, 0, pi * sqrt(3) / 2)[0]
check(f"int_0^(pi sqrt3/2) e^-s^2 s^2 ds = {I044:.5f} > 0.44", I044 > 0.44)
ratio = pi**2 * 3**1.5 / (3 * 4 * 0.44)
check(f"ratio pi^2 3^(3/2)/(3*4*0.44) = {ratio:.4f} <= 9.72", ratio <= 9.72)
fac = (pi / (2 * sqrt(2)))**3
check(f"(pi/(2 sqrt2))^3 = {fac:.4f} < 1.371", fac < 1.371)
check("9.72 * 1.371 < 14", Fr(972, 100) * Fr(1371, 1000) < 14)
check("n_P >= 294/4 |W| = 73.5 |W|", Fr(49 * 6, 4) == Fr(735, 10))
check("14e-24/(512 * 73.5^3) < 7e-32",
      Fr(14, 10**24) / (512 * Fr(735, 10)**3) < Fr(7, 10**32))
check("8 * 4096 * 1.05e7 < 3.45e11", 8 * 4096 * Fr(105, 10) * 10**6
      < Fr(345, 100) * 10**11)
check("3.45e11 >= 1.05e7 * 32^3", Fr(345, 100) * 10**11
      >= Fr(105, 10) * 10**6 * 32**3)
Lg = log(1.05e7 * 32**3) + 6 * log(73.5) + 48 * log(10)
check(f"ln(1.05e7 * 32^3 * 73.5^6 / 1e-48) = {Lg:.3f} > 162.8", Lg > 162.8)
check("512*49*73.5*162.8*1e16/3 > 1e24",
      512 * 49 * Fr(735, 10) * Fr(1628, 10) * 10**16 / 3 > 10**24)

# (C)
okc = True
worstr = 0.0
for beta in [1, 2, 5, 10, 100, 1e3, 1e4, 1e6]:
    def f(t):
        return exp(6 * beta * (cos(t) - 1)) * sin(t)**2
    brk = min(pi, 8 / sqrt(beta))
    den = quad(f, 0, brk, limit=200)[0] + quad(f, brk, pi, limit=200)[0]
    for rp in [1e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0]:
        th = 2 * np.arcsin(rp / (2 * sqrt(2)))
        num = quad(f, 0, min(th, pi), limit=200)[0]
        p = num / den
        b = 14 * (beta * rp * rp)**1.5
        okc &= p <= b
        if b > 0:
            worstr = max(worstr, p / b)
check(f"cap probability <= 14 (beta rho'^2)^(3/2) (largest ratio "
      f"{worstr:.3f})", okc)

print(f"{sum(checks)} of {len(checks)} checks passed")
