"""Independent checks of the sharp quarter (Theorem on L(c,h) <= 1/4).

Exact parts use SymPy rationals; analytic parts use mpmath at 50 digits.
Run: python3 -I check_quarter.py
"""
import sympy as sp
import mpmath as mp
from math import comb, factorial

mp.mp.dps = 50
ok = 0


def check(cond, name):
    global ok
    if not cond:
        raise SystemExit(f"FAILED: {name}")
    ok += 1
    print("passed:", name)


# 1. coefficients of D(k) = (8+2k^2) I1^2 - k^2 I0^2 - k I1(2k)
k = sp.symbols("k")
ORDER = 30
I = lambda nu, x: sum((x / 2) ** (2 * i + nu) / (sp.factorial(i) * sp.factorial(i + nu)) for i in range(ORDER))
D = sp.expand((8 + 2 * k ** 2) * I(1, k) ** 2 - k ** 2 * I(0, k) ** 2 - k * I(1, 2 * k))
poly = sp.Poly(sp.series(D, k, 0, 2 * ORDER - 4).removeO(), k)
N = lambda m: comb(2 * m, m) * (4 * (2 * m + 1) + (m - 1) * (m + 1) * (m + 2)) - 4 ** m * (m + 1) * (m + 2)
for m in range(0, ORDER - 4):
    c = poly.coeff_monomial(k ** (2 * m + 2))
    check(c == sp.Rational(N(m), 4 ** m * factorial(m) * factorial(m + 1) * (m + 1) * (m + 2)), f"D_{m+1} formula")
check([N(m) for m in range(6)] == [0, 0, 0, 80, 1140, 10416], "N_0..N_5")
check(all(16 ** m <= 4 * m * comb(2 * m, m) ** 2 for m in range(1, 400)), "16^m <= 4m C(2m,m)^2")
check(all(N(m) > 0 for m in range(3, 400)), "N_m > 0 for 3 <= m < 400")
check(poly.coeff_monomial(k ** 8) == sp.Rational(1, 2304), "D_4 = 1/2304")

# 2. Bessel sums and the spectral parts
def parts(x):
    x = mp.mpf(x)
    I1 = mp.besseli(1, x)
    d = lambda n: (n + 1) * mp.besseli(n + 1, x) / I1
    dd = lambda n: mp.diff(lambda y: (n + 1) * mp.besseli(n + 1, y) / mp.besseli(1, y), x)
    qx = mp.nsum(lambda n: (n + 2) * d(int(n)) ** 2, [1, 60]) / (3 * x ** 2)
    qz = mp.fsum(dd(n) ** 2 / n for n in range(1, 40))
    return qx, qz


for x in (mp.mpf("0.3"), mp.mpf(1), mp.mpf("2.5")):
    S2 = mp.nsum(lambda j: j ** 2 * mp.besseli(int(j), x) ** 2, [1, mp.inf])
    S3 = mp.nsum(lambda j: j ** 3 * mp.besseli(int(j), x) ** 2, [1, mp.inf])
    S1 = mp.nsum(lambda j: mp.besseli(int(j) - 1, x) * mp.besseli(int(j), x), [1, mp.inf])
    check(abs(S2 - x * mp.besseli(1, 2 * x) / 4) < mp.mpf(10) ** -40, f"S2 at {x}")
    check(abs(S3 - x ** 2 * (mp.besseli(0, x) ** 2 + mp.besseli(1, x) ** 2) / 4) < mp.mpf(10) ** -40, f"S3 at {x}")
    check(abs(mp.besseli(1, 2 * x) - 2 * S1) < mp.mpf(10) ** -40, f"I1(2k) at {x}")
    Dv = (8 + 2 * x ** 2) * mp.besseli(1, x) ** 2 - x ** 2 * mp.besseli(0, x) ** 2 - x * mp.besseli(1, 2 * x)
    qx, qz = parts(x)
    check(abs(mp.mpf(1) / 4 - qx - Dv / (12 * x ** 2 * mp.besseli(1, x) ** 2)) < mp.mpf(10) ** -30, f"1/4 - Q_x = D/(12k^2 I1^2) at {x}")

# 3. both parts below 1/4 on (0, sqrt(88)/3], and the zonal majorant
kmax = mp.sqrt(88) / 3
grid = [kmax * j / 80 for j in range(1, 81)]
worst = 0
for x in grid:
    qx, qz = parts(x)
    A = mp.besseli(2, x) / mp.besseli(1, x)
    Ap = 1 - A ** 2 - 3 * A / x
    maj = mp.mpf(1) / 4 + 4 * Ap ** 2 + A ** 2 / 4 - A / (2 * x) - 2 * A ** 2 / x ** 2
    s = mp.sqrt(x ** 2 + mp.mpf(25) / 4)
    Ab = x / (mp.mpf(3) / 2 + s)
    assert qx <= mp.mpf(1) / 4 and qz <= maj < mp.mpf(1) / 4 and Ab <= A <= x / 4
    worst = max(worst, maj)
check(True, f"Q_x, Q_z <= majorant < 1/4 and A_flat <= A <= k/4 on 80 points (max majorant {mp.nstr(worst, 8)})")
check(4 / (3 * (1 + mp.sqrt(kmax ** 2 + 9))) <= mp.mpf(1) / 4 + mp.mpf(10) ** -40, "Cauchy-Schwarz bound reaches 1/4 at sqrt(88)/3")

# 4. the third Riccati barrier and the factorisation (exact)
kk, s = sp.symbols("kappa s", positive=True)
S = sp.sqrt(kk ** 2 + sp.Rational(25, 4))
phi = kk / (sp.Rational(3, 2) + S)
res = sp.diff(phi, kk) - (1 - phi ** 2 - 3 * phi / kk)
check(sp.simplify(res + sp.Rational(5, 2) * (S - sp.Rational(5, 2)) / (S * (sp.Rational(3, 2) + S) ** 2)) == 0, "residual of phi_flat")
a = kk / (sp.Rational(3, 2) + s)
alpha = 1 - a ** 2 - 3 * a / kk
T = a ** 2 / 4 - a / (2 * kk) - 2 * a ** 2 / kk ** 2
expr = sp.simplify((4 * alpha ** 2 + T).subs(kk, sp.sqrt(s ** 2 - sp.Rational(25, 4))))
p = 8 * s ** 3 + 28 * s ** 2 - 98 * s - 695
check(sp.simplify(expr - (2 * s - 5) * p / (4 * (2 * s + 3) ** 4)) == 0, "factorisation 4 alpha^2 + T(A_flat)")
check(sp.simplify(alpha.subs(kk, sp.sqrt(s ** 2 - sp.Rational(25, 4))) - 16 / (2 * s + 3) ** 2) == 0, "alpha = 16/(2s+3)^2")
check(p.subs(s, sp.Rational(4004, 1000)) < -120 and sp.sqrt(577) / 6 < sp.Rational(4004, 1000), "p(4.004) < -120 and sqrt(577)/6 < 4.004")

# 5. Bessel ratio monotonicity (Lemma ratio), spot values
for x in (mp.mpf("0.1"), mp.mpf(1), mp.mpf(4), mp.mpf(10)):
    r1 = mp.besseli(2, x) / mp.besseli(1, x)
    assert all(mp.besseli(n + 2, x) / mp.besseli(n + 1, x) <= r1 for n in range(1, 30))
check(True, "I_{n+2}/I_{n+1} <= I_2/I_1 at sample points")

print(f"{ok} checks passed, 0 failed")
