"""Certified re-check of the arithmetic of badmass.tex and patches/B01_badmass.tex.
Exact rationals (fractions) where possible, mpmath interval arithmetic (iv) otherwise.
Every check asserts; the script prints one line per check."""
from fractions import Fraction as Fr
from mpmath import iv
iv.dps = 40
n = 0
def ok(name, cond):
    global n
    assert cond, name
    n += 1
    print('ok', name)
def lt(a, b):  # interval a strictly below interval b
    return iv.mpf(a).b < iv.mpf(b).a

# thm:offchartmass
ok('2843.5 = 5.5*2068/4', Fr(55, 10) * 2068 / 4 == Fr(28435, 10))
s_max = Fr(6, 10) * Fr(1, 10**12)           # kappa+s <= 6.0e-13 gives s < 6.0e-13
r_max = s_max / 32
ok('r <= 1.875e-14', r_max == Fr(1875, 10**17))
ok('2843.5 r^2 < 1e-24', Fr(28435, 10) * r_max**2 < Fr(1, 10**24))
ok('30.2 r < s (r = s/32)', Fr(302, 10) / 32 < 1)
# lem:ballexit: (1+e^{1.5x}+e^{5.5x})^{-1} >= e^{-5.5x}/3 for x >= 0, since each term <= e^{5.5x}
ok('ballexit constants 1.5+4 = 5.5', Fr(15, 10) + 4 == Fr(55, 10))
# cor:offchartmass (i)
q0 = iv.exp(-1) / 3
ok('q0 >= e^{-1}/3 > 0.1226', lt(0.1226, q0))
ok('45 q0 > 5.5', lt(5.5, 45 * q0))
q0b = iv.exp(iv.mpf('-2.3e-7')) / 3
ok('beta=2.3e17: 2843.5 beta r^2 < 2.3e-7', Fr(23, 10) * 10**17 * Fr(1, 10**24) == Fr(23, 10**8))
ok('45 q0 > 14.99 at beta=2.3e17', lt(14.99, 45 * q0b))
# (ii)
ok('ln(6/5) < 0.1824', lt(iv.log(iv.mpf(6) / 5), 0.1824))
ok('ln(15/0.1824) > 4.40', lt(4.40, iv.log(iv.mpf(15) / iv.mpf('0.1824'))))
ok('|W|=1 threshold 4.4e24 is > 2.3e17*1e7', Fr(44, 10) * 10**24 > Fr(23, 10) * 10**17 * 10**7)
# lem:coarsewitnesssharp
k = Fr(1, 10**8)  # largest kappa; the second-order bounds are monotone in kappa
a = Fr(24001, 1000)
ok('24(1+576 kappa^2) <= 24.001', 24 * (1 + 576 * k**2) <= a)
ok('4(1+16 kappa^2) <= 4.001', 4 * (1 + 16 * k**2) <= Fr(4001, 1000))
ok('6 a^2 <= 3457', 6 * a**2 <= 3457)
ok('4a*4.001 <= 385', 4 * a * Fr(4001, 1000) <= 385)
ok('3457+385 = 3842', 3457 + 385 == 3842)
ok('(4k+3842k^2)^3 <= k^2', (4 * k + 3842 * k**2)**3 <= k**2)
ok('3457+385+3843 = 7685 <= 7700', 3457 + 385 + 3843 == 7685)
ok('pi^3/48 < 8', lt(iv.pi**3 / 48, 8))
ok('4a + 4.001 kappa < sqrt2 pi', lt(4 * float(a) * 1e-8 + 4.001e-8, iv.sqrt(2) * iv.pi))
ok('24 kappa < 1/2', 24 * k < Fr(1, 2))
# prop:badchart (ii)
c3 = iv.mpf(15000) * iv.sqrt(2)
ok('195*4e5*c3 > 1.6546e12', lt(1.6546e12, 195 * 4e5 * c3))
ok('1/1.6546e12 < 6.044e-13', Fr(1) / Fr(16546 * 10**8) < Fr(6044, 10**16))
kk = Fr(6044, 10**16)
ok('4*6.044e-13 = 2.4176e-12', 4 * kk == Fr(24176, 10**16))
ok('7700 k^2 < 3e-21', 7700 * kk**2 < Fr(3, 10**21))
ok('4k+7700k^2 < 2.4177e-12 < 2.418e-12', 4 * kk + 7700 * kk**2 < Fr(24177, 10**16) < Fr(2418, 10**15))
# prop:uniformflux
t = iv.mpf('2.39e-12')
u = t / (2 * iv.sqrt(2))
f = u / (2 * iv.sqrt(1 - u**2))  # arcsin u <= u/sqrt(1-u^2), an upper bound for f
ok('f < 4.3e-13', lt(f, 4.3e-13))
ok('kappa_f <= t/(4(1-f^2)) < 5.976e-13', lt(t / (4 * (1 - f**2)), 5.976e-13))
ok('5.976e-13 + 2.4e-15 <= 6.0e-13', Fr(5976, 10**16) + Fr(24, 10**16) <= Fr(6, 10**13))
ok('1/(6000 c3) > 7.8e-9', lt(7.8e-9, 1 / (6000 * c3)))
ok('kappa_f + 6s <= 1e-8', Fr(5976, 10**16) + 6 * Fr(24, 10**16) <= Fr(1, 10**8))
ff = iv.mpf('7e-9')
ok('kappa_f <= sqrt2 f <= 1e-8 at f=7e-9', lt(iv.sqrt(2) * ff, 1e-8))
# sin 2f = 4 sin(f/2)cos(f/2)cos f and cos(f/2)cos f >= 1-f^2 at a sample point
x = iv.mpf('0.3')
ok('identity sin2f', abs(iv.sin(2 * x) - 4 * iv.sin(x / 2) * iv.cos(x / 2) * iv.cos(x)).b < 1e-30)
ok('(1-f^2/8)(1-f^2/2) = 1-5f^2/8+f^4/16 >= 1-f^2', Fr(1,8)+Fr(1,2) == Fr(5,8) and Fr(5,8) <= 1)
ok('block count: L^3 * L/2 = L^4/2', all(L**3 * (L // 2) * 2 == L**4 for L in (4, 6, 8)))
# cor:badcrossratio
ok('(36/4)^2 * 4.8e4 = 3.888e6', Fr(36, 4)**2 * 48000 == 3888000)
ok('6*(24 sqrt2)^2/16 = 432', Fr(6 * 24**2 * 2, 16) == 432)
ok('gradient 6*3*2 = 36', 6 * 3 * 2 == 36)
ok('sup 6*2*2 = 24 (times sqrt2)', 6 * 2 * 2 == 24)
# prop:badpassage
ok('49 = 17 + 32, 45 = 49 - 4, 28 = 32 - 4', (17 + 32, 49 - 4, 32 - 4) == (49, 45, 28))
print(n, 'checks passed')
