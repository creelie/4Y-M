"""Certified constants for the ratio tuning (nulower.tex).

Exact rationals: the two-loop exponent b1/(2 b0^2) for SU(N), with
b0 = 11N/(48 pi^2) and b1 = 34N^2/(3 (16 pi^2)^2); the powers of pi cancel,
so b1/(2 b0^2) is computed from the rational parts.

Interval arithmetic (mpmath.iv, 60 digits): for beta_W = 2N/g^2,
  b     = 4 N b0 = 11 N^2/(12 pi^2)      (one-loop slope of beta_W in log(1/a))
  h     = b log 2                        (one-loop step per halving of a)
  gamma = (51/121) b                     (coefficient of log beta_W at two loops)
  c     = gamma h                        (constant of the step-scaling inequality)
for N = 2 and N = 3.
"""
from fractions import Fraction as Fr
from mpmath import iv

iv.dps = 60

print("Exact two-loop exponent b1/(2 b0^2) for SU(N):")
for N in range(2, 11):
    b0 = Fr(11 * N, 48)            # times 1/pi^2
    b1 = Fr(34 * N * N, 3 * 256)   # times 1/pi^4
    e = b1 / (2 * b0 * b0)
    print(f"  N={N}: {e}")
    assert e == Fr(51, 121)
# b = 4 N b0 should equal 11 N^2/12 (times 1/pi^2)
for N in range(2, 11):
    assert 4 * N * Fr(11 * N, 48) == Fr(11 * N * N, 12)
print("  all equal 51/121; b = 4 N b0 = 11 N^2/(12 pi^2) checked for N=2..10")

pi = iv.pi
log2 = iv.log(2)
for N in (2, 3):
    b = iv.mpf(11 * N * N) / (12 * pi**2)
    h = b * log2
    gam = iv.mpf(51) / 121 * b
    c = gam * h
    print(f"SU({N}):")
    for name, x in (("b", b), ("h", h), ("gamma", gam), ("c=gamma*h", c),
                    ("gamma/h", gam / h)):
        print(f"  {name:10s} in {iv.nstr(x, 14)}")
# checks used in the text (SU(2))
N = 2
b = iv.mpf(11 * N * N) / (12 * pi**2)
h = b * log2
gam = iv.mpf(51) / 121 * b
c = gam * h
assert b > iv.mpf('0.37151') and b < iv.mpf('0.37152')
assert h > iv.mpf('0.25751') and h < iv.mpf('0.25752')
assert gam > iv.mpf('0.15658') and gam < iv.mpf('0.15659')
assert c > iv.mpf('0.04032') and c < iv.mpf('0.04033')
print("SU(2) enclosures used in the text: b in [0.37151,0.37152], h in [0.25751,0.25752],")
print("  gamma in [0.15658,0.15659], c in [0.04032,0.04033]: all checked")
