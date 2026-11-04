"""Certified constants for Proposition prop:scalelosses (end of ssec:slabenergy).

All quantities are enclosed with mpmath interval arithmetic (mpmath.iv).

  b = 11/(3 pi^2)   one-loop coefficient of beta_W = 4/g^2 for SU(2)
  h = b log 2       change of the one-loop coupling per dyadic scale
  4/h               exponent of the rate (log 1/a_k)^(-4c/h) in the text,
                    for a loss c g_j^2 = 4c/beta_j per scale

The text uses: h in [0.25751, 0.25752], h in [0.2575118, 0.2575119],
4/h < 15.534.  The assertions below check these claims.
"""
from mpmath import iv
iv.dps = 40

b = iv.mpf(11) / (3 * iv.pi**2)
h = b * iv.log(2)
print("b   =", b)
print("h   =", h)
print("1/h =", 1 / h)
print("4/h =", 4 / h)
assert h.a > iv.mpf('0.25751').b and h.b < iv.mpf('0.25752').a
assert h.a > iv.mpf('0.2575118').b and h.b < iv.mpf('0.2575119').a
assert (4 / h).b < iv.mpf('15.534').a
print("all assertions passed")
