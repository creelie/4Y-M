"""Exact rational enclosures used by the checks (no floating point)."""
from fractions import Fraction as F

def log_gt1(y, K=60):
    """Enclosure (lo, hi) of log y for rational y > 1 by the artanh series
    log y = 2 sum_{k>=0} u^{2k+1}/(2k+1), u=(y-1)/(y+1), tail bounded
    geometrically by 2 u^{2K+1}/((2K+1)(1-u^2))."""
    y = F(y); assert y > 1
    u = (y - 1) / (y + 1)
    s = F(0); p = u
    for k in range(K):
        s += 2 * p / (2 * k + 1)
        p *= u * u
    tail = 2 * p / ((2 * K + 1) * (1 - u * u))
    return s, s + tail

def log_enc(y, K=60):
    y = F(y)
    if y == 1: return F(0), F(0)
    if y > 1: return log_gt1(y, K)
    lo, hi = log_gt1(1 / y, K)
    return -hi, -lo

def exp_upper(x, n=30):
    """Upper bound of e^x for rational 0 <= x <= 1: partial sum plus
    tail bounded by x^n/n! * 1/(1-x/(n+1))."""
    x = F(x); s = F(0); t = F(1)
    for k in range(n):
        s += t; t = t * x / (k + 1)
    return s + t / (1 - x / (n + 1))

def exp_lower(x, n=30):
    x = F(x); s = F(0); t = F(1)
    for k in range(n):
        s += t; t = t * x / (k + 1)
    return s

def log_any(y, K=80):
    """Enclosure of log y, y>0 rational: y = 2^k z with 1 <= z < 2."""
    y = F(y); k = 0
    while y >= 2: y /= 2; k += 1
    while y < 1: y *= 2; k -= 1
    l2lo, l2hi = log_gt1(F(2), K)
    zlo, zhi = log_enc(y, K)
    if k >= 0: return k * l2lo + zlo, k * l2hi + zhi
    return k * l2hi + zlo, k * l2lo + zhi

class I:
    """Closed interval with Fraction end points."""
    def __init__(s, lo, hi=None):
        s.lo = F(lo); s.hi = F(hi if hi is not None else lo)
        assert s.lo <= s.hi
    @staticmethod
    def c(x): return x if isinstance(x, I) else I(x)
    def __add__(s, o): o = I.c(o); return I(s.lo + o.lo, s.hi + o.hi)
    __radd__ = __add__
    def __neg__(s): return I(-s.hi, -s.lo)
    def __sub__(s, o): return s + (-I.c(o))
    def __rsub__(s, o): return I.c(o) - s
    def __mul__(s, o):
        o = I.c(o); p = [s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi]
        return I(min(p), max(p))
    __rmul__ = __mul__
    def inv(s):
        assert s.lo > 0 or s.hi < 0
        return I(1 / s.hi, 1 / s.lo)
    def __truediv__(s, o): return s * I.c(o).inv()
    def __rtruediv__(s, o): return I.c(o) * s.inv()
    def __pow__(s, n):
        r = I(1)
        for _ in range(n): r = r * s
        return r
    def __repr__(s): return "[%.12g, %.12g]" % (float(s.lo), float(s.hi))
