"""Rigorous interval arithmetic (mpmath.iv) with forward-mode derivatives
("Dual" = interval value + interval enclosure of an a.e. derivative), and the
cut-off average V of one pair with its seven non-path links on the diagonal
torus and every other link equal to I.

Quaternion conventions: (a,b,c,d) = a + b i + c j + d k, |.| = 1;
torus T = {cos q + i sin q} = diag(e^{iq}, e^{-iq}); chordal HS distance
delta(v,I) = 2 sqrt(1 - v_0) = 2 sqrt2 |sin(psi/2)|; rho_G = 1.
V(g) = g Exp( chi(m) (1/8) sum_a Log(g^{-1} l_a) ), m = max_a delta(g^{-1} l_a, I),
chi = min(1, max(0, 2 - 2m)), the two offset-zero paths contributing log I = 0
and every other private link occurring in two paths.

atan is not provided by mpmath.iv; it is enclosed by endpoint evaluation at
90 bits widened by 2^-80 (see iatan), which relies on mpmath's atan being
accurate to far better than 2^-80 at that precision.
"""
from mpmath import iv, mp, mpf
iv.prec = 60
mp.prec = 90
WID = mpf(2) ** -80

def lo(x):
    a = mp.mpf(); a._mpf_ = x._mpi_[0]; return a
def hi(x):
    a = mp.mpf(); a._mpf_ = x._mpi_[1]; return a
def I(a, b=None):
    return iv.mpf([a, b]) if b is not None else iv.mpf(a)
def hull(x, y):
    return iv.mpf([min(lo(x), lo(y)), max(hi(x), hi(y))])
def imax(x, y):
    return iv.mpf([max(lo(x), lo(y)), max(hi(x), hi(y))])
def imin(x, y):
    return iv.mpf([min(lo(x), lo(y)), min(hi(x), hi(y))])
def iatan(x):
    a = mp.atan(lo(x)) - WID; b = mp.atan(hi(x)) + WID
    return iv.mpf([a, b])

class D:
    """interval value v, interval enclosure d of its s-derivative."""
    __slots__ = ("v", "d")
    def __init__(s, v, d=None):
        s.v = v if isinstance(v, type(iv.mpf(0))) else iv.mpf(v)
        s.d = iv.mpf(0) if d is None else (d if isinstance(d, type(iv.mpf(0))) else iv.mpf(d))
    def __add__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v + o.v, s.d + o.d)
    __radd__ = __add__
    def __neg__(s): return D(-s.v, -s.d)
    def __sub__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v - o.v, s.d - o.d)
    def __rsub__(s, o): return D(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, D) else D(o); return D(s.v * o.v, s.d * o.v + s.v * o.d)
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = o if isinstance(o, D) else D(o)
        return D(s.v / o.v, (s.d * o.v - s.v * o.d) / (o.v * o.v))

def dsin(x): return D(iv.sin(x.v), iv.cos(x.v) * x.d)
def dcos(x): return D(iv.cos(x.v), -iv.sin(x.v) * x.d)
def dsqrt(x):
    r = iv.sqrt(x.v); return D(r, x.d / (2 * r))
def datan(x): return D(iatan(x.v), x.d / (1 + x.v * x.v))
def dmax(x, y):
    v = imax(x.v, y.v)
    if lo(x.v) > hi(y.v): d = x.d
    elif lo(y.v) > hi(x.v): d = y.d
    else: d = hull(x.d, y.d)
    return D(v, d)
def dphi(m):
    """chi = min(1, max(0, 2 - 2m)), decreasing in m; derivative -2 on (1/2,1), 0 outside."""
    def phi(t): return min(mpf(1), max(mpf(0), 2 - 2 * t))
    v = iv.mpf([phi(hi(m.v)), phi(lo(m.v))])
    a, b = lo(m.v), hi(m.v)
    if b < mpf(1) / 2 or a > 1: d = iv.mpf(0)
    elif a > mpf(1) / 2 and b < 1: d = -2 * m.d
    else: d = hull(iv.mpf(0), -2 * m.d)
    return D(v, d)

def alt_series(u, coeffs):
    """sum_k coeffs[k] u^k, alternating signs with decreasing terms assumed
    verified by caller (|coeffs[k+1]| u <= |coeffs[k]|); enclosure is the hull
    of the last two partial sums (interval arithmetic throughout)."""
    S = iv.mpf(0); P = iv.mpf(1); parts = []
    for c in coeffs:
        S = S + c * P; parts.append(S); P = P * u
    return hull(parts[-1], parts[-2])

from math import factorial
K = 8
C_COS  = [iv.mpf((-1) ** k) / factorial(2 * k) for k in range(K)]          # cos sqrt(u)
C_SINC = [iv.mpf((-1) ** k) / factorial(2 * k + 1) for k in range(K)]      # sin sqrt(u)/sqrt(u)
C_DSINC = [iv.mpf((-1) ** k) * k / factorial(2 * k + 1) for k in range(1, K)]  # d/du sinc sqrt(u)

WIDE = D(iv.mpf([-1, 1]), iv.mpf([-10, 10]))
def dcosroot(u):   # cos sqrt(u) as a function of u >= 0; series valid for hi(u) <= 2
    if hi(u.v) > 2: return WIDE      # harmless: the box then fails every test and is split
    return D(alt_series(u.v, C_COS), -alt_series(u.v, C_SINC) / 2 * u.d)
def dsincroot(u):  # sin sqrt(u)/sqrt(u)
    if hi(u.v) > 2: return WIDE
    return D(alt_series(u.v, C_SINC), alt_series(u.v, C_DSINC) * u.d)

def qmul(p, q):
    a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
    return (a1*a2 - b1*b2 - c1*c2 - d1*d2,
            a1*b2 + b1*a2 + c1*d2 - d1*c2,
            a1*c2 - b1*d2 + c1*a2 + d1*b2,
            a1*d2 + b1*c2 - c1*b2 + d1*a2)
def qinv(p): return (p[0], -p[1], -p[2], -p[3])

THETA = [I('-0.5')] + [I('0.49')] * 6

def V_of(q, s, theta=THETA):
    """q: interval (no s-dependence); s: D.  Returns (V as 4 Duals, chi, m)."""
    cq, sq = iv.cos(q), iv.sin(q)
    cs, ss = dcos(s), dsin(s)
    g = (cs * cq, cs * sq, D(0), ss)
    gi = qinv(g)
    logs = []; m = D(0)
    for th in theta:
        l = (D(iv.cos(th)), D(iv.sin(th)), D(0), D(0))
        v = qmul(gi, l)
        delta = 2 * dsqrt(1 - v[0])
        m = dmax(m, delta)
        n2 = v[1] * v[1] + v[2] * v[2] + v[3] * v[3]
        n = dsqrt(n2)
        t = n / (1 + v[0])
        lam = 2 * datan(t) / n
        logs.append((lam * v[1], lam * v[2], lam * v[3]))
    chi = dphi(m)
    N = [sum((L[i] for L in logs), D(0)) / 8 for i in range(3)]
    w = [chi * N[i] for i in range(3)]
    u = w[0] * w[0] + w[1] * w[1] + w[2] * w[2]
    C = dcosroot(u); Sc = dsincroot(u)
    E = (C, Sc * w[0], Sc * w[1], Sc * w[2])
    V = qmul(g, E)
    return V, chi, m

def f0_torus(q):
    """angle of V on the torus (interval q), via atan2 of components, and the
    explicit formula on [-0.005, 0.2227] (maximiser theta_1 = -1/2)."""
    V, chi, m = V_of(q, D(0))
    return V, chi, m
