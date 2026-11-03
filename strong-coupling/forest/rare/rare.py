"""Explicit constants for the chessboard bound on bad plaquettes (C2, RARE.md).

log z <= -beta*rho^2/16 + 6 - 4*log h(t),  t = (8 beta)^(-1/2),
h(t) = Haar mass of {angle <= t} in SU(2) = (t - sin t cos t)/pi.
Bad base sites: P(all of S bad) <= (6z)^|S|.  Clusters (adjacency of degree D)
have exponential tails once q = e*D*6*z < 1.
"""
from mpmath import mp, mpf, sin, cos, pi, log, e, findroot, sqrt
mp.dps = 120

def logz(beta, rho):
    t = 1/sqrt(8*beta)
    h = (t - sin(t)*cos(t))/pi
    return -beta*rho**2/16 + 6 - 4*log(h)

def logq(beta, rho, D):
    return logz(beta, rho) + log(6*e*D)

def threshold(rhofun, D, lo=1e2, hi=1e40):
    # smallest beta (on a log grid, then bisection) with log q < 0
    import math
    a, b = mpf(lo), mpf(hi)
    assert logq(b, rhofun(b), D) < 0
    for _ in range(200):
        m = sqrt(a*b)
        if logq(m, rhofun(m), D) < 0: b = m
        else: a = m
    return b

print("adjacency D: 8 = nearest neighbour, 80 = l-infinity (*-connected)")
for D in (8, 80):
    print(f"\nD = {D}")
    for name, f in [("rho = beta^(-1/4)", lambda b: b**mpf(-0.25)),
                    ("rho = 0.1", lambda b: mpf('0.1')),
                    ("rho = 0.01", lambda b: mpf('0.01')),
                    ("rho = 3.25e-12 (Cor. localcalm, m=1)", lambda b: mpf(1e-8)/3072)]:
        bt = threshold(f, D)
        print(f"  {name:40s} q<1 for beta_W >= {mp.nstr(bt,4)}")
    for b in (mpf(10)**7, mpf(10)**8, mpf(10)**10):
        r = b**mpf(-0.25)
        print(f"  beta_W={mp.nstr(b,3)}, rho=beta^-1/4: log z <= {mp.nstr(logz(b,r),6)}, q = {mp.nstr(mp.exp(logq(b,r,D)),4)}")
