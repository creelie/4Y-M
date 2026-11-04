"""Certified island thresholds (thm:islandthreshold, thm:islandscales, thm:islandswide) and k_*.
Interval arithmetic: mpmath.iv.  Usage: python3 island_radius.py
Every root beta* is certified by g(lo)>0>g(hi) with lo>96/rho0^2 (g decreasing there)."""
from mpmath import mp, iv, mpf, findroot, log as mlog
mp.dps = 40; iv.dps = 40
I = lambda x: iv.mpf(x)
C0 = I(1)/2 + I(2)/3*iv.log(24*iv.pi**3)
C1 = C0 + iv.log(I(12800))/6
def Delta(s): return I((4*s+3)**4 - 1)
def logq(M, H, s):
    return -I(H)*(iv.log(I(2)) + 2*iv.log(Delta(s))) - (1 + iv.log(I(M)) + iv.log(I(H)))
def g(beta, rho0, lq):   # f(beta)-log q*
    b = I(beta); return 6*(iv.log(b)+C0) - b*I(rho0)**2/16 - lq
def root(rho0, lq):
    lqm = mpf(lq.b); r0 = mpf(rho0)
    fm = lambda b: 6*(mlog(b)+mpf(C0.b)) - b*r0**2/16 - lqm
    a, c = 96/r0**2, 64*abs(lqm)/r0**2      # fm(a)>0>fm(c); fm decreasing on [a,c]
    assert fm(a) > 0 and fm(c) < 0
    for _ in range(200):
        m = (a*c)**0.5
        if fm(m) > 0: a = m
        else: c = m
    b = (a*c)**0.5
    lo, hi = b*(1-mpf(10)**-12), b*(1+mpf(10)**-12)
    assert lo > 96/r0**2
    assert g(lo, I(rho0), lq).a > 0 and g(hi, I(rho0), lq).b < 0, "not certified"
    return lo, hi
def Hj(j): return 8**4 if j <= 2 else 9**4
def Mj(j): return (17*2**j - 4)**4
def betastar_scale(rho, j, s):
    rho0 = I(rho)/(I('1.01')*I('1.023')*I(4)**j)
    rho0m = mpf(rho0.a)   # lower end of rho0 -> larger root: upper bound side
    lq = logq(Mj(j), Hj(j), s)
    # certify with both ends of rho0 interval: take root for rho0.a (largest beta) and rho0.b
    lo1, hi1 = root(mpf(rho0.b), lq); lo2, hi2 = root(mpf(rho0.a), lq)
    return min(lo1, lo2), max(hi1, hi2)
def betastar_one(rho, s):
    return root(mpf(rho), logq(11**4, 6**4, s))
def kstar(beta):
    b = I(beta); rho = iv.sqrt(96*(iv.log(b)+C1)/b)
    t = iv.log(I('1e-4')/(I('1.01')*rho))/iv.log(I(4))
    import math
    fa, fb = math.floor(mpf(t.a)), math.floor(mpf(t.b))
    assert fa == fb, ("k* ambiguous", beta)
    return fa
def run(rho, s, jmax_scan, logb_lo, logb_hi, npts, label):
    print(f"=== {label}: rho={rho}, s={s}")
    b0 = betastar_one(rho, s); print(f"one scale beta*({rho},{s}) in [{mp.nstr(b0[0],6)}, {mp.nstr(b0[1],6)}]")
    bs = [b0]
    worst = mpf(0)
    for j in range(1, jmax_scan+1):
        lo, hi = betastar_scale(rho, j, s); bs.append((lo, hi))
        c = hi/mpf(16)**j; worst = max(worst, c)
        if j <= 4 or j == jmax_scan: print(f" j={j}: beta*_j in [{mp.nstr(lo,6)}, {mp.nstr(hi,6)}], 16^-j beta*_j <= {mp.nstr(c,6)}")
        assert lo/bs[-2][1] >= 16 if j >= 2 else True
    print(f" max_j<= {jmax_scan} of 16^-j beta*_j <= {mp.nstr(worst,6)}")
    hist = {}
    for i in range(npts):
        lb = mpf(logb_lo) + (mpf(logb_hi)-logb_lo)*i/(npts-1)
        beta = mpf(10)**lb
        jm = -1
        for j,(lo,hi) in enumerate(bs):
            if beta >= hi: jm = j
            elif beta >= lo: raise Exception("ambiguous j_max")
            else: break
        assert jm < jmax_scan, "scan more scales"
        if jm < 0: continue
        d = kstar(beta) - jm
        hist[d] = hist.get(d, 0) + 1
    print(f" k*-j_max over {npts} couplings 10^{logb_lo}..10^{logb_hi} (skipping beta<beta*_0):", dict(sorted(hist.items())))
    return bs
if __name__ == "__main__":
    # reproduce the paper
    for s in (1, 5):
        b = betastar_one(mpf('1e-8'), s); print(f"paper check beta*(1e-8,{s}) in [{mp.nstr(b[0],6)},{mp.nstr(b[1],6)}]")
    print("|log q*| s=1:", mp.nstr(-logq(11**4,6**4,1).b,8), " s=5:", mp.nstr(-logq(11**4,6**4,5).b,8))
    run(mpf("1e-8"), 5, 20, 22, 42, 2001, "paper")
    for rho in ('1e-4', '1.023e-4'):
        run(mpf(rho), 5, 30, 22, 42, 2001, "chart radius")
        run(mpf(rho), 5, 30, 14, 42, 2801, "chart radius, wide range")
    for rho in ('1e-4', '0.0837'):
        for s in (1, 5):
            b = betastar_one(mpf(rho), s); print(f"one scale beta*({rho},{s}) in [{mp.nstr(b[0],6)},{mp.nstr(b[1],6)}]")
