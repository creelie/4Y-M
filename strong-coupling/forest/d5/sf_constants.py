"""Certified constants for smallfield.tex (mpmath.iv, 60 digits).

c_2(lam,eta) = 24 log 1.957 + 73 log 14 + 109.5 (1 + log(15979 (1+1/eta) lam))
a_2(lam,kappa) = lam kappa^2 (1 - 216 (1+eta) (kappa'/kappa)^2) - 8 log M - c_2(lam,eta)
with kappa' = kappa/30, eta = 1, M = 1 + 2 (1 - 7 Lambda/8)^{-3} < 23621 (Lemma lem:averagedensity).
"""
from mpmath import iv, mp
iv.dps = 60
M = iv.mpf(23621)                     # upper bound for M
logM = iv.log(M)
def c2(lam, eta):
    return 24*iv.log(iv.mpf('1.957')) + 73*iv.log(iv.mpf(14)) + \
        iv.mpf('109.5')*(1 + iv.log(15979*(1 + 1/iv.mpf(eta))*lam))
def a2(lam, kap, ratio=30, eta=1):
    fac = 1 - 216*(1 + iv.mpf(eta))/iv.mpf(ratio)**2
    return lam*kap**2*fac - 8*logM - c2(lam, eta)
# the rho of the proof: rho^2 = 1/(15978.08 (1+1/eta) lam); 2*216*8100/219
print("2*216*8100/219 <", repr(float(((iv.mpf(2*216*8100)/219)).b)), "<= 15979")
assert (iv.mpf(2*216*8100)/219).b <= 15979
print("8 log M <", repr(float(((8*logM)).b)))
print("c_2(10^12,1) <", repr(float((c2(iv.mpf(10)**12, 1)).b)))
# threshold lam* with a_2 >= 9, by bisection on the left end certified
def threshold(kap):
    kap = iv.mpf(kap)
    lo, hi = iv.mpf(100), iv.mpf(10)**40
    # a_2 increasing for lam > 109.5/(0.52 kap^2); start lo there
    lo = iv.mpf(max(100, float((iv.mpf('109.5')/(iv.mpf('0.52')*kap**2)).b)))
    assert a2(lo, kap).b < 9
    assert a2(hi, kap).a > 9
    for _ in range(400):
        mid = iv.mpf(((lo + hi)/2).mid)
        if a2(mid, kap).a > 9: hi = mid
        else: lo = mid
        if (hi - lo).b / hi.a < 1e-9: break
    return hi
for kap in ['0.08', '1e-2', '1e-4', '1e-8']:
    lam = threshold(kap)
    k = iv.mpf(kap)
    rho = 1/iv.sqrt(iv.mpf('15978.08')*2*lam)
    # rho <= 9.32e-4 and rho <= (kappa-kappa')/90
    assert rho.b <= 9.32e-4 and rho.b <= (k*29/30/90).a
    assert a2(lam, k).a > 9
    print("kappa=%s: a_2 >= 9 for lam >= %.10g;  beta** = 4 lam <= %.7g;  lam kappa^2 <= %.7g"
          % (kap, float(lam.b), float((4*lam).b), float((lam*k**2).b)))
# Z_emptyset lower bound uses c_1 of thm:lfaction with rho^2 = 1/(11520(1+1/eta)lam); rho<=(kappa-kappa')/90
# holds when lam (kappa-kappa')^2 >= 8100/11520
print("8100/11520 =", 8100/11520)
