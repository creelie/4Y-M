"""Island threshold with the joining distance 2s+7 of the factorisation
theorem (two islands and one kernel) in place of 2s+1: only the number of
neighbours Delta changes, from (4s+3)^4-1 to (4s+15)^4-1.  The last block
is the radius at which the calmness hypothesis of the layer theorem holds,
eta + 49 r <= 6.0e-13 with r = 43 rho, eta = 92 rho, at the joining
distance 2s+8 that the two-island corollary needs."""
from mpmath import mp, mpf, log, findroot, e, pi
mp.dps=30
M=mpf(11)**4; H=mpf(6)**4; C0=mpf(1)/2+mpf(2)/3*log(24*pi**3)
def beta_star(rho,s,join):
    D=(2*join+1)**4-1
    logqstar=-H*(log(2)+2*log(D))-log(e*M*H)
    f=lambda b: 6*(log(b)+C0)-b*rho**2/16-logqstar
    lo,hi=16*(-logqstar)/rho**2,16*(-logqstar)/rho**2*mpf('1.1')
    assert f(lo)>0>f(hi)
    for _ in range(200):
        m=(lo+hi)/2
        if f(m)>0: lo=m
        else: hi=m
    b=lo
    return b, logqstar, D
for rho,s in ((mpf('1e-8'),1),(mpf('1e-8'),5),(mpf('1e-9'),1)):
    for join in (2*s+1, 2*s+7):
        b,lq,D=beta_star(rho,s,join)
        print(f"rho={rho} s={s} join={join} Delta={int(D)} log q*={float(lq):.1f} beta*={float(b):.4e}")
rho_c=mpf('6.0e-13')/(92+49*43)
print(f"rho_c={float(rho_c):.4e}")
for s in (1,5):
    b,lq,D=beta_star(rho_c,s,2*s+8)
    print(f"rho=rho_c s={s} join={2*s+8} log q*={float(lq):.1f} beta*={float(b):.4e}")
