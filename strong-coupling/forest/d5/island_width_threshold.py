"""Island threshold with the joining distance 2s+7 of the factorisation
theorem (two islands and one kernel) in place of 2s+1: only the number of
neighbours Delta changes, from (4s+3)^4-1 to (4s+15)^4-1.  The last block
is the radius at which the chart hypothesis kappa + 6*32 r <= 1e-8 holds
on the layer with kappa <= eta + 4r, r = 43 rho, eta = 92 rho."""
from mpmath import mp, mpf, log, findroot, e, pi
mp.dps=30
M=mpf(11)**4; H=mpf(6)**4; C0=mpf(1)/2+mpf(2)/3*log(24*pi**3)
def beta_star(rho,s,join):
    D=(2*join+1)**4-1
    logqstar=-H*(log(2)+2*log(D))-log(e*M*H)
    f=lambda b: 6*(log(b)+C0)-b*rho**2/16-logqstar
    b=findroot(f, 16*(-logqstar)/rho**2*mpf('1.02'))
    return b, logqstar, D
for rho,s in ((mpf('1e-8'),1),(mpf('1e-8'),5),(mpf('1e-9'),1)):
    for join in (2*s+1, 2*s+7):
        b,lq,D=beta_star(rho,s,join)
        print(f"rho={rho} s={s} join={join} Delta={int(D)} log q*={float(lq):.1f} beta*={float(b):.4e}")
rho_c=mpf('1e-8')/(92+196*43)
print(f"rho_c={float(rho_c):.4e}")
for s in (1,5):
    b,lq,D=beta_star(rho_c,s,2*s+7)
    print(f"rho=rho_c s={s} join={2*s+7} log q*={float(lq):.1f} beta*={float(b):.4e}")
