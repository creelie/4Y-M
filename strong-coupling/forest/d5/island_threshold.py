from mpmath import mp, mpf, log, exp, findroot, e, pi
mp.dps=30
M=mpf(11)**4; H=mpf(6)**4; C0=mpf(1)/2+mpf(2)/3*log(24*pi**3)
def Delta(s): return (4*s+3)**4-1
def beta_star(rho,s,Gamma=1):
    D=Delta(s)
    logqstar=-H*(log(2*Gamma)+2*log(D))-log(e*M*H)
    f=lambda b: 6*(log(b)+C0)-b*rho**2/16-logqstar
    b=findroot(f, 16*(-logqstar)/rho**2)
    return b, logqstar, D
for s in (1,5):
  for rho in (mpf('1e-8'), mpf('1e-9')):
    b,lq,D=beta_star(rho,s)
    print(f"s={s} rho={rho} Delta={int(D)} log q*={float(lq):.1f} beta*={float(b):.4e}")
b,lq,D=beta_star(mpf('1e-8'),5,mpf('1e30')); print("Gamma=1e30, s=5, rho=1e-8: beta*=%.4e"%float(b))
print("C0=",C0, " 43 vs 21+7pi:",21+7*pi, " 92 vs 4+28pi:",4+28*pi)
# check: at beta*, A^{1/H} Delta^2 == 1/2
rho=mpf('1e-8'); s=5; b,lq,D=beta_star(rho,s); q=exp(6*(log(b)+C0)-b*rho**2/16); A=e*M*H*q; print("check", A**(1/H)*D**2)
