"""Certified constants (mpmath interval arithmetic) for the block curved.tex."""
from mpmath import iv, mpf
iv.dps=40
pi=iv.pi
def up(x): return x.b
def lo(x): return x.a
print("Lemma lem:smallcurv, kappa (n+2)^2 <= 1/100:")
r432=iv.sqrt(432)
# 2 delta/s_n <= 3 n (n+1) kappa <= 3/100
th=iv.mpf('0.03')
# rho_0 <= C_n^{1/2}/(1-th), C_n^{1/2} <= (n+2)/2
# gamma = sqrt432 * delta * rho_0 <= sqrt432*3(n+1)kappa*(n+2)/(2(1-th)) <= 3 sqrt432/(2(1-th)) * 1/100
gam=3*r432/(2*(1-th))/100
print("  sqrt(432) =",r432)
print("  gamma <=",up(gam))
ratio=1/((1-th)**2*(1-gam)**2)
print("  varrho^2/C_n <=",up(ratio), " (claim 2.31)")
assert up(ratio)<=2.31
print("  2.31*kappa*C_n <= 2.31/4*kappa*(n+2)^2 =",up(iv.mpf('2.31')/4),"kappa (n+2)^2 (claim 0.58)")
assert up(iv.mpf('2.31')/4)<=0.58
print("  plaquette sum: L_kappa+0.58 kappa (n+2)^2 M_kappa with M=24r: coefficient",up(iv.mpf('0.58')*24),"(claim 14)")
assert up(iv.mpf('0.58')*24)<=14
print("Proposition prop:curvedsharp (i): s^2<=1/100, n>=9, t=s^2(n+1)^2<=1")
s2=iv.mpf('0.01'); c2=1-s2
num=c2*(16/pi**2-iv.mpf(20)/100)-4*s2
den=3*pi**2/2+(2+16/pi**2)
print("  numerator coefficient >=",lo(num))
print("  denominator <=",up(den))
k=num/den
print("  ratio/t >=",lo(k)," ratio/(kappa (n+1)^2) >=",lo(k/2)," (claim 0.037)")
assert lo(k/2)>=0.037
print("Proposition prop:curvedsharp (ii): s^2<=1/100, n>=9")
coef=2-32*c2/pi**2+iv.mpf(24)/100
print("  coefficient of t <=",up(coef))
thr=(3*pi**2/2)/(-coef)
print("  t threshold <",up(thr),"  kappa(n+1)^2 = 2t threshold <",up(2*thr)," (claim 31)")
assert up(2*thr)<31
print("  3 pi^2/2 =",3*pi**2/2)
print("Corollary checks: 16/pi^2 =",16/pi**2)
print("R_m: L = eps(m^2+((m+1)^2+(m-1)^2)/2) = eps(2m^2+1) (exact)")
print("Haar term SU(2): (1/12)Tr(ad_x^2) = -(1/12)*4|x|^2 = -|x|^2/3 (exact)")
