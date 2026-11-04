"""Exact-rational re-check of the decimals in the proof of Theorem su3blocks (hand constants).
The text does not depend on this script."""
from fractions import Fraction as F
from math import factorial
m=[F(1),F(0),F(1,2),F(1,4),F(3,4),F(15,16),F(65,32),F(231,64),F(973,128)]
def explow(x): return 1-x+x**2/2-x**3/6
def chain(bw, eps, q, Lam, tau, w, s, C, K, a, gapdec):
    bw=F(bw); r=2*bw; Nb=bw/3; g=F(283,128)-r/2
    eps,q,Lam,tau,w,s,C,K,a=map(F,(eps,q,Lam,tau,w,s,C,K,a))
    assert eps < q*q
    assert q/g <= Lam, float(q/g)          # Lambda(r) <= sqrt(eps)/gamma, and sqrt(tau/gamma) = sqrt(eps)/gamma
    assert eps/g <= tau, float(eps/g)
    assert F(3,2)*Nb**2*tau <= w, float(F(3,2)*Nb**2*tau)
    ew=1/(1-w)                              # e^w <= 1/(1-w)
    assert Lam*(1+(ew-1)*ew) <= s, float(Lam*(1+(ew-1)*ew))
    # M = 2 sqrt3 tau, so sqrt3 * Nb * M * e^w / gamma = 6 Nb tau e^w / gamma
    Cv=5*Nb**2*(tau*s+6*Nb*tau*ew/g); assert Cv<=C, float(Cv)
    Kv=(1+Nb*Lam)*s/(1-C); assert Kv<=K, float(Kv)
    av=10*Nb*K; assert av<=a, float(av)
    y=1-a; assert (y+y*y/2)/2 >= F(gapdec)   # -1/2 log alpha >= (y+y^2/2)/2
    print('beta_W',bw,'alpha_3 <=',a,' gap >=',gapdec)
# beta_W = 0.515: eps(1.03) <= 0.667483 < 0.817^2 from the proof of Theorem su3
chain('0.515','0.667483','0.817','0.48174','0.39358','0.017398','0.49043','0.064287','0.56747','0.9742','0.013')
# beta_W = 0.52: Z_*(1.04)
r=F(104,100)
terms_txt=['1','0','0.2704','0.0468694','0.0365581','0.0095052','0.0035697','0.0009424','0.0002581']
terms=[m[n]*r**n/factorial(n) for n in range(9)]
for t,tt in zip(terms,terms_txt): assert t<=F(tt),(float(t),tt)
S=sum(F(t) for t in terms_txt); assert S<=F('1.3681029'), float(S)
x=3*r; tail=m[8]/3**8*x**9/factorial(9)/(1-x/10); assert tail<F('0.000131'), float(tail)
assert F('1.3681029')+F('0.000131')<=F('1.368234')
assert explow(r/5)>=F('0.812132'), float(explow(r/5))
ratio=F(5,9)*F('0.812132')/F('1.368234'); assert ratio>=F('0.329756'), float(ratio)
assert 1-F('0.329756')<=F('0.670244')
chain('0.52','0.670244','0.81869','0.48417','0.39638','0.017864','0.49314','0.066654','0.5727','0.9927','0.0036')
print('all decimals of the proof of Theorem su3blocks confirmed')
