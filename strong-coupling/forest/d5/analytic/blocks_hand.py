# Exact-rational check of the hand evaluation in patches/A01_blocks.tex.
from fractions import Fraction as F
def le(a,b,msg):
    assert a<=b,(msg,float(a),float(b)); print("ok",msg,float(a),"<=",float(b))
k1=F(2088,1000); k2=F(174,100)
# sqrt bounds by squaring
le(k1**2+9, F(36551,10000)**2, "sqrt(k1^2+9)<=3.6551")
umax=1+F(36551,10000)
lam=F(15435,100000)
le(lam,(8-umax)/umax**2,"lambda*>=0.15435 (8-u)/u^2 at u<=4.6551, decreasing in u")
le(F(3468,1000)**2, k2**2+9, "sqrt(k2^2+9)>=3.468")
Ag=F(38944,100000); le(k2/(1+F(3468,1000)),Ag,"A<=0.38944")
le(F(32572,10000)**2,k1**2+F(25,4),"sqrt(k1^2+25/4)>=3.2572")
le(k1**2+F(25,4),F(32573,10000)**2,"sqrt(k1^2+25/4)<=3.2573")
vmin=F(3,2)+F(32572,10000); vmax=F(3,2)+F(32573,10000)
Gup=(5*vmax+8)*k1/((F(32572,10000)+F(5,2))*vmin**3)
M=F(1071,10000); le(Gup,M,"G(k1)<=0.1071"); le(k1/24,M,"k1/24<=M")
def chain(bw,show=True):
    b=F(bw)/4
    w=8*b*b*(F(1,4)-lam)
    ew=1+w+w*w      # e^w<=1+w+w^2 for 0<=w<=1
    sq=F(288676,1000000); assert sq*sq>=F(1,12)
    s=F(1,4)+(ew-1)*ew*sq
    C=4*b*s*Ag+F(320,3)*b**3*M*ew
    K=(1+b)*s/(1-C)
    a=40*b*K
    print(bw, "w",float(w),"s",float(s),"C",float(C),"K",float(K),"alpha",float(a))
    return w,s,C,K,a
for bw in [F(1,3),F(345,1000),F(348,1000)]:
    chain(bw)
print("--- displayed chain at 0.348 (proof of thm:blocks)")
b=F(87,1000)
le(8*b*b,F(60552,1000000),"8b^2"); le(F(1,4)-lam,F(9565,100000),"tau-lam")
w=F(60552,1000000)*F(9565,100000); le(w,F(57918,10**7),"w")
W_=F(57918,10**7); le(1+W_+W_**2,F(10058254,10**7),"e^w")
le((W_+W_**2)*(1+W_+W_**2),F(58594,10**7),"(e^w-1)e^w")
sq=F(288676,10**6); le(F(1,12),sq*sq,"sqrt(1/12)<=0.288676")
le(F(1,4)+F(58594,10**7)*sq,F(2516915,10**7),"s")
s=F(2516915,10**7)
t1=F(348,1000)*s*Ag; le(t1,F(341107,10**7),"C term1")
t2=F(320,3)*b**3*M*F(10058254,10**7); le(t2,F(75670,10**7),"C term2")
le(F(341107+75670,10**7),F(416777,10**7),"C")
K=F(1087,1000)*s/(1-F(416777,10**7)); le(K,F(285490,10**6),"K")
le(1-F(416777,10**7),F(9583223,10**7),"1-C"); le(F(9583223,10**7),1-F(416777,10**7),"1-C eq")
le(F(348,100)*F(285490,10**6),F(99351,10**5),"alpha")
# m at 1/3: alpha<=0.94554<=0.9456 ; -log(1-y)>=y+y^2/2
y=1-F(9456,10000); le(F(279,10000),(y+y*y/2)/2,"m>=0.0279")
print("--- prop:blockresponse barriers")
import math
def sqrt_up(x,d=10**7):
    r=F(math.isqrt(int(x*d*d))+1,d); assert r*r>=x; return r
def sqrt_dn(x,d=10**7):
    r=F(math.isqrt(int(x*d*d)),d); assert r*r<=x; return r
k=F(351,1000)
Alo=k/(F(3,2)+sqrt_up(k*k+F(25,4))); Aup=k/(1+sqrt_dn(k*k+9))
print("A(0.351) in",float(Alo),float(Aup))
le(F(2718,10000),(1+Alo)/4,"cov>=0.2718"); le((1+Aup)/4,F(27183,100000),"cov<=0.27183")
k=F(3667,10000); Alo=k/(F(3,2)+sqrt_up(k*k+F(25,4)))
le(1,10*k*(1+Alo)/4,"10 beta (1+A)/4 >1 at 0.3667")
print("--- stated bounds at 0.345 and 1/3")
for bw,vals in [(F(345,1000),(F(56924,10**7),F(2516621,10**7),F(41185,10**6),F(285110,10**6),F(98363,10**5))),
                (F(1,3),(F(53139,10**7),F(2515504,10**7),F(393011,10**7),F(2836611,10**7),F(94554,10**5)))]:
    got=chain(bw)
    for g,v,n in zip(got,vals,"w s C K alpha".split()): le(g,v,f"{float(bw)} {n}")
