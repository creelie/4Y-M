# Exact-rational check of patches/A07_constrainedscales.tex (exp bounded above by 1+x+x^2 for tiny x)
from fractions import Fraction as F
import math
def lg(x): return math.log10(x.numerator)-math.log10(x.denominator) if x>0 else float('-inf')
def le(a,b,m):
    a,b=F(a),F(b); assert a<=b,(m,lg(a),lg(b)); print("ok",m,"log10:",round(lg(a),6),"<=",round(lg(b),6))
C=F(10)**26; eta=F(1,10**36); l1=F(1,10**215); s=24576; e=lambda x:1+x+x*x
le((2/eta)**5,F(32)*10**180,"(2/eta)^5")
le(24*C*l1*e(l1)*(1+384*e(eta/2)*F(32)*10**180),F(295,10**6),"first sum")
le(F(17,8)*s*(e(l1)-1),F(1,10**209),"other sums"); le(F(3,10**4),F(1,210),"<1/(2M0)")
x=(eta-l1)/2; le(2/(eta-l1),F(20001,10**4)*10**36,"2/(eta-l1)")
coth=1/x+x/3; le(24*C*e(l1)*coth**4,4*F(10)**172,"f bound")
le((210*(4*F(10)**172+2000))**2,F(10)**350,"first bound"); le(e(l1)*s/16,2000,"s/16")
def Z(x): return 4*e(x)*(1+384*e(x)/x**4+16*e(x)/x**2)
le(e(2*l1),F(1001,1000),"e^{2l1}"); le(Z(l1),4*F(1001,1000)*(1+F(3844,10)*F(10)**860+F(161,10)*F(10)**430),"Z(l1) form"); le(Z(l1),F(10)**864,"Z(l1)")
le(Z(2*l1),4*F(1001,1000)*(1+F(241,10)*F(10)**860+F(41,10)*F(10)**430),"Z(2l1) form"); le(Z(2*l1),F(97,10)*F(10)**861,"Z(2l1)")
le(F(97,10)*F(10)**861,(F(985,100)*F(10)**430)**2,"sqrt"); le(210*s,F(517,100)*10**6,"210s")
le(F(517,100)*10**6*F(985,100)*F(10)**430,F(10)**438,"L_w")
q=l1/2; le(1/q,2*F(10)**215,"1/q")
le(F(10)**607*128*F(1001,1000)/q,F(257,100)*F(10)**824,"first term coef")
le(32*F(257,100)*F(10)**824,F(83,10)*F(10)**825,"n threshold"); le(F(83,10)*F(10)**825,F(10)**827,"<=1e827")
le(F(10)**607*128*F(1001,1000)*48/q**4,F(99,10)*F(10)**1471,"second term"); le(32*F(99,10)*F(10)**1471,F(32,10)*F(10)**1473,"n^4 threshold")
