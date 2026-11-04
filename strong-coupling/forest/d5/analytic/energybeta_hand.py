# Exact-rational check of the hand constants in patches/A02_energybeta.tex
from fractions import Fraction as F
def le(a,b,m): assert a<=b,(m,float(a),float(b)); print("ok",m,float(a),"<=",float(b))
def artanh2(x,m):  # 2*artanh(x): lower = partial sum to x^(2m+1), upper = + next term/(1-x^2)
    s=sum(2*x**(2*k+1)/(2*k+1) for k in range(m+1)); nxt=2*x**(2*m+3)/(2*m+3)
    return s, s+nxt/(1-x*x)
l2=artanh2(F(1,3),8); l32=artanh2(F(1,5),6); l54=artanh2(F(1,9),5)
le(l2[1],F(6931472,10**7),"log2<0.6931472")
le(l2[1]+l32[1],F(10986123,10**7),"log3<1.0986123")
le(F(16094379,10**7),2*l2[0]+l54[0],"log5>1.6094379")
le(F(31415926536,10**10),F(31416,10**4),"pi<3.1416 (pi=3.14159265358979...)")
le(F(31416,10**4)/3,F(10472,10**4),"pi/3<1.0472")
x=F(472,10**4); le(x-x*x/2+x**3/3,F(46122,10**6),"log1.0472<0.046122")
lpi=F(10986123,10**7)+F(46122,10**6); le(lpi,F(1144735,10**6),"logpi")
L2,L3,L5,LP=F(6931472,10**7),F(10986123,10**7),F(16094379,10**7),F(1144735,10**6)
K1=F(5,4)*(2*L2+L3/2+LP/4+F(3,4)); le(K1,F(371474,10**5),"K1"); le(F(371474,10**5),F(37148,10**4),"K1<=3.7148")
K2=F(5,24)*(F(7,2)*L2+L3+LP-F(3,2)*L5); le(K2,F(46984,10**5),"K2"); le(F(46984,10**5),F(4699,10**4),"K2<=0.4699")
c=F(5,8); K3=F(5,4)*F(49,48)*(c/6)/(1-c/6); assert K3==F(1225,8256); le(K3,F(1484,10**4),"K3 term"); le(F(1484,10**4),F(149,1000),"<=0.149")
# symbolic simplification check with floats
import math
lam=0.2
k1=1.25*(0.25*math.log(2/math.pi)-0.75*math.log(lam)-0.5*math.log(2/(3*math.pi))-0.75*math.log(5/8)+1.2*5/8)
print("K1 formula vs closed",k1,1.25*(2*math.log(2)+0.5*math.log(3)+0.25*math.log(math.pi)+0.75))
k2=5/24*(-math.log(2/(3*math.pi))-1.5*math.log(5/8)); print("K2",k2,5/24*(3.5*math.log(2)+math.log(3)+math.log(math.pi)-1.5*math.log(5)))
