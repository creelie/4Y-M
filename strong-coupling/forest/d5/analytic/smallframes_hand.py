# Check of the hand bound 28.04 in patches/A03_smallframessharp.tex
from fractions import Fraction as F
import math
def le(a,b,m): assert a<=b,(m,float(a),float(b)); print("ok",m,float(a),"<=",float(b))
# float reproduction of the paper value
psi=math.asin(1/(4*math.sqrt(2)))/3; kc=2*math.sqrt(2)*math.sin(psi/2)
ph=lambda a:a/2+0.09*a*a; ps=lambda a:a/2*math.exp(a)
pG=lambda p:18*p+(1+math.sqrt(3)/4)*(ps(6*p)*(1+ph(12*p))+ph(12*p))+math.sqrt(2880)/16*p
print("float value pG(psi_c)/kappa_c =",pG(psi)/kc,"psi_c",psi,"kappa_c",kc)
# exact hand chain
y=F(1767767,10**7); le(F(1,32),y*y,"1/(4sqrt2)<=0.1767767")
asn=y+y**3/(6*(1-y*y)); P=F(59243,10**6); le(asn/3,P,"psi_c<=0.059243")
le(F(1414006,10**6),F(14142135,10**7)*(1-P*P/24),"kappa_c/psi_c>=1.414006"); le(F(14142135,10**7)**2,2,"sqrt2>=1.4142135")
a=6*P; b=12*P
def exp_up(x,n=8):
    s=sum(x**k/math.factorial(k) for k in range(n+1)); t=x**(n+1)/math.factorial(n+1)/(1-x/(n+2)); return s+t
E=F(142686,10**5); le(exp_up(a),E,"e^{6p}<=1.42686")
phh=b/2+F(9,100)*b*b; le(phh,F(400945,10**6),"phihat(12p)")
psh=a/2*E; le(psh,F(253595,10**6),"psihat(6p)")
inner=F(253595,10**6)*(1+F(400945,10**6))+F(400945,10**6); le(inner,F(756218,10**6),"inner")
s3=F(17320509,10**7); le(3,s3*s3,"sqrt3<=1.7320509"); c=1+s3/4; le(c,F(14330128,10**7),"1+sqrt3/4")
r=F(5366564,10**5); le(2880,r*r,"sqrt2880<=53.66564")
G=18+F(14330128,10**7)*F(756218,10**6)/P+r/16; le(G,F(396461,10**4),"G(psi_c)<=39.6461")
le(F(396461,10**4)/F(1414006,10**6),F(2804,100),"ratio<28.04")
print("--- kappa_c enclosure (A03, second hunk)")
r2lo,r2hi=F(14142135,10**7),F(14142136,10**7); assert r2lo**2<=2<=r2hi**2
ylo,yhi=1/(4*r2hi),1/(4*r2lo)
c=[F(1),F(1,6),F(3,40),F(15,336)]
def asin_lo(y): return sum(c[k]*y**(2*k+1) for k in range(4))
def asin_hi(y): return sum(c[k]*y**(2*k+1) for k in range(3))+F(1,6)*y**7/(1-y*y)
plo,phi=asin_lo(ylo)/3,asin_hi(yhi)/3
le(F(5923686,10**8),plo,"psi_c>=0.05923686"); le(phi,F(5923710,10**8),"psi_c<=0.05923710")
klo=2*r2lo*(plo/2-(plo/2)**3/6); khi=2*r2hi*(phi/2-(phi/2)**3/6+(phi/2)**5/120)
le(F(8376130,10**8),klo,"kappa_c>=0.0837613"); le(khi,F(837617,10**7),"kappa_c<=0.0837617")
le(8/klo**2,F(114026,100),"8/kappa_c^2<1140.26"); le(1441*F(114026,100),F(16432,10**4)*10**6,"C_*")
