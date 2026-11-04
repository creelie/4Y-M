# Check of patches/A04_islandmass.tex
from fractions import Fraction as F
def le(a,b,m): assert a<=b,(m,float(a),float(b)); print("ok",m,float(a),"<=",float(b))
k=F(5,10**4); d=k-k*k/2
sig=1+384/d**4+16/d**2; le(sig,F(6151,1000)*10**15,"sigma<=6.151e15")
import math
x=math.exp(-5e-4); print("float sigma",1+64*x*(1+4*x+x*x)/(1-x)**4+16*x/(1-x)**2)
N1=21**4-1; N5=37**4-1; assert (N1,N5)==(194480,1874160)
le(2*F(6151,1000)*10**15/N1**2,F(3253,10**3)*10**5,"s=1 bound"); le(2*F(6151,1000)*10**15/N5**2,F(3503,10**3)*10**3,"s=5 bound")
