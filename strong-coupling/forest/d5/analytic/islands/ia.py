"""Exact-rational interval helpers for the hand bounds of patches P1-P7.
Only the operations used in the text: the artanh series for logarithms,
Taylor polynomials with geometric tails for exponentials, and rational
arithmetic."""
from fractions import Fraction as F
import math
def I(a,b=None): a=F(a); return (a,a) if b is None else (a,F(b))
def add(*xs): return (sum(x[0] for x in xs), sum(x[1] for x in xs))
def neg(x): return (-x[1],-x[0])
def sub(x,y): return add(x,neg(y))
def mul(x,y):
    p=[x[0]*y[0],x[0]*y[1],x[1]*y[0],x[1]*y[1]]; return (min(p),max(p))
def sc(c,x): return mul(I(c),x)
def inv(x): assert x[0]>0; return (1/x[1],1/x[0])
def div(x,y): return mul(x,inv(y))
def at2(x,m=14):
    """2 artanh(x) for 0<x<1/2: partial sum and the geometric tail bound"""
    x=F(x); s=sum(2*x**(2*k+1)/(2*k+1) for k in range(m+1))
    return (s, s+2*x**(2*m+3)/(2*m+3)/(1-x*x))
L2=at2(F(1,3)); L3=add(L2,at2(F(1,5))); L5=add(sc(2,L2),at2(F(1,9)))
L7=sub(sc(3,L2),at2(F(1,15))); L11=add(L2,L5,at2(F(1,21)))
L13=add(sc(2,L2),L3,at2(F(1,25))); L17=add(sc(4,L2),at2(F(1,33)))
L23=sub(add(sc(3,L2),L3),at2(F(1,47))); L37=add(sc(2,L2),sc(2,L3),at2(F(1,73)))
LN10=add(L2,L5)
def logq(x):
    """log of a positive rational: x=2^k y, 1<=y<2, log y=2artanh((y-1)/(y+1))"""
    x=F(x); k=0
    while x>=2: x/=2; k+=1
    while x<1: x*=2; k-=1
    return add(sc(k,L2),at2((x-1)/(x+1),40))
def logI(x): return (logq(x[0])[0],logq(x[1])[1])
def log1m(x):
    """log(1-x), 0<x<1/2: -x/(1-x) <= log(1-x) <= -x"""
    x=F(x); return (-x/(1-x),-x)
def expI(x,n=12):
    """e^x for an interval x in [-1,1]"""
    def up(t):
        if t<0: return 1/lo_(-t)
        return sum(t**k/math.factorial(k) for k in range(n+1))+t**(n+1)/math.factorial(n+1)/(1-t/(n+2))
    def lo_(t):
        if t<0: return 1/up(-t)
        return sum(t**k/math.factorial(k) for k in range(n+1))
    return (lo_(x[0]),up(x[1]))
def log10I(x): return div(x,LN10)   # x is a natural-log interval
PI=I(F(314159265,10**8),F(314159266,10**8))
SQ2=I(F(14142135,10**7),F(14142136,10**7))
assert SQ2[0]**2<=2<=SQ2[1]**2
C0=add(I(F(1,2)),sc(F(2,3),add(sc(3,L2),L3,sc(3,logI(PI)))))       # 1/2+2/3 log(24 pi^3)
C1=add(C0,sc(F(1,6),add(sc(9,L2),sc(2,L5))))                         # C0+1/6 log 12800
OK=[0]
def le(a,b,msg):
    """a<=b with a,b numbers or intervals (upper of a, lower of b)"""
    av=a[1] if isinstance(a,tuple) else F(a); bv=b[0] if isinstance(b,tuple) else F(b)
    assert av<=bv,(msg,float(av),float(bv)); OK[0]+=1
    print(f"ok  {msg}: {float(av):.10g} <= {float(bv):.10g}")
def lt(a,b,msg):
    av=a[1] if isinstance(a,tuple) else F(a); bv=b[0] if isinstance(b,tuple) else F(b)
    assert av<bv,(msg,float(av),float(bv)); OK[0]+=1
    print(f"ok  {msg}: {float(av):.10g} < {float(bv):.10g}")
def show(name,x,d=8): print(f"    {name} in [{float(x[0]):.{d}g}, {float(x[1]):.{d}g}]")
