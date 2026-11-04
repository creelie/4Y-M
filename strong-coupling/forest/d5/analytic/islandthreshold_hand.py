# Exact-rational enclosures for patches/A06_islandthreshold.tex.
from fractions import Fraction as F
def at2(x,m=12):  # interval for 2 artanh(x), 0<x<1
    s=sum(2*x**(2*k+1)/(2*k+1) for k in range(m+1)); return (s, s+2*x**(2*m+3)/(2*m+3)/(1-x*x))
def add(*I): return (sum(i[0] for i in I), sum(i[1] for i in I))
def mul(c,I): return (c*I[0],c*I[1]) if c>=0 else (c*I[1],c*I[0])
L2=at2(F(1,3)); L3=add(L2,at2(F(1,5))); L5=add(mul(2,L2),at2(F(1,9)))
L10=add(L2,L5); L11=add(L10,at2(F(1,21))); L23=add(mul(3,L2),L3,mul(-1,at2(F(1,47))))
def logq(x):  # interval log of a positive rational
    j=0
    while x>=2: x/=2; j+=1
    while x<1: x*=2; j-=1
    return add(mul(j,L2),at2((x-1)/(x+1)))
# pi bounds 3.14159265 < pi < 3.14159266
lpi=(logq(F(314159265,10**8))[0],logq(F(314159266,10**8))[1])
C0=add((F(1,2),F(1,2)),mul(F(2,3),add(mul(3,L2),L3,mul(3,lpi))))
H=1296
leMH=add((1,1),mul(4,L11),mul(4,add(L2,L3)))
LD1=add(mul(5,L2),L3,mul(2,L5))                 # log 2400
LD5=add(mul(4,L23),logq(F(279840,279841)))       # log(23^4-1)
def Lq(LD): return add(mul(H,add(L2,mul(2,LD))),leMH)
Ls={1:Lq(LD1),5:Lq(LD5)}
for s in (1,5): print("|log q*| s=%d in [%.4f, %.4f]"%(s,Ls[s][0],Ls[s][1]))
print("C0 in [%.7f,%.7f]"%(float(C0[0]),float(C0[1])))
def g_sign(beta,rho,L,which):  # g = beta rho^2/16 - 6 log beta - 6 C0 - L ; which: 'neg' need upper<0, 'pos' need lower>0
    lb=logq(beta)
    if which=='neg': return beta*rho**2/16-6*lb[0]-6*C0[0]-L[0] < 0
    return beta*rho**2/16-6*lb[1]-6*C0[1]-L[1] > 0
cases=[(F(1,10**8),1,0,F(3426,1000)*10**21,F(3428,1000)*10**21),
       (F(1,10**8),5,0,F(5400,1000)*10**21,F(5402,1000)*10**21),
       (F(1,10**9),1,0,F(3430,1000)*10**23,F(3432,1000)*10**23),
       (F(1,10**9),5,0,F(5404,1000)*10**23,F(5406,1000)*10**23),
       (F(1,10**8),5,30,F(1972,1000)*10**22,F(1974,1000)*10**22)]
for rho,s,gam,lo,hi in cases:
    L=add(Ls[s],mul(H*gam,L10))
    assert lo>96/rho**2
    ok=g_sign(lo,rho,L,'neg') and g_sign(hi,rho,L,'pos')
    B=(16*L[0]/rho**2,16*L[1]/rho**2); eps=(lo/B[1]-1,hi/B[0]-1)
    print(f"rho={float(rho)} s={s} Gamma=10^{gam}: beta* in [{float(lo):.4e},{float(hi):.4e}] {ok}; eps in [{float(eps[0]):.5f},{float(eps[1]):.5f}]")
    assert ok and 0<eps[0] and eps[1]<F(2,100)
print("--- digits and margins")
for s in (1,5): print(s, "%.6f %.6f"%(Ls[s][0],Ls[s][1]), float(Ls[s][1]-Ls[s][0]))
for rho,s,gam,lo,hi in cases:
    L=add(Ls[s],mul(H*gam,L10)); lb=logq(lo); hb=logq(hi)
    gl=lo*rho**2/16-6*lb[0]-6*C0[0]-L[0]; gh=hi*rho**2/16-6*hb[1]-6*C0[1]-L[1]
    print("margins",float(gl),float(gh)); assert gl<-F(1,2) and gh>F(1,2)
assert C0[0]>=F(490816,10**5) and C0[1]<=F(490817,10**5)
assert abs(Ls[1][0]-F(21090194,1000))<F(1,1000) and abs(Ls[5][0]-F(33424872,1000))<F(1,1000)
print("all ok")
