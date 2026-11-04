# Hand bounds of Theorems islandscales and islandswide(ii),(iii) (patch P1, P2).
from ia import *
print("== logarithms of |log q*_j| (s=5) and of the j=0 thresholds")
LD5=add(sc(4,L23),log1m(F(1,23**4)))                 # log(23^4-1)
LD1=add(sc(5,L2),L3,sc(2,L5))                          # log 2400
base=lambda H,LD: sc(H,add(L2,sc(2,LD)))
def Lq(j,LD=LD5):
    if j==0: return add(base(1296,LD),I(1),sc(4,L11),sc(4,add(L2,L3)))
    H=4096 if j<=2 else 6561
    LH=sc(12,L2) if j<=2 else sc(8,L3)
    return add(base(H,LD),I(1),sc(4,logq(17*2**j-4)),LH)
A3=add(base(6561,LD5),I(1),sc(8,L3))                   # |log q*_j| - 4log(17 2^j-4) for j>=3
show("6561(log2+2logDelta)+1+8log3",A3,12)
for j in (1,2,3): show(f"|log q*_{j}|",Lq(j),12)
show("|log q*| s=5,j=0",Lq(0),12); show("|log q*| s=1,j=0",Lq(0,LD1),12)
c=F(101,100)*F(1023,1000); c2=c*c; print("    c^2 =",c2,float(c2))
def g(beta,rho0,L):  # g = beta rho0^2/16 - 6(log beta + C0) - L, interval
    return sub(I(beta*rho0**2/16),add(sc(6,add(logq(beta),C0)),L))
def rho_j(rho,j): return rho if j==0 else rho/(c*4**j)
def bracket(rho,j,lo,hi,L=None,tag=""):
    L=L or Lq(j); r0=rho_j(rho,j)
    assert lo>96/r0**2
    gl=g(lo,r0,L); gh=g(hi,r0,L)
    lt(gl,0,f"g<0 at {float(lo):.5e} {tag}"); lt(0,gh,f"g>0 at {float(hi):.5e} {tag}")
    print(f"      margins {float(gl[1]):.3f} {float(gh[0]):.3f}")
    assert gl[1]<-F(1,2) and gh[0]>F(1,2)
print("== Theorem islandscales, rho=1e-8, s=5")
R8=F(1,10**8)
bracket(R8,1,F(2895,1000)*10**23,F(2897,1000)*10**23,tag="beta*_1")
bracket(R8,2,F(4633,1000)*10**24,F(4635,1000)*10**24,tag="beta*_2")
bracket(R8,3,F(1185,1000)*10**26,F(1187,1000)*10**26,tag="beta*_3")
# 16^{-j} beta*_j <= b  <=>  g_j(16^j b) >= 0 ; g_j(16^j b) = b rho^2/(16c^2) - 6(log b + C0) - (L_j + 24 j log 2)
def G(rho,j,b): return sub(I(b*rho**2/(16*c2)),add(sc(6,add(logq(b),C0)),Lq(j),sc(24*j,L2)))
for j in range(1,13): lt(0,G(R8,j,F(29,10)*10**22),f"16^-j beta*_j<=2.90e22, j={j}")
print("   the j-dependence: L_j+24 j log2 increases in j (H_j, M_j increase)")
# bounds for the k_* comparison
for j in range(1,18): lt(0,G(R8,j,F(3)*10**22),f"16^-j beta*_j<3e22 j={j}"); lt(G(R8,j,F(54,10)*10**21),0,f"16^-j beta*_j>5.4e21 j={j}")
print("   j=0: beta*(1e-8,5) in [5.400,5.402]e21 (Theorem islandthreshold), so 5.4e21<=c_0<3e22")
# D(beta) = 1.0201*96e8 (log beta + C1): 16^{k_*} <= beta/D
Dlo=mul(I(F(10201,10**4)*96*10**8),add(sc(22,LN10),C1)); Dhi=mul(I(F(10201,10**4)*96*10**8),add(sc(42,LN10),C1))
show("D(1e22)",Dlo); show("D(1e42)",Dhi)
le(16**8*Dhi[1],F(54,10)*10**21,"16^8 D(1e42) <= 5.4e21  (k*-jmax>=8)")
lt(F(3)*10**22,16**9*Dlo[0],"3e22 < 16^9 D(1e22)  (k*-jmax<=9)")
lt(10**42,16**17*F(54,10)*10**21,"1e42 < 16^17*5.4e21 (j<=17 suffices)")
print("== Theorem islandswide(ii), rho=1e-4")
R4=F(1,10**4)
bracket(R4,0,F(3408,1000)*10**13,F(3410,1000)*10**13,L=Lq(0,LD1),tag="beta*(1e-4,1)")
bracket(R4,0,F(5382,1000)*10**13,F(5384,1000)*10**13,tag="beta*(1e-4,5)")
bracket(R4,1,F(2892,1000)*10**15,F(2894,1000)*10**15,tag="beta*_1")
bracket(R4,2,F(4628,1000)*10**16,F(4630,1000)*10**16,tag="beta*_2")
bracket(R4,3,F(1185,1000)*10**18,F(1186,1000)*10**18,tag="beta*_3")
for j in range(1,21): lt(0,G(R4,j,F(29,10)*10**14),f"16^-j beta*_j<=2.90e14, j={j}")
print("== Theorem islandswide(iii)")
for j in range(6,24): lt(0,G(R4,j,F(292,100)*10**14),f"c_j<2.92e14 j={j}"); lt(G(R4,j,F(288,100)*10**14),0,f"c_j>2.88e14 j={j}")
lt(16**6*F(292,100)*10**14,10**22,"16^6*2.92e14 < 1e22, so j_max>=6")
lt(10**42,16**23*F(288,100)*10**14,"1e42 < 16^23*2.88e14, so j_max<=22")
le(16**2*Dhi[1],F(288,100)*10**14,"16^2 D(1e42) <= 2.88e14 (k*-jmax>=2)")
lt(F(292,100)*10**14,16**3*Dlo[0],"2.92e14 < 16^3 D(1e22) (k*-jmax<=3)")
print("all",OK[0],"checks passed")
print("== displayed intermediate values of the text")
c216=16*c2; le(c216,F(170810278,10**7),"16c^2<=17.0810278")
A=A3; le(sub(A,I(F(169133301,1000))),F(1,1000),"A within 1e-3 of 169133.301"); le(F(169133301,1000)-F(1,1000),A,"A lower")
for j,v in ((1,F(105605899,1000)),(2,F(105608929,1000)),(3,F(169152832,1000))):
    le(sub(Lq(j),I(v)),F(1,1000),f"L_{j} - {float(v)} <= 1e-3"); le(I(v-F(1,1000)),Lq(j),f"L_{j} lower")
T3=lambda J: add(A,sc(4,L17),sc(28*J,L2))
le(T3(1)[1]-28*L2[1],F(16914464,100),"A+4log17<=169144.64")
def T1(b,rho): return b*rho**2/(16*c2)
def T2(b): return sc(6,add(logq(b),C0))
for rho,b,J,t1,t2,t3 in ((R8,F(29,10)*10**22,12,F(1697790,10),F(3398,10),F(1693776,10)),
                          (R8,F(3)*10**22,17,F(1756334,10),F(3400,10),F(1694746,10)),
                          (R4,F(29,10)*10**14,20,F(1697790,10),F(2293,10),F(1695328,10)),
                          (R4,F(292,100)*10**14,23,F(1709499,10),F(2293,10),F(1695911,10))):
    le(t1,T1(b,rho),f"T1>={float(t1)}"); le(T2(b),t2,f"T2<={float(t2)}"); le(T3(J),t3,f"T3({J})<={float(t3)}")
    lt(t2+t3,t1,"sum")
le(T1(F(54,10)*10**21,R8),F(316141,10),"T1(5.4e21)<=31614.1"); lt(F(316141,10),Lq(1),"<L_1")
le(T1(F(288,100)*10**14,R4),F(1686082,10),"T1(2.88e14)<=168608.2"); lt(F(1686082,10),Lq(3),"<L_3")
le(F(648436,10**5),C1,"C1>=6.48436"); le(C1,F(648437,10**5),"C1<=6.48437")
le(F(5595,10**4)*10**12,Dlo,"D(1e22)>=5.595e11"); le(Dhi,F(10106,10**4)*10**12,"D(1e42)<=1.0106e12")
le(16**8*F(10106,10**4)*10**12,F(435,100)*10**21,"16^8 D<=4.35e21"); le(F(384,100)*10**22,16**9*F(5595,10**4)*10**12,"16^9 D>=3.84e22")
le(16**2*F(10106,10**4)*10**12,F(259,100)*10**14,"16^2 D<=2.59e14"); le(F(229,100)*10**15,16**3*F(5595,10**4)*10**12,"16^3 D>=2.29e15")
print("== displayed values of the text of islandswide (P2)")
v=F(21090194,1000); le(sub(Lq(0,LD1),I(v)),F(1,1000),"|log q*|(s=1,j=0) within 1e-3 of 21090.194"); le(I(v-F(1,1000)),Lq(0,LD1),"(low)")
for j,lo,hi,L in ((0,F(3408,1000)*10**13,F(3410,1000)*10**13,Lq(0,LD1)),(0,F(5382,1000)*10**13,F(5384,1000)*10**13,Lq(0)),
                  (1,F(2892,1000)*10**15,F(2894,1000)*10**15,Lq(1)),(2,F(4628,1000)*10**16,F(4630,1000)*10**16,Lq(2)),(3,F(1185,1000)*10**18,F(1186,1000)*10**18,Lq(3))):
    r0=rho_j(R4,j); lt(g(lo,r0,L),-5,f"margin>5 low j={j}"); lt(5,g(hi,r0,L),f"margin>5 high j={j}")
le(96*c2/R4**2,F(288,100)*10**14,"b>96c^2/rho^2 at rho=1e-4")
le(96*c2/R8**2,F(54,10)*10**21,"b>96c^2/rho^2 at rho=1e-8")
print("all",OK[0],"checks passed")
