# Hand bounds of Theorem islandwindow, Corollaries fibrecovnumbers and labelstretched (patches P3-P5).
from ia import *
rho=F(217,100)*F(1,10**16); r=43*rho; aG=F(25,10**5)
H=1296; M=11**4
le(2758*rho,F(6,10)*F(1,10**12),"2758 rho <= 6.0e-13")
LN={1:add(sc(4,L2),L5,L11,L13,L17),          # log 194480 = log(2^4 5 11 13 17)
    5:add(sc(4,L37),log1m(F(1,37**4)))}      # log(37^4-1)
assert 194480==2**4*5*11*13*17 and 21**4-1==194480
LMH=add(I(1),sc(4,L11),sc(4,add(L2,L3)))     # 1+log M+log H
Lq={s:add(sc(H,add(L2,sc(2,LN[s]))),LMH) for s in (1,5)}
for s in (1,5): show(f"|log q*| s={s}",Lq[s],12)
le(Lq[1],F(3248168,100),"|log q*|(s=1)<=32481.68"); le(F(3248166,100),Lq[1],"|log q*|(s=1)>=32481.66")
le(Lq[5],F(3835408,100),"|log q*|(s=5)<=38354.08"); le(F(3835406,100),Lq[5],"|log q*|(s=5)>=38354.06")
def g(beta,L): return sub(I(beta*rho**2/16),add(sc(6,add(logq(beta),C0)),L))
B={1:(F(11220,10**4)*10**37,F(11221,10**4)*10**37),5:(F(13216,10**4)*10**37,F(13217,10**4)*10**37)}
for s in (1,5):
    lo,hi=B[s]; assert lo>96/rho**2
    a=g(lo,Lq[s]); b=g(hi,Lq[s]); lt(a,-F(1,2),f"g(lo)<-1/2 s={s}"); lt(F(1,2),b,f"g(hi)>1/2 s={s}")
# decimal logarithms
lMq={s:log10I(sub(sc(4,L11),Lq[s])) for s in (1,5)}
lH=log10I(sc(4,add(L2,L3))); lA={s:add(lMq[s],log10I(I(1)),lH) for s in (1,5)}
for s in (1,5): show(f"log10 Mq* s={s}",lMq[s],10); show(f"log10 A* s={s}",lA[s],10)
le(lMq[1],-F(141024,10),"Mq*<=10^-14102.4 (s=1)"); le(lA[1],-F(140988,10),"A<=10^-14098.8")
le(lMq[5],-F(166527,10),"Mq*<=10^-16652.7 (s=5)")
# G_0 = H(Mq* + A* 2^-H/(1-2^-H)) <= H Mq*(1+2eH 2^-H) ; 2eH2^-H < 10^-300
le(2*3*H,F(10)**(-300)*2**H,"2eH 2^-H <= 1e-300")
lG={s:add(lMq[s],lH,I(0,F(1,10**290))) for s in (1,5)}
le(lG[1],-F(1409930,100),"G0<=10^-14099.3 (s=1)"); le(lG[5],-F(1664960,100),"G0<=10^-16649.6 (s=5)")
# sigma(5e-4) <= 6.151e15 (Proposition islandmass); log10 6.151e15
lsig5=log10I(logq(F(6151,1000)*10**15))
le(add(lG[1],lsig5),-F(140835,10),"G0 sigma(2aG) <= 10^-14083.5")
# sigma at alpha_G and alpha_G/2: with u=1-e^{-k}>=k-k^2/2, x(1+4x+x^2)=(1-u)(6-6u+u^2)<=6-12u+7u^2
def sigma_up(k):
    u=k-k*k/2; return 1+64*(6-12*u+7*u*u)/u**4+16/u**2
le(sigma_up(aG),F(9831,1000)*10**16,"sigma(alpha_G)<=9.831e16")
lsigA=log10I(logq(F(9831,1000)*10**16))
lsigA2=log10I(logq(sigma_up(aG/2))); show("log10 sigma(alpha_G/2)",lsigA2)
# monotonicity of 1-12u... : the bound 64(6-12u+7u^2)/u^4 decreases in u on (0,1/2): fine
# Lambda_2 bar at beta*: 1920/b (99.0001 b+4.65e6)^2 r^2 = 1920 b r^2 (99.0001+4.65e6/b)^2, increasing in b
def lbstar(s):
    b=B[s][1]; return log10I(logq(1920*b*r*r*(F(990001,10**4)+F(465,100)*10**6/B[s][0])**2))
nu={1:45*11**4,5:45*19**4}
lbn={s:add(lbstar(s),log10I(logq(nu[s]))) for s in (1,5)}
for s in (1,5): show(f"log10 b_* nu_s s={s}",lbn[s],10)
le(lbn[1],F(22084,1000),"b_* nu_1 <= 10^22.084"); le(lbn[5],F(23104,1000),"b_* nu_5 <= 10^23.104")
# c_2 (theta=0, kappa=2aG): b_* 45(2s+12)^2(2s+7)^2 e^{kappa(2s+14)}
def lc2(s,kap): return add(lbstar(s),log10I(add(logq(45*(2*s+12)**2*(2*s+7)**2),I(kap*(2*s+14)))))
le(lc2(1,2*aG),F(2213,100),"c_2<=10^22.13 (s=1,theta=0)")
# epsilon
def lsum(*xs): # log10 of a sum of positive terms given log10 intervals: upper bound
    m=max(x[1] for x in xs); return m+F(1,10**6) if sum(1 for x in xs if x[1]>m-20)==1 else None
l11=log10I(sc(4,L11))
# R=15, theta=0: eps <= |W| Mq* (11^4 31^4 Mq* + c2 G0 tau/sigma), tau<=sigma(2aG)
t1=add(lMq[1],l11,log10I(logq(31**4)),lMq[1]); show("R=15 first term",t1,10)
le(t1,-28194,"first term <= 10^-28194")
t2=add(lMq[1],lc2(1,2*aG),lG[1],lsig5); show("R=15 second coefficient",t2,10)
le(t2,-F(281638,10),"second <= 10^-28163.8/sigma")
# R=2e8: tau <= e^{-kappa R/2} sigma(kappa/2), kappa R/2 = 5e4
R=2*10**8
t1=add(lMq[1],l11,log10I(logq((2*R+1)**4)),lMq[1]); show("R=2e8 first term",t1,10)
ltau=add(log10I(I(-F(5,10**4)*R/2)),lsigA)
t2=add(lMq[1],lc2(1,2*aG),lG[1],ltau,I(1000)); show("R=2e8 second term at sigma=1e-1000",t2,10)
le(add(t1,I(0,F(1,10**6))),-28166,"eps<=1e-28166|W| (theta=0,R=2e8)"); lt(t2,t1[0]-50,"second term negligible")
# theta=1/2, kappa=aG, R=4e8: tau <= e^{-aG R/2} sigma(aG/2)
R=4*10**8
t1=add(lMq[1],l11,log10I(logq((2*R+1)**4)),lMq[1]); show("theta=1/2 R=4e8 first",t1,10)
ltau=add(log10I(I(-aG*R/2)),lsigA2)
t2=add(lMq[1],lc2(1,aG),lG[1],ltau,I(1000)); lt(t2,t1[0]-50,"second negligible")
le(add(t1,I(0,F(1,10**6))),-28165,"eps<=1e-28165|W| (theta=1/2)")
# s=5, theta=0, R=2e8
R=2*10**8
t1=add(lMq[5],l11,log10I(logq((2*R+1)**4)),lMq[5]); show("s=5 R=2e8 first",t1,10)
ltau=add(log10I(I(-F(5,10**4)*R/2)),lsigA)
t2=add(lMq[5],lc2(5,2*aG),lG[5],ltau,I(1000)); lt(t2,t1[0]-50,"second negligible")
le(add(t1,I(0,F(1,10**6))),-33267,"eps<=1e-33267|W| (s=5)")
print("== Corollary fibrecovnumbers")
lg={}
for s in (1,5):
    R=2*s+9
    # -log g_R = -1/2 log2 - 4 log(2R+1) - 1/4 log eps >= 1/4(L - log M) - 1/2 log 2 - 4 log(2R+1)
    lg[s]=sub(sc(F(1,4),sub(Lq[s],sc(4,L11))),add(sc(F(1,2),L2),sc(4,logq(2*R+1))))
    show(f"-log g_(2s+9) s={s}",lg[s],10)
le(F(81051,10),lg[1],"-log g_11 >= 8105.1"); le(F(95710,10),lg[5],"-log g_19 >= 9571.0")
lnten=LN10
cs={1:F(2626,100),5:F(2729,100)}; lam={1:F(736,10**4),5:F(503,10**4)}; fs={1:F(70118,10),5:F(82859,10)}
as_={1:F(2275,100),5:F(2377,100)}; aps={1:F(140626,10),5:F(166118,10)}
for s in (1,5):
    R=2*s+9
    c=add(l11,lbn[s],log10I(I(2*aG*(8*s+41)))); show(f"c_{s}",c); le(c,cs[s],f"c_{s}")
    le(lam[s]*10**4*R,lg[s],f"lambda_{s}")
    f_=add(log10I(L2),lbn[s],lsigA,sc(F(1,2),lMq[s]),log10I(I(aG*(2*s+8)))); show(f"-f_{s}",f_); le(f_,-fs[s],f"f_{s}")
    a=add(log10I(logq(F(9,2))),lbn[s],log10I(I(4*aG*(2*s+9)))); show(f"a_{s}",a); le(a,as_[s],f"a_{s}")
    ap=add(log10I(logq(F(17,4))),lbn[s],lsigA,lMq[s],log10I(I(aG*(2*s+8)))); show(f"-a'_{s}",ap); le(ap,-aps[s],f"a'_{s}")
# (iii): R<=2000*1e50+1, 2R+1<=4e53+3
lgR={}
for s in (1,5):
    lgR[s]=sub(sc(F(1,4),sub(Lq[s],sc(4,L11))),add(sc(F(1,2),L2),sc(4,logq(4*10**53+3))))
    show(f"-log g_R (iii) s={s}",lgR[s])
le(F(76239,10),lgR[1],"-log g_R>=7623.9"); le(F(90920,10),lgR[5],"-log g_R>=9092.0")
# at d=1e8: chain ratio 10 g_R^{sqrt d/8002} e^{sqrt d/2} <= 10 e^{-(lambda/8002-1/2) 1e4}
for s in (1,5):
    ch=add(I(1),log10I(sc(-1,sub(sc(F(10**4,8002),lgR[s]),I(5000))))); show(f"chain ratio s={s}",ch)
    le(ch,-F(19652,10),f"chain <=10^-1965.2 s={s}")
    # main: 9/4 b nu e^{-2aG(d/2-2s-9)} e^{sqrt d/2} at d=1e8
    mn=add(log10I(logq(F(9,4))),lbn[s],log10I(I(-2*aG*(5*10**7-2*s-9)+5000))); show(f"main ratio s={s}",mn)
    le(mn,-8662,f"main <= 10^-8662 s={s}")
    th=add(log10I(logq(F(17,4))),lbn[s],lsigA,lMq[s],log10I(I(aG*(4*s+17)))); show(f"third coefficient s={s}",th)
    le(th,-14062,f"third <= 10^-14062 s={s}")
le(F(1,2),F(76239,10)/8002,"lambda/8002>1/2 (chain ratio decreasing)")
print("== Corollary labelstretched")
lgS={}
for s in (1,5):
    lgS[s]=sub(sc(F(1,4),sub(Lq[s],sc(4,L11))),add(sc(F(1,2),L2),sc(4,logq(4*10**52+3))))
    show(f"-log g_R (stretched) s={s}",lgS[s])
le(F(76331,10),lgS[1],"-log g_R>=7633.1"); le(F(91012,10),lgS[5],"-log g_R>=9101.2")
le(F(1,20)+F(4,10**4),F(76331,10)/401,"lambda/401>1/20+4e-4")
us={1:F(2111422,100),5:F(2493864,100)}
for s in (1,5):
    far=add(log10I(L2),lbn[s],lsigA,sc(F(3,2),lMq[s]),log10I(I(aG*(4*s+17)))); show(f"far coefficient x eps s={s}",far)
    le(far,-us[s],f"far <= 10^-u_{s}")
    mn=add(l11,lbn[s],log10I(I(2*aG*(2*s+14)-aG*10**8+500)),lMq[s]); show(f"main ratio x eps s={s}",mn); le(mn,-24000,"main")
    # chain: 4 401^4 d^2 e^{3 lam} e^{-lam sqrt d/401} e^{sqrt d/20} eps, at d=1e8
    ch=add(log10I(add(logq(4*401**4*10**16),sc(3,sc(-1,sc(-1,lgS[s]))) )),I(0))
    # careful: e^{3 lambda} with lambda=-log g_R is an upper factor; e^{-lambda sqrt d/401} decreases with lambda;
    # combined exponent 3 lam - lam 1e4/401 = -lam (1e4/401-3) <= -lgS_low (1e4/401-3)
    ex=sub(I(500),sc(F(10**4,401)-3,lgS[s]))
    ch=add(log10I(add(logq(4*401**4*10**16),ex)),lMq[s]); show(f"chain ratio x eps s={s}",ch); le(ch,-86000,"chain")
print("== displayed intermediate values of the text of islandwindow (P3)")
for s_,v in ((1,-F(14102446,1000)),(5,-F(16652797,1000))):
    le(sub(lMq[s_],I(v)),F(1,1000),f"log10 Mq* within 1e-3 of {float(v)} (up)"); le(I(v-F(1,1000)),lMq[s_],"(low)")
le(log10I(add(I(1),sc(4,add(L2,L3)))),F(3547,1000),"log10(eH)<=3.547")
le(lH,F(31127,10**4),"log10 H<=3.1127"); le(add(lMq[1],I(F(31127,10**4))),-F(140993,10),"G0 via 3.1127 (s=1)"); le(add(lMq[5],I(F(31127,10**4))),-F(166496,10),"G0 via 3.1127 (s=5)")
le(lbstar(1),F(162645,10**4),"b(beta*)<=10^16.2645 s=1"); le(lbstar(5),F(163356,10**4),"b(beta*)<=10^16.3356 s=5")
lt(F(465,100)*10**6/F(990001,10**4),B[1][0],"beta*>4.65e6/99.0001: b(beta) increasing")
le(lc2(1,2*aG),F(22122,1000),"c2<=10^22.122 s=1"); le(lc2(5,2*aG),F(23140,1000),"c2<=10^23.140 s=5")
le(lc2(1,aG),lc2(1,2*aG),"c2 smaller at kappa=aG")
le(sigma_up(aG/2),F(1573,1000)*10**18,"sigma(aG/2)<=1.573e18")
lsA=log10I(logq(F(9831,1000)*10**16)); lsA2=log10I(logq(F(1573,1000)*10**18))
l5e4=-21714; le(log10I(I(-50000)),l5e4,"e^-5e4<=10^-21714")
R=15; le(add(lMq[1],l11,log10I(logq((2*R+1)**4)),lMq[1]),-F(2819476,100),"first R=15 <= -28194.76")
R=2*10**8; le(add(lMq[1],l11,log10I(logq((2*R+1)**4)),lMq[1]),-F(2816631,100),"first R=2e8 <= -28166.31")
le(add(lMq[5],l11,log10I(logq((2*R+1)**4)),lMq[5]),-F(3326701,100),"first s=5 R=2e8 <= -33267.01")
R=4*10**8; le(add(lMq[1],l11,log10I(logq((2*R+1)**4)),lMq[1]),-F(2816511,100),"first R=4e8 <= -28165.11")
le(add(lMq[1],I(F(22122,1000)),lG[1],lsig5),-F(2816386,100),"R=15 second <= -28163.86")
for tag,lm,c2b,ls in (("s=1 R=2e8",lMq[1],F(22122,1000),lsA),("s=1 theta=1/2",lMq[1],F(22122,1000),lsA2),("s=5 R=2e8",lMq[5],F(23140,1000),lsA)):
    lg0=lG[1] if lm is lMq[1] else lG[5]
    le(add(lm,I(c2b),lg0,I(l5e4),ls,I(1000)),-48000,f"second <= 10^-48000 ({tag})")
# 10^v+10^-48000 <= 10^b when v-b <= -0.01: 10^-0.01 <= 0.9772 and 10^(-48000-b) <= 10^-14000
for v,b in ((-F(2816631,100),-28166),(-F(2816511,100),-28165),(-F(3326701,100),-33267)):
    le(v-b,-F(1,100),f"margin >= 0.01 below 10^{b}")
le(F(9772,10**4)**100,F(1,10),"10^-0.01<=0.9772")
print("all",OK[0],"checks passed")
