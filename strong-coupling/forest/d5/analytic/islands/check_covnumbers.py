# Hand chain of Corollaries fibrecovnumbers and labelstretched (patches P4, P5):
# every number recomputed from the rounded values stated in the text of
# Theorem islandwindow (patch P3), exact rationals with interval logarithms.
from ia import *
aG=F(25,10**5)
L19=add(L2,sc(2,L3),at2(F(1,37)))          # log 19 = log 18 + 2 artanh(1/37)
L45=add(sc(2,L3),L5)
log10=lambda x: log10I(x)
Lq={1:I(F(3248166,100),F(3248168,100)),5:I(F(3835406,100),F(3835408,100))}   # P3: within 1e-2
lMq={1:I(-F(14102447,1000),-F(14102445,1000)),5:I(-F(16652798,1000),-F(16652796,1000))}  # P3: within 1e-3
lb={1:I(16,F(162645,10**4)),5:I(16,F(163356,10**4))}                          # P3: b_*(beta*)
lnu={1:log10(add(L45,sc(4,L11))),5:log10(add(L45,sc(4,L19)))}
lbn={s:add(lb[s],lnu[s]) for s in (1,5)}
le(lbn[1],F(220833,10**4),"b_* nu_1 <= 10^22.0833"); le(lbn[5],F(231039,10**4),"b_* nu_5 <= 10^23.1039")
lbn={1:I(16,F(220833,10**4)),5:I(16,F(231039,10**4))}
lsA=log10(logq(F(9831,1000)*10**16))
l11=log10(sc(4,L11))
def mlg(s,R2):  # -log g_R at eps=Mq*: 1/4(|log q*|-4log11)-1/2log2-4log(2R+1)
    return sub(sc(F(1,4),sub(Lq[s],sc(4,L11))),add(sc(F(1,2),L2),sc(4,logq(R2))))
lg={1:mlg(1,23),5:mlg(5,39)}
for s,v in ((1,F(810513,100)),(5,F(957112,100))):
    le(sub(lg[s],I(v)),F(1,100),f"-log g_(2s+9) within 1e-2 of {float(v)}"); le(I(v-F(1,100)),lg[s],"(low)")
le(F(81051,10),lg[1],"-log g_11>=8105.1"); le(F(95710,10),lg[5],"-log g_19>=9571.0")
c={1:F(2626,100),5:F(2729,100)}; lam={1:F(736,10**4),5:F(503,10**4)}; f={1:F(70118,10),5:F(82859,10)}
a={1:F(2275,100),5:F(2377,100)}; ap={1:F(140626,10),5:F(166118,10)}
print("== fibrecovnumbers (i), (ii)")
for s in (1,5):
    R=2*s+9
    le(add(l11,lbn[s],log10(I(2*aG*(8*s+41)))),c[s],f"c_{s}")
    le(lam[s],sc(F(1,10**4*R),lg[s]),f"lambda_{s} <= -log g_R/(1e4 R)")
    le(2*aG,sc(F(1,10**4*R),lg[s]),"lambda/(1e4R)>=2 alpha_G")
    le(add(log10(L2),lbn[s],lsA,sc(F(1,2),lMq[s]),log10(I(aG*(2*s+8)))),-f[s],f"2 Phi_R <= 10^-f_{s}")
    le(add(log10(logq(F(9,2))),lbn[s],log10(I(4*aG*(2*s+9)))),a[s],f"a_{s}")
    le(add(log10(logq(F(17,4))),lbn[s],lsA,lMq[s],log10(I(aG*(2*s+8)))),-ap[s],f"a'_{s}")
le(1,lb[1][0]-1,"b_* nu_s >= 10^16 > 10/(9/4): b_*(beta) increasing, b_*(beta*)>=1920 r^2 99^2 beta*")
r=43*F(217,100)/10**16; lt(10**16,1920*r*r*99**2*F(1122,1000)*10**37,"b_*(beta*)>=10^16")
print("== fibrecovnumbers (iii)")
lgR={s:mlg(s,4*10**53+3) for s in (1,5)}
le(F(76239,10),lgR[1],"-log g_R>=7623.9"); le(F(90920,10),lgR[5],"-log g_R>=9092.0")
le(2000*10**50+1,F(4*10**53+3-1,2),"R<=2000e50+1 gives 2R+1<=4e53+3")
lt(F(1,2),F(76239,10)/8002,"lambda/8002>1/2")
lt(1000,10**4,"sqrt d >= 1e4 > 1/(4 alpha_G)=1000: main ratio decreasing")
LG={1:F(76239,10),5:F(90920,10)}
chain={1:-F(19652,10),5:-F(2761)}; main={1:-F(86634,10),5:-F(86624,10)}; third={1:-F(140627,10),5:-F(16612)}
for s in (1,5):
    ch=add(I(1),log10(I(-(LG[s]/8002-F(1,2))*10**4)))
    le(ch,chain[s],f"chain ratio s={s}")
    mn=add(log10(logq(F(9,4))),lbn[s],log10(I(-2*aG*(5*10**7-2*s-9)+5000)))
    le(mn,main[s],f"main ratio s={s}")
    th=add(log10(logq(F(17,4))),lbn[s],lsA,lMq[s],log10(I(aG*(4*s+17))))
    le(th,third[s],f"third coefficient s={s}")
print("== labelstretched")
lgS={s:mlg(s,4*10**52+3) for s in (1,5)}
le(F(76331,10),lgS[1],"-log g_R>=7633.1"); le(F(91012,10),lgS[5],"-log g_R>=9101.2")
le(200*10**50+1,F(4*10**52+3-1,2),"R<=200e50+1")
lt(F(1,20)+F(4,10**4),F(76331,10)/401,"lambda/401>1/20+4e-4")
LS={1:F(76331,10),5:F(91012,10)}
far={1:-F(2111422,100),5:-F(2493864,100)}
for s in (1,5):
    fa=add(log10(L2),lbn[s],lsA,sc(F(3,2),lMq[s]),log10(I(aG*(4*s+17))))
    le(fa,far[s],f"far coefficient x eps s={s} <= 10^-u_s")
    mn=add(l11,lbn[s],log10(I(2*aG*(2*s+14)-aG*10**8+500)),lMq[s]); le(mn,-24000,f"main ratio x eps s={s}")
    ex=sub(I(500),I((F(10**4,401)-3)*LS[s]))
    ch=add(log10(add(logq(4*401**4*10**16),ex)),lMq[s]); le(ch,-86000,f"chain ratio x eps s={s}")
    show(f"  values s={s}: far",fa); show("  main",mn); show("  chain",ch)
le(F(1,20),200*aG,"200 alpha_G=1/20")
print("== displayed values of the text (P4, P5)")
T={1:(F(262596,10**4),-F(701184,100),F(22742,1000),-F(1406273,100)),5:(F(272871,10**4),-F(828599,100),F(23766,1000),-F(1661206,100))}
for s in (1,5):
    R=2*s+9; cT,fT,aT,apT=T[s]
    le(add(l11,lbn[s],log10(I(2*aG*(8*s+41)))),cT,f"text c s={s}")
    le(add(log10(L2),lbn[s],lsA,sc(F(1,2),lMq[s]),log10(I(aG*(2*s+8)))),fT,f"text f s={s}")
    le(add(log10(logq(F(9,2))),lbn[s],log10(I(4*aG*(2*s+9)))),aT,f"text a s={s}")
    le(add(log10(logq(F(17,4))),lbn[s],lsA,lMq[s],log10(I(aG*(2*s+8)))),apT,f"text a' s={s}")
le(F(7368,10**5),sc(F(1,110000),lg[1]),"text 0.07368"); le(F(5037,10**5),sc(F(1,190000),lg[5]),"text 0.05037")
for s,fr,mn_,ch_ in ((1,-F(2111428,100),-24716,-86581),(5,-F(2493879,100),-27265,-103119)):
    le(add(log10(L2),lbn[s],lsA,sc(F(3,2),lMq[s]),log10(I(aG*(4*s+17)))),fr,f"text far s={s}")
    le(add(l11,lbn[s],log10(I(2*aG*(2*s+14)-aG*10**8+500)),lMq[s]),mn_,f"text main s={s}")
    le(add(log10(add(logq(4*401**4*10**16),sub(I(500),I((F(10**4,401)-3)*LS[s])))),lMq[s]),ch_,f"text chain s={s}")
print("all",OK[0],"checks passed")
