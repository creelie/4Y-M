"""Exact-rational re-check of the hand bounds of Theorem su3 (SU(3) single link).
Every quantity is an upper or lower bound computed in Fractions; the text does not depend on this script."""
from fractions import Fraction as F
from math import factorial, comb
# Haar moments int t^a tbar^b over SU(3), from the hook-length counts in the text
f = {(1,1):1,(3,0):1,(2,2):2,(4,1):3,(6,0):5,(3,3):6,(5,2):11,(7,1):21,(4,4):23,(0,0):1}
def mom(a,b):
    if (a-b)%3: return 0
    return f.get((a,b), f.get((b,a)))
m=[]
for n in range(9):
    s=sum(comb(n,a)*mom(a,n-a) for a in range(n+1))
    m.append(F(s,2**n))
print('m_n =',[str(x) for x in m])
assert m[2:]==[F(1,2),F(1,4),F(3,4),F(15,16),F(65,32),F(231,64),F(973,128)]
def zeta(r):
    head=sum(m[n]*r**n/factorial(n) for n in range(9))
    x=3*r; t9=x**9/factorial(9)
    assert x<10
    tail=m[8]/3**8 * t9/(1-x/10)   # sum_{n>=9} x^n/n! <= t9/(1-x/10)
    return head,tail
def explow(x):  # e^{-x} >= 1-x+x^2/2-x^3/6 for x>=0
    return 1-x+x**2/2-x**3/6
def sqrtup(y, d=6):  # decimal upper bound for sqrt(y)
    s=F(int(float(y)**0.5*10**d)+1,10**d); assert s*s>=y; return s
for bw in [F(1,2),F(515,1000),F(9,29)]:
    r=2*bw; head,tail=zeta(r); Z=head+tail
    e=1-F(5,9)*explow(r/5)/Z
    g=F(283,128)-r/2
    s=sqrtup(e)
    a=2*r*s/g
    print('beta_W',bw,float(bw),'head',float(head),'tail',float(tail),'Z<=',float(Z),'e(r)<=',float(e),'sqrt<=',float(s),'gamma',g,float(g),'alpha<=',float(a))

print('--- checks of the decimals written in the proof of Theorem su3 ---')
def up(x,d): # decimal upper rounding
    from math import ceil
    return F(ceil(x*10**d),10**d)
for r,terms_txt,sumtxt,tailtxt,Ztxt,exptxt,ratiotxt,epstxt,sq,gtxt in [
   (F(1),['1','0','0.25','0.0416667','0.03125','0.0078125','0.0028212','0.0007162','0.0001886'],'1.3344552','0.00009','1.334546','0.818666','0.340800','0.659200','0.81192',F(219,128)),
   (F(103,100),['1','0','0.265225','0.0455303','0.0351722','0.0090569','0.0033687','0.0008808','0.0002389'],'1.3594728','0.00012','1.359593','0.813760','0.332517','0.667483','0.817',F(5427,3200))]:
    terms=[m[n]*r**n/factorial(n) for n in range(9)]
    for t,tt in zip(terms,terms_txt): assert t<=F(tt), (t,tt)
    assert sum(F(tt) for tt in terms_txt)<=F(sumtxt)
    x=3*r; tail=m[8]/3**8*x**9/factorial(9)/(1-x/10); assert tail<=F(tailtxt),float(tail)
    assert F(sumtxt)+F(tailtxt)<=F(Ztxt)
    assert explow(r/5)>=F(exptxt)
    assert F(5,9)*F(exptxt)/F(Ztxt)>=F(ratiotxt)
    assert 1-F(ratiotxt)<=F(epstxt) and F(epstxt)<F(sq)**2
    assert F(283,128)-r/2==gtxt
    print(r,'alpha<=',float(2*r*F(sq)/gtxt))
assert 2*F(1)*F('0.81192')/F(219,128)<=F('0.9491')
assert 2*F(103,100)*F('0.817')/F(5427,3200)<=F('0.9924')
y=1-F('0.9491'); assert y+y*y/2>=F('0.052')
y=1-F('0.9924'); assert y+y*y/2>=F('0.0076')
print('all decimals of the proof confirmed')
