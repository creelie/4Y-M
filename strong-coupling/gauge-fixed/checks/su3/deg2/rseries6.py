"""Tilted SU(3) moments of degree <= DEG at diagonal b with rigorous error bounds (degree-6 version of
rseries.py, same expansion and bookkeeping; series products in C, the sum over (p, q) regrouped by the
powers of det D and det D').
|computed - exact| <= gamma(n) * ab_exact coefficientwise; truncation by Cauchy estimates."""
import numpy as np, math
from ring6 import FRing, EX
from series import hookprod, dim
from rseries import gam, U_, trunc_tail
from cmul import cmul
NMUL=80          # a degree-6 output coefficient collects at most 2^6 = 64 products
class TS:
    __slots__=('v','a','n')
    def __init__(s,v,a,n): s.v=v; s.a=a; s.n=n
def mkconst(R,c):
    v=np.zeros(R.N,complex); v[0]=c; a=np.zeros(R.N); a[0]=abs(c); return TS(v,a,0)
def mkvar(R,k):
    v=np.zeros(R.N,complex); v[R.idx[tuple(1 if i==k else 0 for i in range(18))]]=1; return TS(v,np.abs(v),0)
def add(x,y): return TS(x.v+y.v,x.a+y.a,max(x.n,y.n)+2)
def sub(x,y): return TS(x.v-y.v,x.a+y.a,max(x.n,y.n)+2)
def mul(R,x,y):
    v,a=cmul(R,x.v,y.v,x.a,y.a); return TS(v,a,x.n+y.n+NMUL)
def det3(R,M):
    m=lambda a,b: mul(R,a,b)
    t1=m(M[0][0],sub(m(M[1][1],M[2][2]),m(M[1][2],M[2][1])))
    t2=m(M[0][1],sub(m(M[1][0],M[2][2]),m(M[1][2],M[2][0])))
    t3=m(M[0][2],sub(m(M[1][0],M[2][1]),m(M[1][1],M[2][0])))
    return add(sub(t1,t2),t3)
def clam(l): return 1.0/(2.0**sum(l)*float(hookprod(l)))
def lincomb(coefs,xs):
    """sum_i coefs[i] xs[i] for real float coefficients (relative error gamma(8) each), any summation order"""
    V=np.array([x.v for x in xs]); A=np.array([x.a for x in xs]); c=np.asarray(coefs,float)
    return TS(c@V,np.abs(c)@A,max(x.n for x in xs)+8+2*len(xs))
def Zseries_rig(R,b,P=26,K=9):
    assert 3*(K+1)>P
    D=[[add(mkconst(R,b[j] if i==j else 0),mkvar(R,3*j+i)) for i in range(3)] for j in range(3)]
    Dp=[[add(mkconst(R,np.conj(b[j]) if i==j else 0),mkvar(R,9+3*j+i)) for i in range(3)] for j in range(3)]
    def ssum(L):
        out=L[0]
        for x in L[1:]: out=add(out,x)
        return out
    m=lambda a,b: mul(R,a,b)
    X=[[ssum([m(D[a][k],Dp[c][k]) for k in range(3)]) for c in range(3)] for a in range(3)]
    e1=ssum([X[0][0],X[1][1],X[2][2]])
    e2=ssum([sub(m(X[0][0],X[1][1]),m(X[0][1],X[1][0])),sub(m(X[0][0],X[2][2]),m(X[0][2],X[2][0])),sub(m(X[1][1],X[2][2]),m(X[1][2],X[2][1]))])
    e3=det3(R,X)
    h=[mkconst(R,1.0),e1,sub(m(e1,e1),e2)]
    for k in range(3,P+2): h.append(add(sub(m(e1,h[k-1]),m(e2,h[k-2])),m(e3,h[k-3])))
    prods={}
    def hh(i,j):
        key=(min(i,j),max(i,j))
        if key not in prods: prods[key]=h[key[1]] if key[0]==0 else m(h[key[0]],h[key[1]])
        return prods[key]
    pairs=[(p,q) for p in range(P+1) for q in range(p+1)]
    S=[hh(p,0) if q==0 else sub(hh(p,q),hh(p+1,q-1)) for p,q in pairs]
    dD=det3(R,D); dDp=det3(R,Dp)
    pw=[mkconst(R,1.0)]; pwp=[mkconst(R,1.0)]
    for k in range(1,K+1): pw.append(m(pw[-1],dD)); pwp.append(m(pwp[-1],dDp))
    Z=None
    for k in range(K+1):
        inner=None
        for l in range(K+1):
            co=[clam((p+k,q+k,k))*clam((p+l,q+l,l))/dim(p,q) for p,q in pairs]
            Skl=lincomb(co,S)
            t=Skl if l==0 else m(pwp[l],Skl)
            inner=t if inner is None else add(inner,t)
        t=inner if k==0 else m(pw[k],inner)
        Z=t if Z is None else add(Z,t)
    return Z
def moments_rig(R,b,P=26,K=9):
    """(w, ew): w[k] approximates E_H[u^alpha ubar^beta e^{Re Tr bU}] (u_ij <-> var 3i+j, ubar_ij <-> 9+3i+j),
    ew[k] bounds the absolute error. Monomials of nonzero torus weight have exact moment 0 (conjugation by the
    diagonal torus preserves Haar measure and Re Tr(Ub)); they are set to 0 with error 0."""
    b=np.asarray(b,complex); DEG=R.deg
    Z=Zseries_rig(R,b,P,K)
    g=gam(Z.n); errZ=g*Z.a/(1-g)
    perm=[0]*18
    for i in range(3):
        for j in range(3): perm[3*j+i]=3*i+j; perm[9+3*j+i]=9+3*i+j
    bo=np.abs(b).max(); V=np.abs(b).sum()
    best=None
    for rh in np.linspace(0.01,0.6,591):
        val=rh**-DEG*trunc_tail(P,bo+3*rh)
        if best is None or val<best[0]: best=(val,rh)
    rho=best[1]; T=trunc_tail(P,bo+3*rho)
    ehalf=math.exp(V/2)*(1+1e-12)
    w=np.zeros(R.N,complex); ew=np.zeros(R.N)
    for k,mo in enumerate(R.mons):
        e=[0]*18
        for v,x in enumerate(mo): e[perm[v]]+=x
        kk=R.idx[tuple(e)]
        if not ZW[kk]: continue
        sc=2.0**R.degs[k]*R.fact[k]
        w[kk]=Z.v[k]*sc
        au=mo[:9]; ab_=mo[9:]
        def eps(alpha):
            d=sum(alpha); return 2.0**d*float(np.prod([math.factorial(x) for x in alpha]))*rho**(-d)*T
        ea=eps(au); eb=eps(ab_)
        ew[kk]=sc*errZ[k]+2*U_*sc*abs(Z.v[k])+ea*ehalf+(ehalf+ea)*eb
    return w,ew*(1+1e-10)
def weight(m):
    w=[0,0,0]
    for v,x in enumerate(m):
        if x==0: continue
        if v<9: i,j=divmod(v,3); s=1
        else: i,j=divmod(v-9,3); s=-1
        w[i]+=s*x; w[j]-=s*x
    return tuple(w)
ZW=None
def ring(DEG=6):
    global ZW
    R=FRing(DEG); ZW=np.array([weight(m)==(0,0,0) for m in R.mons]); return R
if __name__=='__main__':
    import time, sys
    R=ring(6)
    from series import Zseries, moments
    for b in [[1.5,1.125,1.125],[1.5,1.5,1.5],[0.9,0.75,-0.3+0.5j]]:
        t=time.time(); w,ew=moments_rig(R,b); dt=time.time()-t
        w2=moments(R,Zseries(R,np.asarray(b,complex),P=34,KMAX=16))
        nz=ZW
        print(b,f"{dt:.1f}s  max err bound {ew.max():.3e}  max actual diff {np.abs(w-w2)[nz].max():.3e}  max nonzero-weight float {np.abs(w2[~nz]).max():.2e}  Z0 {w[0].real:.12f}",flush=True)
