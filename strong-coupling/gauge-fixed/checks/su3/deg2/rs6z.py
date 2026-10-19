"""Rigorous tilted SU(3) moments of degree <= 6 (rseries6.py) computed in the subring of torus-weight-zero
monomials. e1, e2, e3, det D, det D' are invariant under conjugation by the diagonal torus, so their exact
coefficients vanish off weight zero; they are formed in the full ring and projected, and every later product
is a product of invariant elements, whose weight-zero coefficients involve only weight-zero factors."""
import numpy as np, math, sys, time
from ring6 import FRing
from series import hookprod, dim
from rseries import gam, U_, trunc_tail
from cmul import cmul
import rseries6 as R6m
from rseries6 import TS, mkconst, mkvar, add, sub, det3, clam, weight, NMUL
class SubRing:
    def __init__(s,R):
        zw=np.array([weight(m)==(0,0,0) for m in R.mons]); s.zw=zw
        s.full=np.nonzero(zw)[0]; s.N=len(s.full)
        loc=-np.ones(R.N,np.int64); loc[s.full]=np.arange(s.N); s.loc=loc
        keep=zw[R.I]&zw[R.J]
        s.I=loc[R.I[keep]].astype(np.int32); s.J=loc[R.J[keep]].astype(np.int32); s.K=loc[R.K[keep]].astype(np.int32)
        assert (s.K>=0).all()
        s.mons=[R.mons[i] for i in s.full]; s.degs=R.degs[s.full]; s.fact=R.fact[s.full]
def proj(S,x): return TS(x.v[S.full].copy(),x.a[S.full].copy(),x.n)
def smul(S,x,y):
    v,a=cmul(S,x.v,y.v,x.a,y.a); return TS(v,a,x.n+y.n+NMUL)
def sconst(S,c):
    v=np.zeros(S.N,complex); v[0]=c; a=np.zeros(S.N); a[0]=abs(c); return TS(v,a,0)
def Zseries_rig(R,S,b,P=34,K=11,CH=64):
    assert 3*(K+1)>P
    D=[[add(mkconst(R,b[j] if i==j else 0),mkvar(R,3*j+i)) for i in range(3)] for j in range(3)]
    Dp=[[add(mkconst(R,np.conj(b[j]) if i==j else 0),mkvar(R,9+3*j+i)) for i in range(3)] for j in range(3)]
    def ssum(L):
        out=L[0]
        for x in L[1:]: out=add(out,x)
        return out
    m=lambda a,b: R6m.mul(R,a,b)
    X=[[ssum([m(D[a][k],Dp[c][k]) for k in range(3)]) for c in range(3)] for a in range(3)]
    e1=ssum([X[0][0],X[1][1],X[2][2]])
    e2=ssum([sub(m(X[0][0],X[1][1]),m(X[0][1],X[1][0])),sub(m(X[0][0],X[2][2]),m(X[0][2],X[2][0])),sub(m(X[1][1],X[2][2]),m(X[1][2],X[2][1]))])
    e3=det3(R,X); dD=det3(R,D); dDp=det3(R,Dp)
    e1,e2,e3,dD,dDp=[proj(S,x) for x in (e1,e2,e3,dD,dDp)]
    ms=lambda a,b: smul(S,a,b)
    h=[sconst(S,1.0),e1,sub(ms(e1,e1),e2)]
    for k in range(3,P+2): h.append(add(sub(ms(e1,h[k-1]),ms(e2,h[k-2])),ms(e3,h[k-3])))
    prods={}
    def hh(i,j):
        key=(min(i,j),max(i,j))
        if key not in prods: prods[key]=h[key[1]] if key[0]==0 else ms(h[key[0]],h[key[1]])
        return prods[key]
    pairs=[(p,q) for p in range(P+1) for q in range(p+1)]
    Sv=np.empty((len(pairs),S.N),complex); Sa=np.empty((len(pairs),S.N)); nS=0
    for i,(p,q) in enumerate(pairs):
        x=hh(p,0) if q==0 else sub(hh(p,q),hh(p+1,q-1)); Sv[i]=x.v; Sa[i]=x.a; nS=max(nS,x.n)
    KL=[(k,l) for k in range(K+1) for l in range(K+1)]
    CO=np.array([[clam((p+k,q+k,k))*clam((p+l,q+l,l))/dim(p,q) for p,q in pairs] for k,l in KL])
    SKv=CO@Sv; SKa=np.abs(CO)@Sa; nL=nS+8+2*len(pairs)
    pw=[sconst(S,1.0)]; pwp=[sconst(S,1.0)]
    for k in range(1,K+1): pw.append(ms(pw[-1],dD)); pwp.append(ms(pwp[-1],dDp))
    Z=None
    for k in range(K+1):
        inner=None
        for l in range(K+1):
            r=KL.index((k,l)); Skl=TS(SKv[r],SKa[r],nL)
            t=Skl if l==0 else ms(pwp[l],Skl)
            inner=t if inner is None else add(inner,t)
        t=inner if k==0 else ms(pw[k],inner)
        Z=t if Z is None else add(Z,t)
    return Z
_R=None; _S=None
def rings():
    global _R,_S
    if _R is None: _R=R6m.ring(6); _S=SubRing(_R)
    return _R,_S
def moments_rig(b,P=34,K=11):
    """(w, ew) on the full degree-6 index set; nonzero-weight monomials exactly 0 with error 0."""
    R,S=rings(); b=np.asarray(b,complex); DEG=6
    Z=Zseries_rig(R,S,b,P,K)
    g=gam(Z.n); errZ=g*Z.a/(1-g)
    perm=[0]*18
    for i in range(3):
        for j in range(3): perm[3*j+i]=3*i+j; perm[9+3*j+i]=9+3*i+j
    bo=np.abs(b).max(); V=np.abs(b).sum()
    best=None
    for rh in np.linspace(0.01,0.8,791):
        val=rh**-DEG*trunc_tail(P,bo+3*rh)
        if best is None or val<best[0]: best=(val,rh)
    rho=best[1]; T=trunc_tail(P,bo+3*rho)
    ehalf=math.exp(V/2)*(1+1e-12)
    w=np.zeros(R.N,complex); ew=np.zeros(R.N)
    for k,mo in enumerate(S.mons):
        e=[0]*18
        for v,x in enumerate(mo): e[perm[v]]+=x
        kk=R.idx[tuple(e)]
        sc=2.0**S.degs[k]*S.fact[k]
        w[kk]=Z.v[k]*sc
        au=mo[:9]; ab_=mo[9:]
        def eps(alpha):
            d=sum(alpha); return 2.0**d*float(np.prod([math.factorial(x) for x in alpha]))*rho**(-d)*T
        ea=eps(au); eb=eps(ab_)
        ew[kk]=sc*errZ[k]+2*U_*sc*abs(Z.v[k])+ea*ehalf+(ehalf+ea)*eb
    return w,ew*(1+1e-10)
if __name__=='__main__':
    t=time.time(); R,S=rings(); print('rings',S.N,len(S.I),time.time()-t,flush=True)
    for b in [[1.5,1.125,1.125],[0.9,0.75,-0.3+0.5j]]:
        t=time.time(); w,ew=moments_rig(b); dt=time.time()-t
        w1,ew1=R6m.moments_rig(R,b)
        print(b,f"{dt:.1f}s  max err {ew.max():.3e}  vs P=26 version diff {np.abs(w-w1).max():.3e} (its bound {ew1.max():.2e})  Z0 {w[0].real:.15f}",flush=True)
