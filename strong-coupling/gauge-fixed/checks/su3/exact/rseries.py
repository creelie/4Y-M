"""Tilted SU(3) moments at diagonal b with rigorous error bounds.

Z(D, D') = E_H exp(Tr(UD)/2 + Tr(Ubar D')/2) = sum_{p>=q>=0} s_(p,q)(D D'^T)/d_(p,q) G_(p,q)(det D) G_(p,q)(det D'),
truncated to p <= P and powers of det at most K (3(K+1) > P), as a Taylor series of total degree 4 in the
18 entries of D - b, D' - conj(b). Rounding: every series carries the same computation on absolute values
(ab) and an operation count n, with |computed - exact| <= gamma(n) * ab_exact coefficientwise.
Truncation: Cauchy estimates (see trunc_bound)."""
import numpy as np, math
from series import Ring, hookprod, dim
U_=2.0**-53
def gam(n): return n*U_/(1-n*U_)
class TS:
    __slots__=('v','a','n')
    def __init__(s,v,a,n): s.v=v; s.a=a; s.n=n
def mkconst(R,c):
    v=np.zeros(R.N,complex); v[0]=c; a=np.zeros(R.N); a[0]=abs(c); return TS(v,a,0)
def mkvar(R,k):
    v=np.zeros(R.N,complex); v[R.idx[tuple(1 if i==k else 0 for i in range(18))]]=1; return TS(v,np.abs(v),0)
def add(x,y): return TS(x.v+y.v,x.a+y.a,max(x.n,y.n)+2)
def sub(x,y): return TS(x.v-y.v,x.a+y.a,max(x.n,y.n)+2)
def mulR(R,a,b): return np.bincount(R.K,a[R.I]*b[R.J],R.N)
def mul(R,x,y): return TS(R.mul(x.v,y.v),mulR(R,x.a,y.a),x.n+y.n+40)
def scal(c,x,nc=0): return TS(c*x.v,abs(c)*x.a,x.n+nc+4)
def divint(x,d): return TS(x.v/d,x.a/d,x.n+2)
def det3(R,M):
    m=lambda a,b: mul(R,a,b)
    t1=m(M[0][0],sub(m(M[1][1],M[2][2]),m(M[1][2],M[2][1])))
    t2=m(M[0][1],sub(m(M[1][0],M[2][2]),m(M[1][2],M[2][0])))
    t3=m(M[0][2],sub(m(M[1][0],M[2][1]),m(M[1][1],M[2][0])))
    return add(sub(t1,t2),t3)
def clam(l): return 1.0/(2.0**sum(l)*float(hookprod(l)))      # relative error <= gamma(2)
def Zseries_rig(R,b,P=24,K=8):
    assert 3*(K+1)>P
    D=[[add(mkconst(R,b[j] if i==j else 0),mkvar(R,3*j+i)) for i in range(3)] for j in range(3)]
    Dp=[[add(mkconst(R,np.conj(b[j]) if i==j else 0),mkvar(R,9+3*j+i)) for i in range(3)] for j in range(3)]
    def ssum(L):
        out=L[0]
        for x in L[1:]: out=add(out,x)
        return out
    X=[[ssum([mul(R,D[a][k],Dp[c][k]) for k in range(3)]) for c in range(3)] for a in range(3)]
    m=lambda a,b: mul(R,a,b)
    e1=ssum([X[0][0],X[1][1],X[2][2]])
    e2=ssum([sub(m(X[0][0],X[1][1]),m(X[0][1],X[1][0])),sub(m(X[0][0],X[2][2]),m(X[0][2],X[2][0])),sub(m(X[1][1],X[2][2]),m(X[1][2],X[2][1]))])
    e3=det3(R,X)
    h=[mkconst(R,1.0),e1,sub(m(e1,e1),e2)]
    for k in range(3,P+2): h.append(add(sub(m(e1,h[k-1]),m(e2,h[k-2])),m(e3,h[k-3])))
    dD=det3(R,D); dDp=det3(R,Dp)
    pw=[mkconst(R,1.0)]; pwp=[mkconst(R,1.0)]
    for k in range(1,K+1): pw.append(m(pw[-1],dD)); pwp.append(m(pwp[-1],dDp))
    Z=None
    for p in range(P+1):
        for q in range(p+1):
            s_=h[p] if q==0 else sub(m(h[p],h[q]),m(h[p+1],h[q-1]))
            cs=[clam((p+k,q+k,k)) for k in range(K+1)]
            G=ssum([scal(c,x,2) for c,x in zip(cs,pw)]); Gp=ssum([scal(c,x,2) for c,x in zip(cs,pwp)])
            t=divint(m(s_,m(G,Gp)),dim(p,q))
            Z=t if Z is None else add(Z,t)
    return Z
def trunc_tail(P,x,terms=200):
    """upper bound for sum_{n>P} (3x/2)^n/n! (ratio test)"""
    y=1.5*x; t=y**(P+1)/math.factorial(P+1); s=0.0; n=P+1
    while True:
        s+=t; n+=1; t*=y/n
        if t<1e-30*s and y/n<0.5: break
    return s*(1+1e-12)+t*2
def moments_rig(R,b,P=24,K=8):
    """Returns (w, ew): w[k] approximates E[u^alpha ubar^beta e^{Re Tr bU}] for monomial k (u_{ij} <-> var 3i+j,
    ubar_{ij} <-> 9+3i+j), and ew[k] bounds the absolute error."""
    b=np.asarray(b,complex)
    Z=Zseries_rig(R,b,P,K)
    g=gam(Z.n); errZ=g*Z.a/(1-g)
    perm=[0]*18
    for i in range(3):
        for j in range(3): perm[3*j+i]=3*i+j; perm[9+3*j+i]=9+3*i+j
    bo=np.abs(b).max(); V=np.abs(b).sum()
    rho=min(0.25,max(0.02,2.0*(1.5*bo)/(1.5*3*(P+1))))
    # choose rho minimising rho^-4 (3(bo+3rho)/2)^(P+1)
    best=None
    for rh in np.linspace(0.01,0.5,491):
        val=rh**-4*trunc_tail(P,bo+3*rh)
        if best is None or val<best[0]: best=(val,rh)
    rho=best[1]; T=trunc_tail(P,bo+3*rho)
    ehalf=math.exp(V/2)*(1+1e-12)
    w=np.zeros(R.N,complex); ew=np.zeros(R.N)
    for k,mo in enumerate(R.mons):
        e=[0]*18
        for v,x in enumerate(mo): e[perm[v]]+=x
        kk=R.idx[tuple(e)]
        sc=2.0**R.degs[k]*R.fact[k]
        w[kk]=Z.v[k]*sc
        au=mo[:9]; ab_=mo[9:]
        def eps(alpha):
            d=sum(alpha); return 2.0**d*float(np.prod([math.factorial(x) for x in alpha]))*rho**(-d)*T
        ea=eps(au); eb=eps(ab_)
        ew[kk]=sc*errZ[k]+2*U_*sc*abs(Z.v[k])+ea*ehalf+(ehalf+ea)*eb
    return w,ew*(1+1e-10)
if __name__=='__main__':
    import time
    R=Ring(4)
    for b in [[1.36,1.36,1.36],[1.36,1.02,1.02],[0.5,0.25,-0.375+0.5j]]:
        t=time.time(); w,ew=moments_rig(R,b); dt=time.time()-t
        from series import Zseries, moments
        w2=moments(R,Zseries(R,np.asarray(b,complex),P=34,KMAX=16))
        print(b,f"{dt:.1f}s  max err bound {ew.max():.3e}  max actual diff {np.abs(w-w2).max():.3e}  Z0 {w[0].real:.12f}")
