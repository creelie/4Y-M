"""Exact tilted moments of SU(3) Haar measure at a general (diagonal) b by the character expansion
Z(D,D') = E_H exp(Tr(UD)/2 + Tr(Ubar D')/2) = sum_{p>=q>=0} s_(p,q)(D D'^T)/d_(p,q) G_(p,q)(det D) G_(p,q)(det D'),
G_(p,q)(z) = sum_k c_(p+k,q+k,k) z^k, c_lambda = 2^-|lambda| / prod hooks.
Truncated multivariate Taylor series in the 18 entries of D - b, D' - conj(b), total degree <= DEG."""
import numpy as np, itertools, math
from functools import lru_cache
NV=18
class Ring:
    def __init__(s,deg=4):
        s.deg=deg
        mons=[]
        for d in range(deg+1):
            for c in itertools.combinations_with_replacement(range(NV),d):
                e=[0]*NV
                for v in c: e[v]+=1
                mons.append(tuple(e))
        s.mons=mons; s.N=len(mons); s.idx={m:i for i,m in enumerate(mons)}
        s.degs=np.array([sum(m) for m in mons])
        I=[];J=[];K=[]
        for i,a in enumerate(mons):
            for j,b in enumerate(mons):
                if s.degs[i]+s.degs[j]<=deg:
                    I.append(i);J.append(j);K.append(s.idx[tuple(x+y for x,y in zip(a,b))])
        s.I=np.array(I);s.J=np.array(J);s.K=np.array(K)
        s.fact=np.array([float(np.prod([math.factorial(x) for x in m])) for m in mons])
    def mul(s,a,b):
        w=a[s.I]*b[s.J]
        return np.bincount(s.K,w.real,s.N)+1j*np.bincount(s.K,w.imag,s.N)
    def const(s,c):
        a=np.zeros(s.N,complex); a[0]=c; return a
    def var(s,v,c=1.0):
        a=np.zeros(s.N,complex); e=[0]*NV; e[v]=1; a[s.idx[tuple(e)]]=c; return a
def hookprod(l):
    l1,l2,l3=l
    return math.factorial(l1+2)*math.factorial(l2+1)*math.factorial(l3)//((l1-l2+1)*(l1-l3+2)*(l2-l3+1))
def clam(l): return 1.0/(2.0**sum(l)*hookprod(l))
def dim(p,q): return (p-q+1)*(q+1)*(p+2)//2
def det3(R,M):
    m=lambda a,b: R.mul(a,b)
    return (m(M[0][0],m(M[1][1],M[2][2])-m(M[1][2],M[2][1]))-m(M[0][1],m(M[1][0],M[2][2])-m(M[1][2],M[2][0]))
            +m(M[0][2],m(M[1][0],M[2][1])-m(M[1][1],M[2][0])))
def Zseries(R,b,P=30,KMAX=14):
    """b: diagonal entries (3 complex). Variable v=3*j+i is D_{ji} (pairs with u_{ij}); 9+v is D'_{ji} (pairs with ubar_{ij})."""
    D=[[R.const(b[j] if i==j else 0)+R.var(3*j+i) for i in range(3)] for j in range(3)]
    Dp=[[R.const(np.conj(b[j]) if i==j else 0)+R.var(9+3*j+i) for i in range(3)] for j in range(3)]
    X=[[sum(R.mul(D[a][k],Dp[c][k]) for k in range(3)) for c in range(3)] for a in range(3)]
    e1=X[0][0]+X[1][1]+X[2][2]
    m=R.mul
    e2=(m(X[0][0],X[1][1])-m(X[0][1],X[1][0]))+(m(X[0][0],X[2][2])-m(X[0][2],X[2][0]))+(m(X[1][1],X[2][2])-m(X[1][2],X[2][1]))
    e3=det3(R,X)
    h=[R.const(1.0),e1]
    h.append(m(e1,h[1])-e2)
    for k in range(3,P+2): h.append(m(e1,h[k-1])-m(e2,h[k-2])+m(e3,h[k-3]))
    dD=det3(R,D); dDp=det3(R,Dp)
    pw=[R.const(1.0)]; pwp=[R.const(1.0)]
    for k in range(1,KMAX+1): pw.append(m(pw[-1],dD)); pwp.append(m(pwp[-1],dDp))
    Z=np.zeros(R.N,complex)
    for p in range(P+1):
        for q in range(p+1):
            s_=h[p] if q==0 else m(h[p],h[q])-m(h[p+1],h[q-1])
            cs=[clam((p+k,q+k,k)) for k in range(KMAX+1)]
            G=sum(c*x for c,x in zip(cs,pw)); Gp=sum(c*x for c,x in zip(cs,pwp))
            Z+=m(s_,m(G,Gp))/dim(p,q)
    return Z
def moments(R,Z):
    """L(f) = sum_alpha f_alpha E[u^alpha e^{Re Tr bU}] with the monomial basis of the same ring,
    variable v=3*i+j <-> u_{ij}, 9+v <-> ubar_{ij}. Returns the weight vector w with L(f) = f . w."""
    # coefficient of D_{ji} monomials -> u_{ij} moments: remap variable index 3j+i -> 3i+j
    perm=[0]*NV
    for i in range(3):
        for j in range(3): perm[3*j+i]=3*i+j; perm[9+3*j+i]=9+3*i+j
    w=np.zeros(R.N,complex)
    for k,mo in enumerate(R.mons):
        e=[0]*NV
        for v,x in enumerate(mo): e[perm[v]]+=x
        w[R.idx[tuple(e)]]=Z[k]*2.0**R.degs[k]*R.fact[k]
    return w
