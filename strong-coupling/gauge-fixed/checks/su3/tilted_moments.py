"""Exact tilted Haar moments on SU(3).

For a sequence X = (X_1, ..., X_n) with X_s in {U, U*} (encoded +1, -1) the tensor
E_H[X_1 (x) ... (x) X_n exp(r Re Tr U)] commutes with V^{(x)n} for V in U(3) (conjugation invariance),
so it is a combination of permutation operators (Schur-Weyl). Its coefficients are obtained from the
traces E_H[prod_cycles Tr U^{m_c} exp(r Re Tr U)], which are class functions; by the Weyl integration
formula and exp(r cos t) = sum_k I_k(r) e^{ikt} each is an exact rational combination of the sums
S(n1,n2,n3) = sum_k I_{k-n1}(r) I_{k-n2}(r) I_{k-n3}(r). Every quantity below is therefore an exact
rational (Gaussian rational) combination of class integrals, evaluated at the end in ball arithmetic."""
from fractions import Fraction as Fr
import itertools
from functools import lru_cache
import flint

# ---------- Gaussian rationals ----------
def gq(re=0,im=0): return (Fr(re),Fr(im))
def gadd(a,b): return (a[0]+b[0],a[1]+b[1])
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gconj(a): return (a[0],-a[1])
def gscale(a,c): c=Fr(c); return (a[0]*c,a[1]*c)
ZERO=gq(0); ONE=gq(1); I_=gq(0,1)

# ---------- Laurent polynomials in z1,z2,z3 ----------
def lmul(A,B):
    C={}
    for a,x in A.items():
        for b,y in B.items():
            k=(a[0]+b[0],a[1]+b[1],a[2]+b[2]); C[k]=C.get(k,0)+x*y
    return {k:v for k,v in C.items() if v!=0}
def trpow(m): return {(m,0,0):1,(0,m,0):1,(0,0,m):1} if m!=0 else {(0,0,0):3}
VDM={(0,0,0):1}
for i,j in ((0,1),(0,2),(1,2)):
    a=[0,0,0]; a[i]=1; b=[0,0,0]; b[j]=1
    VDM=lmul(VDM,lmul({tuple(a):1,tuple(b):-1},{tuple(-x for x in a):1,tuple(-x for x in b):-1}))
@lru_cache(None)
def class_poly(key):
    """(1/6) * prod Tr U^m * |Delta|^2 as a dict exponent -> Fraction."""
    P=dict(VDM)
    for m in key: P=lmul(P,trpow(m))
    return {k:Fr(v,6) for k,v in P.items()}

# ---------- Schur-Weyl ----------
def cycles(p):
    seen=[False]*len(p); out=[]
    for i in range(len(p)):
        if not seen[i]:
            c=[]; j=i
            while not seen[j]: seen[j]=True; c.append(j); j=p[j]
            out.append(c)
    return out
@lru_cache(None)
def sw(n):
    perms=list(itertools.permutations(range(n)))
    G=flint.fmpq_mat(len(perms),len(perms))
    def comp(p,q): return tuple(p[q[i]] for i in range(n))
    sg=[]
    for p in perms:
        inv=sum(1 for i in range(n) for j in range(i+1,n) if p[i]>p[j]); sg.append(-1 if inv%2 else 1)
    for a,p in enumerate(perms):
        for b,q in enumerate(perms):
            # Tr(P_p P_q) = 3^{#cycles(p o q)}; for n >= 4 add the sign vector to remove the kernel
            G[a,b]=3**len(cycles(comp(p,q)))+(sg[a]*sg[b] if n>=4 else 0)
    Gi=G.inv()
    return perms,[[Fr(int(Gi[a,b].p),int(Gi[a,b].q)) for b in range(len(perms))] for a in range(len(perms))]
@lru_cache(None)
def mcoef(X):
    """coefficients m_sigma of E[X_1(x)...(x)X_n e^{r ReTr U}] as dicts classkey -> Fraction."""
    n=len(X); perms,Gi=sw(n)
    tkeys=[tuple(sorted(sum(X[s] for s in c) for c in cycles(p))) for p in perms]
    out=[]
    for a in range(len(perms)):
        d={}
        for b,k in enumerate(tkeys):
            if Gi[a][b]!=0: d[k]=d.get(k,0)+Gi[a][b]
        out.append({k:v for k,v in d.items() if v!=0})
    return perms,out
@lru_cache(None)
def mentry(X,I,J):
    """E[prod_s (X_s)_{I_s J_s} e^{r ReTr U}] as dict classkey -> Fraction.
    P_p maps e_J to e_I with I_s = J_{p(s)} (the convention of the Gram matrix is immaterial)."""
    if len(X)==0: return {():Fr(1)}
    perms,ms=mcoef(X); d={}
    for p,m in zip(perms,ms):
        if all(I[s]==J[p[s]] for s in range(len(X))):
            for k,v in m.items(): d[k]=d.get(k,0)+v
    return {k:v for k,v in d.items() if v!=0}

# ---------- sparse polynomials in the entries of U and U* ----------
# term: (X, idx, coeff) with idx = (i1,j1,...,in,jn); value coeff * prod (X_s)_{i_s j_s}
def dag(A): return {(j,i):gconj(v) for (i,j),v in A.items()}
def tr1(A,x): return [((x,),(i,a),v) for (a,i),v in A.items()]
def tr2(A,x,B,y): return [((x,y),(i,j,k,a),gmul(va,vb)) for (a,i),va in A.items() for (j,k),vb in B.items()]
def scal(p,c): return [(X,idx,gmul(v,c) if isinstance(c,tuple) else gscale(v,c)) for X,idx,v in p]
def pmul(p,q): return [(X+Y,i1+i2,gmul(a,b)) for X,i1,a in p for Y,i2,b in q]
def const(c): return [((),(),c)]
def retr(A,B):  # Re Tr(A B*)
    s=ZERO
    for k,v in A.items():
        if k in B: s=gadd(s,gmul(v,gconj(B[k])))
    return gq(s[0],0)
def G_poly(W,h):
    """G_{W,h} = <grad phi_W, grad phi_h> = (1/2)ReTr(W h*) - (1/2)ReTr(WUhU) - (1/3)ImTr(WU)ImTr(hU)."""
    p=const(gscale(retr(W,h),Fr(1,2)))
    p+=scal(tr2(W,1,h,1),Fr(-1,4))+scal(tr2(dag(h),-1,dag(W),-1),Fr(-1,4))
    a=tr1(W,1)+scal(tr1(dag(W),-1),-1); b=tr1(h,1)+scal(tr1(dag(h),-1),-1)
    p+=scal(pmul(a,b),Fr(1,12))
    return p
def g_poly(W,h):
    """G_{W,h} minus its Haar mean (4/9)ReTr(W h*)."""
    p=G_poly(W,h)
    return p+const(gscale(retr(W,h),Fr(-4,9)))
def expect(p):
    """E_H[p e^{r ReTr U}] as dict classkey -> Gaussian rational."""
    d={}
    for X,idx,v in p:
        n=len(X)
        I=tuple(idx[2*s] for s in range(n)); J=tuple(idx[2*s+1] for s in range(n))
        for k,c in mentry(X,I,J).items():
            # E[prod Tr U^m e^{r ReTr U}] is invariant under U -> conj(U), i.e. under m -> -m
            k=min(k,tuple(sorted(-x for x in k)))
            d[k]=gadd(d.get(k,ZERO),gscale(v,c))
    return {k:v for k,v in d.items() if v!=ZERO}
def basis18():
    B=[]
    for i in range(3):
        for j in range(3):
            for c in (ONE,I_): B.append({(i,j):c})
    return B
def basis6():
    return [{(k,k):c} for k in range(3) for c in (ONE,I_)]

# ---------- evaluation in ball arithmetic ----------
_S={}
def S(n,r,K=40):
    key=(n,str(r))
    if key not in _S:
        tot=flint.arb(0)
        for k in range(-K,K+1):
            tot+=r.bessel_i(k-n[0])*r.bessel_i(k-n[1])*r.bessel_i(k-n[2])
        # tail |k| > K (valid for 0 <= r <= 2): every order is at least |k| - 12 >= K - 11 = mm, and
        # I_m(r) <= (r/2)^m e^{r^2/4}/m! decreases in m, so with R = upper bound of r the tail is at most
        # 2 * sum_{j >= 0} [(R/2)^(mm+j) e^{R^2/4}/(mm+j)!]^3 <= 4 [(R/2)^mm e^{R^2/4}/mm!]^3
        assert 0<=r and r<=2
        mm=K-11; R=flint.arb(r.upper())
        t=4*((R/2)**mm*(R*R/4).exp()/flint.arb.fac_ui(mm))**3
        tot+=flint.arb(0,t.upper())
        _S[key]=tot
    return _S[key]
_CI={}
def CI(key,r):
    k2=(key,str(r))
    if k2 not in _CI:
        P=class_poly(key); tot=flint.arb(0)
        for n,c in P.items():
            assert max(abs(x) for x in n)<=12
            tot+=flint.arb(c.numerator)/c.denominator*S(n,r)
        _CI[k2]=tot
    return _CI[k2]
def ev(d,r):
    re=flint.arb(0); im=flint.arb(0)
    for k,(a,b) in d.items():
        v=CI(k,r)
        re+=flint.arb(a.numerator)/a.denominator*v; im+=flint.arb(b.numerator)/b.denominator*v
    return re,im
