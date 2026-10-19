"""Exact (Gaussian rational) coefficients of phi_h, G_{c,h} in the monomials of u_{ij} (var 3i+j) and
ubar_{ij} (var 9+3i+j), for h in the real basis E_ij, iE_ij of M_3(C) and c in the diagonal units."""
from fractions import Fraction as Fr
import numpy as np
from series import Ring
def gmul(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def gadd(a,b): return (a[0]+b[0],a[1]+b[1])
def gconj(a): return (a[0],-a[1])
ZERO=(Fr(0),Fr(0))
def padd(P,Q):
    out=dict(P)
    for k,v in Q.items(): out[k]=gadd(out.get(k,ZERO),v)
    return {k:v for k,v in out.items() if v!=ZERO}
def pscale(P,c): return {k:gmul(v,c) for k,v in P.items() if gmul(v,c)!=ZERO}
def pmul(P,Q):
    out={}
    for a,x in P.items():
        for b,y in Q.items():
            k=tuple(i+j for i,j in zip(a,b)); out[k]=gadd(out.get(k,ZERO),gmul(x,y))
    return {k:v for k,v in out.items() if v!=ZERO}
def mono(v):
    e=[0]*18; e[v]=1; return tuple(e)
ONE=(Fr(1),Fr(0)); I=(Fr(0),Fr(1))
def u(i,j): return {mono(3*i+j):ONE}
def ub(i,j): return {mono(9+3*i+j):ONE}
def cmat(M):  # complex numpy 3x3 with entries in {0, 1, i, ...} -> Gaussian rationals
    return [[(Fr(M[i][j].real).limit_denominator(),Fr(M[i][j].imag).limit_denominator()) for j in range(3)] for i in range(3)]
def lin(c,Ufun):  # (cU)_{ab}
    return [[ (lambda a,b: (lambda acc: acc)(sum_list([pscale(Ufun(k,b),c[a][k]) for k in range(3)])))(a,b) for b in range(3)] for a in range(3)]
def sum_list(L):
    out={}
    for P in L: out=padd(out,P)
    return out
def tr(A): return sum_list([A[0][0],A[1][1],A[2][2]])
def conjm(c): return [[gconj(x) for x in row] for row in c]
def phi(h):  # Re Tr(hU) = (Tr(hU) + Tr(conj(h) Ubar))/2
    return pscale(padd(tr(lin(h,u)),tr(lin(conjm(h),ub))),(Fr(1,2),Fr(0)))
def G(c,h):
    A=lin(c,u);B=lin(h,u); T=sum_list([pmul(A[a][k],B[k][a]) for a in range(3) for k in range(3)])
    A=lin(conjm(c),ub);B=lin(conjm(h),ub); Tb=sum_list([pmul(A[a][k],B[k][a]) for a in range(3) for k in range(3)])
    imc=pscale(padd(tr(lin(c,u)),pscale(tr(lin(conjm(c),ub)),(Fr(-1),Fr(0)))),(Fr(0),Fr(-1,2)))
    imh=pscale(padd(tr(lin(h,u)),pscale(tr(lin(conjm(h),ub)),(Fr(-1),Fr(0)))),(Fr(0),Fr(-1,2)))
    # 1/2 Re Tr(c h*)
    const=Fr(0)
    for i in range(3):
        for j in range(3):
            const+=gmul(c[i][j],gconj(h[i][j]))[0]
    out={tuple([0]*18):(const/2,Fr(0))} if const!=0 else {}
    out=padd(out,pscale(padd(T,Tb),(Fr(-1,4),Fr(0))))
    out=padd(out,pscale(pmul(imc,imh),(Fr(-1,3),Fr(0))))
    return out
E=[]
for i in range(3):
    for j in range(3):
        M=np.zeros((3,3),complex); M[i,j]=1; E.append(M); M=np.zeros((3,3),complex); M[i,j]=1j; E.append(M)
FB=[np.diag([1,0,0]).astype(complex),np.diag([0,1,0]).astype(complex),np.diag([0,0,1]).astype(complex),np.diag([0,0,1j])]
def tovec(R,P):
    v=np.zeros(R.N,complex); e=np.zeros(R.N)
    for k,(re,im) in P.items():
        x=complex(float(re),float(im)); v[R.idx[k]]=x
        e[R.idx[k]]=abs(x)*2**-52   # one rounding per component
    return v,e
def build(R2):
    Eg=[cmat(M) for M in E]; Fg=[cmat(M) for M in FB]
    PHI=[tovec(R2,phi(h)) for h in Eg]
    GG=[[tovec(R2,G(Eg[p],Eg[q])) for q in range(18)] for p in range(18)]
    GAM=[[tovec(R2,G(f,h)) for h in Eg] for f in Fg]
    A=lambda L: (np.array([x[0] for x in L]),np.array([x[1] for x in L]))
    PHIv,PHIe=A(PHI)
    GGv=np.array([[x[0] for x in row] for row in GG]); GGe=np.array([[x[1] for x in row] for row in GG])
    GAMv=np.array([[x[0] for x in row] for row in GAM]); GAMe=np.array([[x[1] for x in row] for row in GAM])
    return PHIv,PHIe,GGv,GGe,GAMv,GAMe
def save(path):
    PHIv,PHIe,GGv,GGe,GAMv,GAMe=build(Ring(2))
    np.savez(path,PHIv=PHIv,PHIe=PHIe,GGv=GGv,GGe=GGe,GAMv=GAMv,GAMe=GAMe)
if __name__=='__main__':
    import os; save(os.path.join(os.path.dirname(os.path.abspath(__file__)),'epoly.npz'))
