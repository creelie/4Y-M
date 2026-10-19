"""Exact polynomial data for the degree-2 corrector. su(3) is spanned by A_ij = E_ij - E_ji, S_ij = i(E_ij + E_ji)
(i < j), H1 = i(E11 - E22) and H3 = i(E11 + E22 - 2 E33); with weights 1/2 (off-diagonal, H1) and 1/6 (H3) the sum
sum_a w_a X_a (x) X_a is the Casimir tensor of the metric -tr(XY), so |grad f|^2 = sum_a w_a (X_a f)^2 and
Delta = sum_a w_a X_a^2. All entries are Gaussian integers, so the scaled polynomials below have Gaussian
integer coefficients; they are computed in double precision with a bound on every partial sum below 2^50,
hence exactly. U_ij <-> variable 3i+j, Ubar_ij <-> 9+3i+j, left-invariant fields X f(U) = d/dt f(U e^{tX})."""
import numpy as np, scipy.sparse as sp
import rs6z
from series import NV
R,SR=rs6z.rings()
degs=R.degs; N3=int(np.searchsorted(degs,4)); N2=int(np.searchsorted(degs,3)); N1=int(np.searchsorted(degs,2))
mons=R.mons; idx=R.idx
P3=np.array([[idx[tuple(x+y for x,y in zip(mons[i],mons[j]))] for j in range(N3)] for i in range(N3)],np.int32)
def generators():
    G=[]
    for (i,j) in [(0,1),(0,2),(1,2)]:
        A=np.zeros((3,3),complex);A[i,j]=1;A[j,i]=-1;G.append(A)
        S=np.zeros((3,3),complex);S[i,j]=1j;S[j,i]=1j;G.append(S)
    G.append(np.diag([1j,-1j,0])); G.append(np.diag([1j,1j,-2j]))
    return G
GEN=generators(); W6=np.array([3,3,3,3,3,3,3,1],float)   # 6 w_a
def derivation(X):
    rows=[];cols=[];vals=[]
    for k in range(N3):
        m=mons[k]
        for v,x in enumerate(m):
            if x==0: continue
            base=list(m); base[v]-=1
            if v<9: i,j=divmod(v,3); terms=[(3*i+kk, X[kk,j]) for kk in range(3)]
            else: i,j=divmod(v-9,3); terms=[(9+3*i+kk, np.conj(X[kk,j])) for kk in range(3)]
            for (u,c) in terms:
                if c==0: continue
                e=list(base); e[u]+=1; rows.append(idx[tuple(e)]); cols.append(k); vals.append(x*c)
    return sp.csr_matrix((vals,(rows,cols)),shape=(N3,N3),dtype=complex)
DX=[derivation(X) for X in GEN]; ADX=[abs(D) for D in DX]
CJ=np.array([idx[tuple(m[9:]+m[:9])] for m in mons[:N3]])
def lin2(h):
    """2 phi_h, phi_h = Re Tr(U h) = sum_ij (U_ij h_ji + conj)/2"""
    f=np.zeros(N3,complex)
    for i in range(3):
        for j in range(3):
            e=[0]*NV; e[3*i+j]=1; f[idx[tuple(e)]]+=h[j,i]
            e=[0]*NV; e[9+3*i+j]=1; f[idx[tuple(e)]]+=np.conj(h[j,i])
    return f
EH=[]
for i in range(3):
    for j in range(3):
        M=np.zeros((3,3),complex); M[i,j]=1; EH.append(M); M=np.zeros((3,3),complex); M[i,j]=1j; EH.append(M)
PHI2=np.array([lin2(h) for h in EH])                       # 18 x N3, 2 phi_{E_m}
FB=[np.diag([1,0,0]),np.diag([0,1,0]),np.diag([0,0,1]),np.diag([0,0,1j])]
PHB2=[lin2(np.asarray(h,complex)) for h in FB]              # V = sum_t beta_t phi_{FB_t}, beta = (a1, a2, x3, y3)
def mulmat(f):
    """sparse matrix of g -> f g on polys of degree <= 3 (f of degree <= 1, g of degree <= 2)"""
    rows=[];cols=[];vals=[]
    nz=np.nonzero(f)[0]
    for i in nz:
        for j in range(N2):
            k=P3[i,j]
            if k<N3: rows.append(k);cols.append(j);vals.append(f[i])
    return sp.csr_matrix((vals,(rows,cols)),shape=(N3,N3),dtype=complex)
# 6 Gamma_t = sum_a 6 w_a (X_a 2phi_t)(X_a .) : equals 12 * (grad phi_t . grad .)
GT6=[];AGT6=[]
for f in PHB2:
    M=None;A=None
    for D,AD,w in zip(DX,ADX,W6):
        Mf=mulmat(D@f); t=w*(Mf@D); ta=w*(abs(Mf)@AD)
        M=t if M is None else M+t; A=ta if A is None else A+ta
    GT6.append(M); AGT6.append(A)
def lap6(f,fa):
    """6 Delta f and an entrywise bound of the partial sums"""
    v=sum(w*(D@(D@f)) for D,w in zip(DX,W6)); a=sum(w*(AD@(AD@fa)) for AD,w in zip(ADX,W6)); return v,a
# real basis of degree 1..2: 2B = e + ebar, -i(e - ebar); self-conjugate monomials 2e
B2=[]
for k in range(1,N2):
    c=CJ[k]
    if c==k: e=np.zeros(N3,complex); e[k]=2; B2.append(e)
    elif k<c:
        e=np.zeros(N3,complex); e[k]=1; e[c]=1; B2.append(e)
        e=np.zeros(N3,complex); e[k]=-1j; e[c]=1j; B2.append(e)
B2=np.array(B2); NB=len(B2)
LAPB=np.array([lap6(b,np.abs(b))[0] for b in B2])/12        # Delta B
GTB=[np.array([G@b for b in B2])/24 for G in GT6]          # grad phi_t . grad B
DB=[np.array([D@b for b in B2])/2 for D in DX]             # X_a B (degree <= 2)
EXACT=2.0**50
def scaled(Cint):
    """Cint: NB x 18 integers (C = Cint / 2^30). Returns integer-valued arrays (18 x N3 complex):
    Qt = 2^31 q, G_a = X_a Qt, F0 = 12 s phi + 12 s Delta q, F_t = 12 s grad phi_t . grad q, s = 2^31."""
    Cint=np.asarray(Cint,float); assert np.all(Cint==np.round(Cint))
    Qt=Cint.T@B2; Qa=np.abs(Cint.T)@np.abs(B2)
    s=2.0**31
    G=[];bound=Qa.max()
    for D,AD in zip(DX,ADX):
        g=(D@Qt.T).T; ga=(AD@Qa.T).T; G.append(g); bound=max(bound,ga.max())
    L6=[];La=[]
    for q,qa in zip(Qt,Qa):
        v,a=lap6(q,qa); L6.append(v); La.append(a)
    L6=np.array(L6); La=np.array(La)
    F0=6*s*PHI2+2*L6; F0a=6*s*np.abs(PHI2)+2*La
    Ft=[];Fta=[]
    for M,A in zip(GT6,AGT6):
        Ft.append((M@Qt.T).T); Fta.append((A@Qa.T).T)
    bound=max(bound,La.max(),F0a.max(),max(x.max() for x in Fta))
    assert bound<EXACT, bound
    for arr in [Qt,F0]+G+Ft:
        assert np.all(arr.real==np.round(arr.real)) and np.all(arr.imag==np.round(arr.imag))
    return dict(Qt=Qt,G=G,F0=F0,Ft=Ft,s=s,bound=bound)
if __name__=='__main__':
    # checks: Delta phi = -(8/3) phi; grad phi_t . grad phi_h against the formula with orthonormal Gell-Mann basis
    f=PHI2[0]; v,_=lap6(f,np.abs(f)); print('Delta phi/phi', (v/6)[np.nonzero(f)[0]]/f[np.nonzero(f)[0]])
    print('NB',NB,'N3',N3)
