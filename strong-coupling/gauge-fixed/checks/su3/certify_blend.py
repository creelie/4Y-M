"""Certificate for the blended single-link bound of four-dimensional SU(3) in the forest gauge (su3.md, Section 5).

With r = 2 beta_W, for every b in rK (K the convex hull of SU(3)), |h|_HS = 1 and s >= 0:

  (M)  E_H[|grad phi_h|^2 e^{phi_b}] <= mt(r) := max_h E_H[ |grad phi_h|^2 e^{r Re Tr U} ],
  (A)  E_H[|grad phi_h|^2 U] = adj(h)/9 at r = 0, so nu_b|grad phi_h|^2 <= mh := mt/(mt + (5/9) e^{-r/5}),
  (V)  Var_{nu_b}((1 - 8s/3) phi_h + s <grad V, grad phi_h>) <= ks(r),
       ks(r) := max_{W in T, |h|=1} E_H[ ((1 - 8s/3) phi_h + s r g_{W,h})^2 e^{r Re Tr U} ],
       g_{W,h} = <grad phi_W, grad phi_h> - (4/9) Re Tr(W h*),
  (D)  the spectral gap of nu_b is at least gap2(r) = (39/64)(7/6) + 3/2 - r/2 = 283/128 - r/2,

and the Lipschitz rate is at most s sqrt(mh) + sqrt(ks/gap2). This script proves (M), the identity (A)
exactly, (V) for one s, and 2 r Lambda < 1, in ball arithmetic. The maximum over the torus
T = {diag(e^{i t1}, e^{i t2}, e^{-i(t1+t2)})} is covered by dyadic boxes, each checked by a second-order
bound (convexity in the linear term).

Usage: python3 certify_blend.py [r] [s] [MT] [KS]"""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from tilted_moments import *
from fractions import Fraction as Fr
import numpy as np, time
R=sys.argv[1] if len(sys.argv)>1 else '1.14'
SB=sys.argv[2] if len(sys.argv)>2 else '0.322'
MT=sys.argv[3] if len(sys.argv)>3 else '0.8051'
KS=sys.argv[4] if len(sys.argv)>4 else '0.0508'
r=flint.arb(R); s=flint.arb(SB)
c1=1-8*s/3; c2=s*r
t0=time.time()
B18=basis18(); B6=basis6()
BLOCKS=[[0,1,8,9,16,17],[2,3,6,7],[4,5,12,13],[10,11,14,15]]
blockof={p:i for i,Bk in enumerate(BLOCKS) for p in Bk}

def pd(A,t):
    """True if t I - A is certified positive definite (A symmetric list of lists of arb), by LDL^T."""
    n=len(A); M=[[ (t if i==j else 0)-A[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        if not M[k][k]>0: return False
        for i in range(k+1,n):
            f=M[i][k]/M[k][k]
            for j in range(k+1,n): M[i][j]=M[i][j]-f*M[k][j]
    return True
def haar0(d):
    """the exact Haar mean (r = 0) of an expectation dict: S(n, 0) = 1 if n1 = n2 = n3, else 0."""
    tot=ZERO
    for k,v in d.items():
        c=sum((x for n,x in class_poly(k).items() if n[0]==n[1]==n[2]),Fr(0))
        tot=gadd(tot,gscale(v,c))
    return tot

# ---- (0) self-tests at r = 0, exact ----
for p in range(18):
    for q in range(18):
        assert haar0(expect(G_poly(B18[p],B18[q])))==gq(Fr(4,9) if p==q else 0), (p,q)
print("(0) at r = 0 the moment matrix of |grad phi_h|^2 is exactly (4/9) I")
def mat(A): return [[A.get((i,j),ZERO) for j in range(3)] for i in range(3)]
def madd(X,Y): return [[gadd(X[i][j],Y[i][j]) for j in range(3)] for i in range(3)]
def adj(X):
    """adjugate of a 3 x 3 matrix of Gaussian rationals: adj(X)_{ij} = (-1)^{i+j} minor_{ji}."""
    out=[[ZERO]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rr=[k for k in range(3) if k!=j]; cc=[k for k in range(3) if k!=i]
            m=gadd(gmul(X[rr[0]][cc[0]],X[rr[1]][cc[1]]),gscale(gmul(X[rr[0]][cc[1]],X[rr[1]][cc[0]]),-1))
            out[i][j]=gscale(m,(-1)**(i+j))
    return out
for p in range(18):
    for q in range(p,18):
        Bp,Bq=mat(B18[p]),mat(B18[q])
        A2,Ap,Aq=adj(madd(Bp,Bq)),adj(Bp),adj(Bq)
        for i in range(3):
            for j in range(3):
                pol=gscale(gadd(A2[i][j],gscale(gadd(Ap[i][j],Aq[i][j]),-1)),Fr(1,2))
                lhs=haar0(expect(pmul(G_poly(B18[p],B18[q]),[((1,),(i,j),ONE)])))
                assert lhs==gscale(pol,Fr(1,9)), (p,q,i,j,lhs,pol)
print("(A) E_H[|grad phi_h|^2 U] = adj(h)/9 exactly (all 171 x 9 polarised coefficients)")

# ---- (M) the gradient moment ----
Q=[[None]*18 for _ in range(18)]
for p in range(18):
    for q in range(18):
        d=expect(G_poly(B18[p],B18[q]))
        assert all(v[1]==0 for v in d.values()), 'imaginary part'
        Q[p][q]=ev(d,r)[0]
mt_ok=pd(Q,flint.arb(MT))
print(f"(M) max_h E[|grad phi_h|^2 e^(r ReTrU)] < {MT} at r = {R}: {mt_ok}")
qm=np.array([[float(x.mid()) for x in row] for row in Q]); print('    float value',np.linalg.eigvalsh(qm).max())
mt=flint.arb(MT); mh=mt/(mt+flint.arb(5)/9*(-r/5).exp())
print(f"    mh = mt/(mt + (5/9) e^(-r/5)) <= {mh}")

# ---- (V) exact moments of phi_p, g_{d_a,p} ----
phi=[scal(tr1(B18[p],1)+tr1(dag(B18[p]),-1),Fr(1,2)) for p in range(18)]   # Re Tr(e_p U)
g=[[g_poly(B6[a],B18[p]) for p in range(18)] for a in range(6)]
def sym_entry(x,y,p,q):
    d1=expect(pmul(x[p],y[q])); d2=expect(pmul(x[q],y[p]))
    assert all(v[1]==0 for v in d1.values()) and all(v[1]==0 for v in d2.values()), 'imaginary part'
    if blockof[p]!=blockof[q]:
        assert not d1 and not d2, 'off-block entry'
        return None
    return ev(d1,r)[0]+ev(d2,r)[0]
def blockmats(x,y,scale):
    """block matrices of scale * (E[x_p y_q e] + E[x_q y_p e]) (symmetric), exact block structure checked."""
    full=[[None]*18 for _ in range(18)]
    for p in range(18):
        for q in range(18): full[p][q]=sym_entry(x,y,p,q)
    return [flint.arb_mat([[full[p][q]*scale for q in Bk] for p in Bk]) for Bk in BLOCKS]
A0=blockmats(phi,phi,c1*c1/2)                                  # c1^2 E[phi phi e]
A=[blockmats(phi,g[a],c1*c2) for a in range(6)]                 # c1 c2 (PG_a + PG_a^T)
C=[[blockmats(g[a],g[b],c2*c2/2) for b in range(6)] for a in range(6)]   # c2^2 K_ab
A0b=A0; Ab=[[A[a][bi] for a in range(6)] for bi in range(4)]; Cb=[[[C[a][b][bi] for b in range(6)] for a in range(6)] for bi in range(4)]
print(f"    exact moments built, block structure verified ({time.time()-t0:.1f} s)")
def frob_upper(Mt):
    t=flint.arb(0)
    for row in Mt.tolist():
        for x in row: t+=x*x
    return t.sqrt()
# remainder: |d^2/dt^2 w_a| <= 4 rho^2, |d^2/dt^2 w_a w_b| <= 16 rho^2, halved in Taylor's formula
M2=[2*sum((frob_upper(Ab[bi][a]) for a in range(6)),flint.arb(0))
    +8*sum((frob_upper(Cb[bi][a][b]) for a in range(6) for b in range(6)),flint.arb(0)) for bi in range(4)]
M2f=[float(m.upper()) for m in M2]
def tof(Mt): return np.array([[float(x.mid()) for x in row] for row in Mt.tolist()])
A0f=[tof(A0b[bi]) for bi in range(4)]; Abf=[[tof(Ab[bi][a]) for a in range(6)] for bi in range(4)]
Cbf=[[[tof(Cb[bi][a][b]) for b in range(6)] for a in range(6)] for bi in range(4)]
KSf=float(flint.arb(KS).mid())
def wfloat(t1,t2):
    t=(t1,t2,-t1-t2)
    w=np.array([f(tk) for tk in t for f in (np.cos,np.sin)])
    d1=np.array([-np.sin(t[0]),np.cos(t[0]),0,0,np.sin(t[2]),-np.cos(t[2])])
    d2=np.array([0,0,-np.sin(t[1]),np.cos(t[1]),np.sin(t[2]),-np.cos(t[2])])
    return w,d1,d2
def ub_float(x1,x2,rx):
    t1,t2=2*np.pi*float(x1),2*np.pi*float(x2); rho=2*np.pi*float(rx)
    w,d1,d2=wfloat(t1,t2); ub=0
    for bi in range(4):
        Ma=[sum(w[b]*Cbf[bi][a][b] for b in range(6)) for a in range(6)]
        N=A0f[bi]+sum(w[a]*(Abf[bi][a]+Ma[a]) for a in range(6))
        D1=sum(d1[a]*(Abf[bi][a]+2*Ma[a]) for a in range(6)); D2=sum(d2[a]*(Abf[bi][a]+2*Ma[a]) for a in range(6))
        l=max(np.linalg.eigvalsh(N+e1*rho*D1+e2*rho*D2)[-1] for e1 in (-1,1) for e2 in (-1,1))
        ub=max(ub,l+M2f[bi]*rho**2)
    return ub
fl=max(ub_float(Fr(i,90),Fr(j,90),0) for i in range(90) for j in range(90))
print(f"    float maximum over a 90 x 90 grid of the torus: {fl:.6f}")
n0=32
todo=[(Fr(2*i+1,2*n0),Fr(2*j+1,2*n0),Fr(1,2*n0)) for i in range(n0) for j in range(n0)]
boxes=[]
while todo:
    x1,x2,rx=todo.pop()
    if ub_float(x1,x2,rx)<KSf*(1-1e-6): boxes.append((x1,x2,rx)); continue
    assert rx>Fr(1,2**20), 'subdivision too deep'
    hh=rx/2
    todo+=[(x1+e1*hh,x2+e2*hh,hh) for e1 in (-1,1) for e2 in (-1,1)]
area=sum(4*rx*rx for _,_,rx in boxes)
assert area==1, area
print(f"    {len(boxes)} dyadic boxes cover the torus exactly (area {area}); float search {time.time()-t0:.1f} s")
twopi=2*flint.arb.pi(); KSa=flint.arb(KS); fails=0
for (x1,x2,rx) in boxes:
    t1=twopi*flint.arb(x1.numerator)/x1.denominator; t2=twopi*flint.arb(x2.numerator)/x2.denominator; t3=-t1-t2
    rho=twopi*flint.arb(rx.numerator)/rx.denominator
    w=[t1.cos(),t1.sin(),t2.cos(),t2.sin(),t3.cos(),t3.sin()]
    d1=[-t1.sin(),t1.cos(),flint.arb(0),flint.arb(0),t3.sin(),-t3.cos()]
    d2=[flint.arb(0),flint.arb(0),-t2.sin(),t2.cos(),t3.sin(),-t3.cos()]
    ok=True
    for bi in range(4):
        Ma=[sum((Cb[bi][a][b]*w[b] for b in range(1,6)),Cb[bi][a][0]*w[0]) for a in range(6)]
        Ea=[Ab[bi][a]+Ma[a] for a in range(6)]; Fa=[Ab[bi][a]+Ma[a]*2 for a in range(6)]
        N=sum((Ea[a]*w[a] for a in range(6)),A0b[bi])
        D1=sum((Fa[a]*d1[a] for a in range(1,6)),Fa[0]*d1[0])
        D2=sum((Fa[a]*d2[a] for a in range(1,6)),Fa[0]*d2[0])
        tt=KSa-M2[bi]*rho*rho
        for e1 in (-1,1):
            for e2 in (-1,1):
                if not pd((N+D1*(e1*rho)+D2*(e2*rho)).tolist(),tt): ok=False
    if not ok: fails+=1
ks_ok=(fails==0)
print(f"(V) max over the torus of the blended variance < {KS} at r = {R}, s = {SB}: {ks_ok}  ({fails} failed boxes, {time.time()-t0:.1f} s)")

# ---- (D) the gap and the final inequality ----
gap2=flint.arb(283)/128-r/2
Lam=s*mh.sqrt()+(KSa/gap2).sqrt()
print(f"(D) gap2(r) = {gap2}")
print(f"    Lambda <= {Lam}   2 r Lambda = {2*r*Lam}   beta_W = {r/2}")
print(f"    rate -log(2 r Lambda) = {-(2*r*Lam).log()}")
ok=mt_ok and ks_ok and (2*r*Lam<1)
print("ALL OK" if ok else "FAILED")
