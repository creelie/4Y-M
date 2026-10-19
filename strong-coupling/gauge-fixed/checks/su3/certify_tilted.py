"""Certificate for the tilted single-link bound of four-dimensional SU(3) in the forest gauge.

With r = 2 beta_W, for every b in rK (K the convex hull of SU(3)) and every h with |h|_HS = 1:

  (M)  nu_b |grad phi_h|^2 <= mt(r) := max_h E_H[ |grad phi_h|^2 e^{r Re Tr U} ],
  (V)  Var_{nu_b} <grad V, grad phi_h> <= r^2 kt(r),
       kt(r) := max_{W in SU(3), |h|=1} E_H[ (G_{W,h} - E_H G_{W,h})^2 e^{r Re Tr U} ],
       G_{W,h} = <grad phi_W, grad phi_h>,

and the Lipschitz rate is at most (3/8) sqrt(mt) + (3/8) r sqrt(kt/gap(r)) with gap(r) of Lemma D.
This script proves mt(0.98) <= MT and kt(0.98) <= KT and that 2 r Lambda < 1, in ball arithmetic.
W may be taken diagonal, W = diag(e^{i t1}, e^{i t2}, e^{-i(t1+t2)}); the maximum over the torus is
covered by dyadic boxes, each checked by a second-order bound (convexity in the linear term).

Usage: python3 certify_tilted.py [r] [MT] [KT]"""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from tilted_moments import *
from fractions import Fraction as Fr
import numpy as np, time
R=sys.argv[1] if len(sys.argv)>1 else '0.98'
MT=sys.argv[2] if len(sys.argv)>2 else '0.7123'
KT=sys.argv[3] if len(sys.argv)>3 else '0.1830'
r=flint.arb(R)
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

# ---- (0) self-test: at r = 0 the moment matrix is (4/9) I, the Haar value of Lemma A(4) ----
z=flint.arb(0)
for p in range(18):
    for q in range(18):
        v=ev(expect(G_poly(B18[p],B18[q])),z)[0]-(flint.arb(4)/9 if p==q else 0)
        assert v.contains(0) and v.rad()<1e-12, (p,q,v)
print("(0) at r = 0 the exact moment matrix is (4/9) I")

# ---- (M) the gradient moment ----
Q=[[None]*18 for _ in range(18)]
for p in range(18):
    for q in range(18):
        d=expect(G_poly(B18[p],B18[q]))
        re,im=ev(d,r)
        assert all(v[1]==0 for v in d.values()), 'imaginary part'
        Q[p][q]=re
mt_ok=pd(Q,flint.arb(MT))
print(f"(M) max_h E[|grad phi_h|^2 e^(r ReTrU)] < {MT} at r = {R}: {mt_ok}")
qm=np.array([[float(x.mid()) for x in row] for row in Q]); print('    float value',np.linalg.eigvalsh(qm).max())

# ---- (V) the tilted variance: exact tensor T[a,b,p,q] = E[g(d_a,e_p) g(d_b,e_q) e^(r ReTrU)] ----
g=[[g_poly(B6[a],B18[p]) for p in range(18)] for a in range(6)]
K=[[None]*6 for _ in range(6)]
for a in range(6):
    for b in range(6):
        Kab=[[flint.arb(0)]*18 for _ in range(18)]
        for p in range(18):
            for q in range(18):
                d1=expect(pmul(g[a][p],g[b][q])); d2=expect(pmul(g[b][p],g[a][q]))
                assert all(v[1]==0 for v in d1.values()) and all(v[1]==0 for v in d2.values())
                if blockof[p]!=blockof[q]:
                    assert not d1 and not d2, 'off-block entry'
                    continue
                Kab[p][q]=(ev(d1,r)[0]+ev(d2,r)[0])/2
        K[a][b]=Kab
print(f"    exact tensor built, block structure verified ({time.time()-t0:.1f} s)")
Kb=[[[flint.arb_mat([[K[a][b][p][q] for q in Bk] for p in Bk]) for b in range(6)] for a in range(6)] for Bk in BLOCKS]
def frob_upper(Mt):
    s=flint.arb(0)
    for row in Mt.tolist():
        for x in row: s+=x*x
    return s.sqrt()
M2=[8*sum((frob_upper(Kb[bi][a][b]) for a in range(6) for b in range(6)),flint.arb(0)) for bi in range(4)]
M2f=[float(m.upper()) for m in M2]
Kmid=[[[np.array([[float(x.mid()) for x in row] for row in Kb[bi][a][b].tolist()]) for b in range(6)] for a in range(6)] for bi in range(4)]
KTf=float(flint.arb(KT).mid())
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
        Ma=[sum(w[b]*Kmid[bi][a][b] for b in range(6)) for a in range(6)]
        N=sum(w[a]*Ma[a] for a in range(6)); D1=2*sum(d1[a]*Ma[a] for a in range(6)); D2=2*sum(d2[a]*Ma[a] for a in range(6))
        l=max(np.linalg.eigvalsh(N+s1*rho*D1+s2*rho*D2)[-1] for s1 in (-1,1) for s2 in (-1,1))
        ub=max(ub,l+M2f[bi]*rho**2)
    return ub
# float search for a box cover of the torus [0,1)^2 (theta = 2 pi x), dyadic boxes (center, half-width)
n0=32
todo=[(Fr(2*i+1,2*n0),Fr(2*j+1,2*n0),Fr(1,2*n0)) for i in range(n0) for j in range(n0)]
boxes=[]
while todo:
    x1,x2,rx=todo.pop()
    if ub_float(x1,x2,rx)<KTf*(1-1e-6): boxes.append((x1,x2,rx)); continue
    assert rx>Fr(1,2**20), 'subdivision too deep'
    h=rx/2
    todo+=[(x1+s1*h,x2+s2*h,h) for s1 in (-1,1) for s2 in (-1,1)]
area=sum(4*rx*rx for _,_,rx in boxes)
assert area==1, area
print(f"    {len(boxes)} dyadic boxes cover the torus exactly (area {area}); float search {time.time()-t0:.1f} s")
# rigorous check of every box
twopi=2*flint.arb.pi(); KTa=flint.arb(KT); fails=0
for (x1,x2,rx) in boxes:
    t1=twopi*flint.arb(x1.numerator)/x1.denominator; t2=twopi*flint.arb(x2.numerator)/x2.denominator; t3=-t1-t2
    rho=twopi*flint.arb(rx.numerator)/rx.denominator
    w=[t1.cos(),t1.sin(),t2.cos(),t2.sin(),t3.cos(),t3.sin()]
    d1=[-t1.sin(),t1.cos(),flint.arb(0),flint.arb(0),t3.sin(),-t3.cos()]
    d2=[flint.arb(0),flint.arb(0),-t2.sin(),t2.cos(),t3.sin(),-t3.cos()]
    ok=True
    for bi in range(4):
        Ma=[sum((Kb[bi][a][b]*w[b] for b in range(1,6)),Kb[bi][a][0]*w[0]) for a in range(6)]
        N=sum((Ma[a]*w[a] for a in range(1,6)),Ma[0]*w[0])
        D1=sum((Ma[a]*(2*d1[a]) for a in range(1,6)),Ma[0]*(2*d1[0]))
        D2=sum((Ma[a]*(2*d2[a]) for a in range(1,6)),Ma[0]*(2*d2[0]))
        tt=KTa-M2[bi]*rho*rho
        for s1 in (-1,1):
            for s2 in (-1,1):
                A=(N+D1*(s1*rho)+D2*(s2*rho)).tolist()
                if not pd(A,tt): ok=False
    if not ok: fails+=1
kt_ok=(fails==0)
print(f"(V) max over the torus of the tilted variance < {KT} at r = {R}: {kt_ok}  ({fails} failed boxes, {time.time()-t0:.1f} s)")

# ---- final inequality ----
g2=flint.arb('2.06548'); mu0=flint.arb(7)/6
gap=(mu0.sqrt()-g2.sqrt()*r/2)**2+flint.arb(3)/2-r
assert gap>0 and g2.sqrt()*r<2*mu0.sqrt()
Lam=flint.arb(3)/8*flint.arb(MT).sqrt()+flint.arb(3)/8*r*(flint.arb(KT)/gap).sqrt()
print(f"    gap(r) = {gap}   Lambda <= {Lam}   2 r Lambda = {2*r*Lam}   beta_W = {r/2}")
print(f"    rate -log(2 r Lambda) = {-(2*r*Lam).log()}")
ok=mt_ok and kt_ok and (2*r*Lam<1)
print("ALL OK" if ok else "FAILED")
