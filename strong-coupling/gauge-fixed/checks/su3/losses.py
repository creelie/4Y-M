"""The three losses of the tilted bound at r = 0.98 (su3.md, section 5); not used by the proof.

1. Z(rI) and the normalised gradient moment max_h nu_{rI}|grad phi_h|^2 (exact moments);
2. the decoupled bound of Lemma H at W = I and the true variance max_h Var_{nu_{rI}} G_{I,h} (exact);
3. Galerkin upper bounds (polynomials of degree two, 10^5 Haar points) for the gap of nu_b."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from tilted_moments import *
from lip import haar, features
import numpy as np
r=0.98; ra=flint.arb(str(r)); Z=float(CI((),ra).mid()); B18=basis18(); Id={(k,k):ONE for k in range(3)}
def val(p): return float(ev(expect(p),ra)[0].mid())
Q=np.array([[val(G_poly(B18[p],B18[q])) for q in range(18)] for p in range(18)])
print(f"1. Z(rI) = {Z:.4f}; max_h E[|grad phi_h|^2 e] = {np.linalg.eigvalsh(Q).max():.4f}, normalised {np.linalg.eigvalsh(Q).max()/Z:.4f}")
G=[G_poly(Id,B18[p]) for p in range(18)]; g=[g_poly(Id,B18[p]) for p in range(18)]
mean=np.array([val(G[p]) for p in range(18)])/Z
S2=np.array([[val(pmul(G[p],G[q])) for q in range(18)] for p in range(18)])/Z
T=np.array([[val(pmul(g[p],g[q])) for q in range(18)] for p in range(18)])
print(f"2. bound of Lemma H at W = I: {np.linalg.eigvalsh(T).max():.4f}; true max_h Var G_(I,h) under nu_rI: {np.linalg.eigvalsh(S2-np.outer(mean,mean)).max():.4f}")
U0=haar(100000)
def galerkin_gap(b,D=2):
    X,F,Gr=features(U0,D)
    w=np.exp(np.einsum('nij,ji->n',U0,b).real); w/=w.sum()
    mF=w@F; C=(F-mF).T@((F-mF)*w[:,None]); Gw=(Gr*np.sqrt(w)[None,:,None]).reshape(-1,Gr.shape[2]); A=Gw.T@Gw
    e,V=np.linalg.eigh(C); keep=e>1e-8*e.max(); P=V[:,keep]/np.sqrt(e[keep])
    return np.linalg.eigvalsh(P.T@A@P).min()
gl=(np.sqrt(7/6)-np.sqrt(2.06548)*r/2)**2+1.5-r
for name,b in [('b = 0',0*np.eye(3)),('b = rI',r*np.eye(3)),('b = r diag(1,1,0)',np.diag([r,r,0]))]:
    print(f"3. {name}: Galerkin gap {galerkin_gap(b):.4f}   (Lemma D: {gl:.4f}; exact at b = 0: 8/3)")
