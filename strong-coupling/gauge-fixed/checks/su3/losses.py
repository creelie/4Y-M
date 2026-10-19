"""The losses of the blended bound at r = 1.14, s = 0.322 (su3.md, Section 6); not used by the proof.

1. Z(rI), the bounds mt and mh of Lemmas G and G', and the normalised moment max_h nu_{rI}|grad phi_h|^2;
2. the decoupled quadratic form of Lemma H' at W = I and at W = omega I, against the true
   max_h Var_{nu_{rI}}((1 - 8s/3) phi_h + s Gamma) (all exact moments);
3. the gap: see gap_lemma.py (Galerkin values against Lemma D')."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from tilted_moments import *
import numpy as np
r=1.14; s=0.322; ra=flint.arb(str(r)); Z=float(CI((),ra).mid()); B18=basis18()
c1=1-8*s/3; c2=s*r
def val(p): return float(ev(expect(p),ra)[0].mid())
Q=np.array([[val(G_poly(B18[p],B18[q])) for q in range(18)] for p in range(18)])
mt=np.linalg.eigvalsh(Q).max(); mh=mt/(mt+5/9*np.exp(-r/5))
print(f"1. Z(rI) = {Z:.4f}; mt = {mt:.4f}, mh = {mh:.4f}; true max_h nu_rI|grad phi_h|^2 = {mt/Z:.4f}")
phi=[scal(tr1(B18[p],1)+tr1(dag(B18[p]),-1),Fr(1,2)) for p in range(18)]
def X(W,p): return [(a,b,gscale(v,Fr(c1).limit_denominator(10**9))) for a,b,v in phi[p]]+scal(g_poly(W,B18[p]),Fr(c2).limit_denominator(10**9))
Id={(k,k):ONE for k in range(3)}
w3=Fr(866025403784439,10**15)   # sqrt(3)/2, rational approximation
Om={(k,k):gq(Fr(-1,2),w3) for k in range(3)}
for name,W in [('W = I',Id),('W = omega I (rational approximation of omega)',Om)]:
    T=np.array([[val(pmul(X(W,p),X(W,q))) for q in range(18)] for p in range(18)])
    print(f"2. decoupled form of Lemma H' at {name}: {np.linalg.eigvalsh(T).max():.5f}")
Xs=[X(Id,p) for p in range(18)]
mean=np.array([val(Xs[p]) for p in range(18)])/Z
S2=np.array([[val(pmul(Xs[p],Xs[q])) for q in range(18)] for p in range(18)])/Z
print(f"   true max_h Var_(nu_rI)((1 - 8s/3) phi_h + s Gamma) = {np.linalg.eigvalsh(S2-np.outer(mean,mean)).max():.5f}")
