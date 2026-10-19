"""Floating-point value of 2 r sup_h of the bound of the corrector lemma on the face b = diag(r, rho, rho),
with the blended corrector s phi_h (best s) and with the degree-2 corrector (best C), at the exact moments.
Usage: python3 face2.py > face2.txt"""
import numpy as np, os, sys
from multiprocessing import Pool
def work(arg):
    r,q=arg
    import box2, rs6z
    w,ew=rs6z.moments_rig([r,q*r,q*r]); w=w/w[0].real
    C,kap,v2=box2.choose_C(w,np.array([r,q*r,q*r,0.0]),r)
    # blended corrector: q_h = s phi_h, i.e. C restricted to the span of phi_h
    gap=box2.GAPC-r/2
    W3=w[box2.P3]
    PH=box2.PHI2/2; Lphi=-(8/3)*PH+sum(b*G for b,G in zip([r,q*r,q*r,0.0],[np.array([g@p for p in box2.PHI2])/24 for g in box2.GT6]))
    M=sum((w6/6)*((D@PH.T).T@W3@(D@PH.T)).real for D,w6 in zip(box2.DX,box2.W6))
    best=9
    for s in np.linspace(0,0.6,121):
        X=PH+s*Lphi; mu=(X@w[:box2.N3]).real; V=(X@W3@X.T).real-np.outer(mu,mu)
        best=min(best,np.sqrt(max(np.linalg.eigvalsh(s*s*M).max(),0))+np.sqrt(max(np.linalg.eigvalsh(V).max(),0)/gap))
    return r,q,2*r*best,v2,kap
if __name__=='__main__':
    pts=[(r,q) for r in (1.44,1.50,1.54) for q in np.round(np.linspace(0.5,1.0,11),3)]
    with Pool(int(os.environ.get('NPROC','4'))) as p:
        for r,q,v1,v2,k in p.imap(work,pts):
            print(f'r={r:.2f} rho/r={q:.3f}: blended={v1:.4f} degree2={v2:.4f} kappa={k:.3f}',flush=True)
