"""Haar Monte Carlo check of the rigorous moments of Section 5 at non-scalar diagonal b (not used by the proof).
Prints, for a few monomials u^alpha ubar^beta, the value of moments_rig with its error bound and the Monte
Carlo estimate of E_H[u^alpha ubar^beta e^{Re Tr(Ub)}] with its standard error."""
import sys, numpy as np
from series import Ring
from rseries import moments_rig
n=int(sys.argv[1]) if len(sys.argv)>1 else 400000
rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 1)
def haar(m):
    Z=(rng.standard_normal((m,3,3))+1j*rng.standard_normal((m,3,3)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z); d=np.diagonal(R,axis1=1,axis2=2); Q=Q*(d/np.abs(d))[:,None,:]
    return Q/(np.linalg.det(Q)**(1/3))[:,None,None]
R4=Ring(4); U=haar(n)
for b in [np.array([1.375,1.03125,1.03125],complex),np.array([1.03125,1.03125,0.515625+0.6875j]),np.array([0.5,0.25,-0.375+0.5j])]:
    w,ew=moments_rig(R4,b)
    wt=np.exp(np.einsum('nii,i->n',U,b).real)
    print('b =',b)
    for e in [{},{0:1},{0:1,9:1},{8:1,17:1},{1:1,10:1},{0:1,4:1,8:1},{0:2,9:2},{1:1,10:1,2:1,11:1},{3:1,12:1,8:1,17:1}]:
        m=[0]*18
        for k,v in e.items(): m[k]=v
        k=R4.idx[tuple(m)]
        f=np.ones(len(U),complex)
        for v,x in e.items():
            i,j=divmod(v%9,3); f*=(U[:,i,j] if v<9 else np.conj(U[:,i,j]))**x
        g=f*wt; mc=g.mean(); se=g.std()/np.sqrt(n)
        print(f'  {str(e):28s} exact {w[k].real:+.6f}{w[k].imag:+.6f}i (error <= {ew[k]:.1e})   MC {mc.real:+.6f}{mc.imag:+.6f}i (s.e. {se:.1e})')
