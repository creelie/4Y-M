"""Route R1 for the SU(3) single-link rate: psi = f/lambda_F, residual R = <grad V, grad f>/lambda_F - nu f.
L(b) <= sqrt(lmax A)/lambda_F + sqrt(lmax B / Lambda(r)), A = nu|grad f_h|^2, B = Var R_h (quadratic forms in h)."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from lip import haar, T, coords
import numpy as np
lamF=8/3; mu0=7/6
def Lam(r):
    G=np.sqrt(8/3)*r
    return max(np.sqrt(mu0)-G/2,0)**2+1.5-r
def basis():
    H=[]
    for k in range(18):
        h=np.zeros(18); h[k]=1
        H.append((h[:9]+1j*h[9:]).reshape(3,3))
    return H
# f_h = ReTr(U h); with h = E_ij (real) or iE_ij: grad_a f_h = ReTr(T_a U h)
def run(b,n=400000,seed=1):
    rng=np.random.default_rng(seed)
    import lip; lip.rng=rng
    U=haar(n)
    w=np.exp(np.einsum('nij,ji->n',U,b).real); w/=w.sum()
    TU=np.einsum('aij,njk->naik',T,U)            # n x 8 x 3 x 3
    gV=np.einsum('naij,ji->na',TU,b).real         # n x 8
    H=basis()
    gF=np.stack([np.einsum('naij,ji->na',TU,h).real for h in H],2)   # n x 8 x 18
    F=np.stack([np.einsum('nij,ji->n',U,h).real for h in H],1)       # n x 18
    A=np.einsum('n,nap,naq->pq',w,gF,gF)
    R=np.einsum('na,nap->np',gV,gF)/lamF
    mR=w@R; Rc=R-mR
    B=np.einsum('n,np,nq->pq',w,Rc,Rc)
    return np.linalg.eigvalsh(A).max(), np.linalg.eigvalsh(B).max()
for r in (0.0,0.25,0.5,0.75,1.0,1.25):
    for name,b in (('rI',r*np.eye(3)),('r(1,1,-1)',r*np.diag([1,1,-1.0]))):
        a,bb=run(b)
        Lr=Lam(r); ell=np.sqrt(a)/lamF+np.sqrt(bb/Lr)
        print(f"r={r} {name}: lmaxA={a:.4f} lmaxB={bb:.5f} Lambda={Lr:.4f} ell={ell:.4f} 2r*ell={2*r*ell:.3f}",flush=True)
