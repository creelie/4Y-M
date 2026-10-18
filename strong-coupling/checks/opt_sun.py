"""Joint maximisation of Hess S0(Q)(v,v)/|v|^2 over SU(N) fields on (Z/L)^d (lower bounds for kappa_N(d)).
Q_e = matrix_exp(A_e), A_e in su(N); v_e = X_e in su(N); metric -Tr(XY)."""
import torch, itertools, sys, math
torch.set_default_dtype(torch.float64)
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from opt_su2 import lattice
def su(N,raw):   # raw (n,N,N) complex -> anti-Hermitian traceless
    A=raw-raw.conj().transpose(-1,-2)
    tr=torch.diagonal(A,dim1=-2,dim2=-1).sum(-1)/N
    return A-tr[:,None,None]*torch.eye(N,dtype=A.dtype)
def quad(N,ra,rx,P):
    Q=torch.matrix_exp(su(N,ra)); X=su(N,rx); Qs=Q.conj().transpose(1,2)
    e1,e2,e3,e4=P[:,0],P[:,1],P[:,2],P[:,3]
    F=[Q[e1],Q[e2],Qs[e3],Qs[e4]]
    D1=[X[e1]@Q[e1], X[e2]@Q[e2], -Qs[e3]@X[e3], -Qs[e4]@X[e4]]
    D2=[X[e1]@X[e1]@Q[e1], X[e2]@X[e2]@Q[e2], Qs[e3]@X[e3]@X[e3], Qs[e4]@X[e4]@X[e4]]
    def pr(L):
        M=L[0]
        for A in L[1:]: M=M@A
        return M
    tot=0
    for i in range(4):
        L=list(F); L[i]=D2[i]; tot=tot+pr(L)
    for i in range(4):
        for j in range(i+1,4):
            L=list(F); L[i]=D1[i]; L[j]=D1[j]; tot=tot+2*pr(L)
    val=torch.diagonal(tot,dim1=1,dim2=2).sum(-1).real.sum()
    nrm=-torch.diagonal(X@X,dim1=1,dim2=2).sum(-1).real.sum()
    return val/nrm
def run(N,d,L,restarts,steps=3000,sign=1,seed=0):
    nE,P=lattice(d,L); g=torch.Generator().manual_seed(seed); best=-1e9
    for r in range(restarts):
        ra=(torch.randn(nE,N,N,generator=g)+1j*torch.randn(nE,N,N,generator=g)).requires_grad_(True)
        rx=(torch.randn(nE,N,N,generator=g)+1j*torch.randn(nE,N,N,generator=g)).requires_grad_(True)
        opt=torch.optim.Adam([ra,rx],lr=0.03)
        for t in range(steps):
            opt.zero_grad(); f=-sign*quad(N,ra,rx,P); f.backward(); opt.step()
        v=sign*quad(N,ra,rx,P).item(); best=max(best,v); print(f"N={N} d={d} L={L} r={r} sign={sign}: {v:.6f}",flush=True)
    return best
if __name__=="__main__":
    N,d,L,R=map(int,sys.argv[1:5]); sign=int(sys.argv[5]) if len(sys.argv)>5 else 1
    print("best",run(N,d,L,R,sign=sign))
