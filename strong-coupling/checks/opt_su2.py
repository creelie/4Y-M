"""Maximise the Rayleigh quotient  Hess S0(Q)(v,v)/|v|^2  jointly over SU(2) fields Q and
tangent vectors v on the torus (Z/L)^d (S0 = sum_p Re Tr Q_p). Gives LOWER bounds for
kappa(d) = sup_Q lambda_max(Hess S0) (also minimise for the bottom)."""
import torch, itertools, sys, math
torch.set_default_dtype(torch.float64)
def lattice(d,L):
    sites=list(itertools.product(range(L),repeat=d)); idx={}
    for x in sites:
        for mu in range(d): idx[(x,mu)]=len(idx)
    P=[]
    for x in sites:
        for mu in range(d):
            for nu in range(mu+1,d):
                xm=tuple((x[k]+(k==mu))%L for k in range(d)); xn=tuple((x[k]+(k==nu))%L for k in range(d))
                P.append((idx[(x,mu)],idx[(xm,nu)],idx[(xn,mu)],idx[(x,nu)]))
    return len(idx), torch.tensor(P)
s1=torch.tensor([[0,1],[1,0]],dtype=torch.complex128); s2=torch.tensor([[0,-1j],[1j,0]]); s3=torch.tensor([[1,0],[0,-1]],dtype=torch.complex128)
T=torch.stack([1j*s1,1j*s2,1j*s3])/math.sqrt(2)
def su2(q):  # q: (n,4) -> unit quaternion matrices
    q=q/q.norm(dim=1,keepdim=True)
    I=torch.eye(2,dtype=torch.complex128)
    return q[:,0,None,None]*I + math.sqrt(2)*torch.einsum('nk,kab->nab',q[:,1:].to(torch.complex128),T)
def alg(x): return torch.einsum('nk,kab->nab',x.to(torch.complex128),T)
def quad(q,x,P):
    Q=su2(q); X=alg(x)
    Qs=Q.conj().transpose(1,2)
    e1,e2,e3,e4=P[:,0],P[:,1],P[:,2],P[:,3]
    F=[Q[e1],Q[e2],Qs[e3],Qs[e4]]
    D1=[X[e1]@Q[e1], X[e2]@Q[e2], -Qs[e3]@X[e3], -Qs[e4]@X[e4]]
    D2=[X[e1]@X[e1]@Q[e1], X[e2]@X[e2]@Q[e2], Qs[e3]@X[e3]@X[e3], Qs[e4]@X[e4]@X[e4]]
    tot=0
    def pr(L):
        M=L[0]
        for A in L[1:]: M=M@A
        return M
    for i in range(4):
        L=list(F); L[i]=D2[i]; tot=tot+pr(L)
    for i in range(4):
        for j in range(i+1,4):
            L=list(F); L[i]=D1[i]; L[j]=D1[j]; tot=tot+2*pr(L)
    val=torch.diagonal(tot,dim1=1,dim2=2).sum(-1).real.sum()
    return val/(x*x).sum()
def run(d,L,restarts=6,steps=3000,sign=+1,seed=0):
    nE,P=lattice(d,L); best=-1e9; g=torch.Generator().manual_seed(seed)
    for r in range(restarts):
        q=torch.randn(nE,4,generator=g,requires_grad=True); x=torch.randn(nE,3,generator=g,requires_grad=True)
        opt=torch.optim.Adam([q,x],lr=0.03)
        for t in range(steps):
            opt.zero_grad(); f=-sign*quad(q,x,P); f.backward(); opt.step()
        v=sign*quad(q,x,P).item(); best=max(best,v)
        print(f"d={d} L={L} restart {r}: {'max' if sign>0 else '-min'} = {v:.6f}",flush=True)
    return best
if __name__=="__main__":
    d=int(sys.argv[1]); L=int(sys.argv[2]); R=int(sys.argv[3]) if len(sys.argv)>3 else 6
    print("best", run(d,L,R))
