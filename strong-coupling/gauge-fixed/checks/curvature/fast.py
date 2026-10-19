import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from thess import *
import scipy.sparse as sp, scipy.sparse.linalg as sla
def kblocks(lat,Q):
    P=lat.P; Q1=Q[P[:,0]];Q2=Q[P[:,1]];Q3=Q[P[:,2]];Q4=Q[P[:,3]]
    W=Q1@Q2@H(Q3)@H(Q4); npl=len(P)
    Tb=T[None].expand(npl,3,2,2)
    TT=[Tb, Q1[:,None]@Tb@H(Q1)[:,None], Q4[:,None]@Tb@H(Q4)[:,None], Tb]
    Wb=W[:,None,None]; z=torch.zeros(npl,3,3)
    blocks=[[z for _ in range(4)] for _ in range(4)]
    for i in range(4):
        A=rt(TT[i][:,:,None]@TT[i][:,None,:]@Wb); blocks[i][i]=0.5*(A+A.transpose(1,2))
    for (i,j) in [(0,1),(2,3)]:
        A=rt(TT[i][:,:,None]@TT[j][:,None,:]@Wb); blocks[i][j]=blocks[i][j]+A; blocks[j][i]=blocks[j][i]+A.transpose(1,2)
    for i in (0,1):
        for j in (2,3):
            A=-rt(TT[i][:,:,None]@Wb@TT[j][:,None,:]); blocks[i][j]=blocks[i][j]+A; blocks[j][i]=blocks[j][i]+A.transpose(1,2)
    return torch.stack([torch.cat([blocks[i][j] for j in range(4)],2) for i in range(4)],1).reshape(npl,12,12)
class Fast:
    def __init__(s,lat):
        s.lat=lat; n3=lat.nl*3
        mask=np.ones(lat.nl,bool); mask[list(lat.F)]=False
        s.mask3=np.repeat(mask,3)
        s.rows=lat.rows.numpy(); s.cols=lat.cols.numpy(); s.n3=n3
        P=lat.P; s.r12=(P[:,:,None]*3+torch.arange(3)[None,None]).reshape(len(P),12)
    def top(s,Q,sign,v0=None):
        K=kblocks(s.lat,Q).detach().numpy().reshape(-1)
        Hs=sp.coo_matrix((K,(s.rows,s.cols)),shape=(s.n3,s.n3)).tocsr()
        Ha=Hs[s.mask3][:,s.mask3]
        w,V=sla.eigsh(sign*Ha,k=1,which='LA',v0=v0,tol=1e-10)
        v=np.zeros(s.n3); v[s.mask3]=V[:,0]
        return w[0],v,V[:,0]
    def grad_obj(s,Q,v):
        K=kblocks(s.lat,Q); vt=torch.tensor(v)[s.r12]
        return torch.einsum('pi,pij,pj->',vt,K,vt)
def run(L,seed,nstarts,steps,lr,init=None):
    torch.manual_seed(seed); np.random.seed(seed)
    lat=Lat(L,'forest'); f=Fast(lat); best_all=0
    for k in range(nstarts):
        sc=[0.3,1.0,2.0,3.0,0.6,1.5][k%6]; sign=1 if k%2==0 else -1
        Q0=expm(torch.randn(lat.nl,3)*sc)
        th=torch.zeros(lat.nl,3,requires_grad=True); opt=torch.optim.Adam([th],lr=lr); best=-1; v0=None
        for it in range(steps):
            Q=expm(th)@Q0
            val,v,v0=f.top(Q.detach(),sign,v0)
            best=max(best,val)
            g=f.grad_obj(Q,v)*sign
            opt.zero_grad(); (-g).backward(); opt.step()
        best_all=max(best_all,best)
        print(L,k,sign,sc,best,flush=True)
        torch.save((expm(th)@Q0).detach(),f'fast_{L}_{seed}_{k}.pt')
    print('best',best_all,flush=True)
if __name__=="__main__":
    run(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),float(sys.argv[5]))
