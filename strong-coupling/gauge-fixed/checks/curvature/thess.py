import torch, itertools, numpy as np, sys
torch.set_default_dtype(torch.float64)
SIG=[0,0,2,3,3,1,2,1,1,2,1,3,3,2,0,0]
d=4
s1=torch.tensor([[0,1],[1,0]],dtype=torch.complex128); s2=torch.tensor([[0,-1j],[1j,0]],dtype=torch.complex128); s3=torch.tensor([[1,0],[0,-1]],dtype=torch.complex128)
T=torch.stack([1j*s1,1j*s2,1j*s3])/np.sqrt(2)
def H(A): return A.conj().transpose(-1,-2)
def expm(v):
    X=torch.einsum('...a,aij->...ij',v.to(torch.complex128),T)
    return torch.linalg.matrix_exp(X)
def rt(A): return torch.diagonal(A,dim1=-2,dim2=-1).sum(-1).real
class Lat:
    def __init__(s,L,gauge):
        s.L=L; sites=list(itertools.product(range(L),repeat=d)); idx={x:i for i,x in enumerate(sites)}
        lk=lambda x,mu: idx[tuple(x)]*d+mu
        s.nl=len(sites)*d; P=[]
        for x in sites:
            for mu in range(d):
                for nu in range(mu+1,d):
                    xm=list(x);xm[mu]=(xm[mu]+1)%L; xn=list(x);xn[nu]=(xn[nu]+1)%L
                    P.append((lk(x,mu),lk(xm,nu),lk(xn,mu),lk(x,nu)))
        s.P=torch.tensor(P)
        if gauge=='forest': F={lk(x,SIG[(x[0]%2)+2*(x[1]%2)+4*(x[2]%2)+8*(x[3]%2)]) for x in sites}
        elif gauge=='temporal': F={lk(x,0) for x in sites}
        else: F=set()
        s.F=F
        s.act=torch.tensor([l for l in range(s.nl) if l not in F])
        s.ii=(s.act[:,None]*3+torch.arange(3)[None]).reshape(-1)
        P=s.P; npl=len(P)
        r=(P[:,:,None]*3+torch.arange(3)[None,None]).reshape(npl,12)
        s.rows=r[:,:,None].expand(npl,12,12).reshape(-1); s.cols=r[:,None,:].expand(npl,12,12).reshape(-1)
    def hess(s,Q,full=False):
        P=s.P; Q1=Q[P[:,0]];Q2=Q[P[:,1]];Q3=Q[P[:,2]];Q4=Q[P[:,3]]
        W=Q1@Q2@H(Q3)@H(Q4); npl=len(P)
        Tb=T[None].expand(npl,3,2,2)
        TT=[Tb, Q1[:,None]@Tb@H(Q1)[:,None], Q4[:,None]@Tb@H(Q4)[:,None], Tb]
        Wb=W[:,None,None]
        z=torch.zeros(npl,3,3)
        blocks=[[z for _ in range(4)] for _ in range(4)]
        for i in range(4):
            A=rt(TT[i][:,:,None]@TT[i][:,None,:]@Wb); blocks[i][i]=0.5*(A+A.transpose(1,2))
        for (i,j) in [(0,1),(2,3)]:
            A=rt(TT[i][:,:,None]@TT[j][:,None,:]@Wb); blocks[i][j]=blocks[i][j]+A; blocks[j][i]=blocks[j][i]+A.transpose(1,2)
        for i in (0,1):
            for j in (2,3):
                A=-rt(TT[i][:,:,None]@Wb@TT[j][:,None,:]); blocks[i][j]=blocks[i][j]+A; blocks[j][i]=blocks[j][i]+A.transpose(1,2)
        K=torch.stack([torch.cat([blocks[i][j] for j in range(4)],2) for i in range(4)],1).reshape(npl,12,12)
        n3=s.nl*3
        Hm=torch.zeros(n3*n3).index_add_(0,s.rows*n3+s.cols,K.reshape(-1)).reshape(n3,n3)
        return Hm if full else Hm[s.ii][:,s.ii]
def optimise(lat,sign,Q0,steps=400,lr=0.05):
    th=torch.zeros(lat.nl,3,requires_grad=True)
    opt=torch.optim.Adam([th],lr=lr); best=-1e9
    for it in range(steps):
        Q=expm(th)@Q0
        ev=torch.linalg.eigvalsh(lat.hess(Q))
        val=ev[-1] if sign>0 else -ev[0]
        best=max(best,val.item())
        opt.zero_grad(); (-val).backward(); opt.step()
    return best, (expm(th)@Q0).detach()
if __name__=="__main__":
    g=sys.argv[1]; L=int(sys.argv[2]); n=int(sys.argv[3]); seed=int(sys.argv[4]); steps=int(sys.argv[5])
    torch.manual_seed(seed); lat=Lat(L,g)
    I=torch.eye(2,dtype=torch.complex128).expand(lat.nl,2,2)
    ev=torch.linalg.eigvalsh(lat.hess(I)); print('identity',ev[0].item(),ev[-1].item(),flush=True)
    allbest=0
    for k in range(n):
        sc=float(torch.rand(1))*3
        Q0=expm(torch.randn(lat.nl,3)*sc)
        sign=1 if k%2==0 else -1
        b,Q=optimise(lat,sign,Q0,steps)
        allbest=max(allbest,b); print(k,sign,round(sc,2),b,flush=True)
        torch.save(Q,f'best_{g}_{L}_{seed}_{k}.pt')
    print('best',allbest)
