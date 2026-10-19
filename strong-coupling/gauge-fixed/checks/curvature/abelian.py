"""Is the extremal gauge-fixed Hessian reached by an abelian (U(1) inside SU(2)) field? Period-2 angles."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fast import *
def run(L,seed,nst,steps,lr):
    torch.manual_seed(seed); lat=Lat(L,'forest'); f=Fast(lat); best_all=0
    sites=list(itertools.product(range(L),repeat=d))
    # period-2 angle per (residue, direction): 64 parameters, extended periodically
    res=[(x[0]%2)+2*(x[1]%2)+4*(x[2]%2)+8*(x[3]%2) for x in sites]
    rix=torch.tensor([res[i]*4+mu for i in range(len(sites)) for mu in range(4)])
    for k in range(nst):
        sign=1 if k%2==0 else -1
        th=(torch.rand(64)*2*np.pi).requires_grad_(True); opt=torch.optim.Adam([th],lr=lr); best=-1; v0=None
        for it in range(steps):
            v3=torch.zeros(lat.nl,3); v3=v3.clone(); 
            ang=th[rix]
            Q=expm(torch.stack([torch.zeros_like(ang),torch.zeros_like(ang),ang],1))
            val,v,v0=f.top(Q.detach(),sign,v0); best=max(best,val)
            g=f.grad_obj(Q,v)*sign; opt.zero_grad(); (-g).backward(); opt.step()
        print(k,sign,best,flush=True); best_all=max(best_all,best)
        torch.save(th.detach(),f'abel_{seed}_{k}.pt')
    print('best abelian period-2',best_all)
run(4,5,8,500,0.05)
