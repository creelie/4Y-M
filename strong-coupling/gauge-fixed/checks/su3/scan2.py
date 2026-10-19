import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from lip import *
rng2=np.random.default_rng(5)
U0=haar(100000)
for r in (0.25,0.5,0.75,1.0,1.25,1.5):
    best=(0,None)
    for s in [(r,r,r),(r,r,0),(r,0,0),(r,r,r/2)]:
        for th in (0,np.pi/3,2*np.pi/3,np.pi):
            b=np.diag([s[0],s[1],s[2]*np.exp(1j*th)])
            w=np.exp(np.einsum('nij,ji->n',U0,b).real); ess=w.sum()**2/(w**2).sum()
            L,V=Lrate(b,D=2,n=100000,U0=U0)
            if L>best[0]: best=(L,(s,th))
            print(f"r={r} s={s} th={th:.3f} L={L:.4f} Var={V:.4f} 2rL={2*r*L:.3f} ess={ess:.0f}",flush=True)
    print('BEST',r,best,flush=True)
