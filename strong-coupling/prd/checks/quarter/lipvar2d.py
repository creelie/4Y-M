import numpy as np, scipy.sparse as sp, sys, time
from scipy.optimize import linprog
from scipy.special import iv
n=int(sys.argv[1]); 
# Fibonacci points on the full sphere, keep r>=0 hemisphere (coords x1,x2,r)
N=2*n; i=np.arange(N)+0.5; z=1-2*i/N; ph=np.pi*(1+5**0.5)*i
P=np.c_[np.sqrt(1-z**2)*np.cos(ph),np.sqrt(1-z**2)*np.sin(ph),z]; P=P[P[:,2]>0]
D=np.linalg.norm(P[:,None]-P[None],axis=2); n=len(P)
I,J=np.where(~np.eye(n,dtype=bool)); m=len(I)
A=sp.csr_matrix((np.r_[np.ones(m),-np.ones(m)],(np.r_[np.arange(m),np.arange(m)],np.r_[I,J])),shape=(m,n))
rng=np.random.default_rng(0)
for k in [float(x) for x in sys.argv[2:]]:
    w=P[:,2]*np.exp(k*P[:,0]); w/=w.sum()
    best=0
    for s in range(3):
        g=[P[:,1],P[:,0],rng.normal(size=n)][s]
        for it in range(8):
            a=w*(g-w@g); r=linprog(-a,A_ub=A,b_ub=D[I,J],bounds=[(-2,2)]*n,method='highs'); g=r.x
        v=w@(g-w@g)**2; best=max(best,v)
    lin=max(w@(P[:,1])**2, w@(P[:,0]-w@P[:,0])**2)
    print(f"n={n} kappa={k}: max Var={best:.5f}, linear max={lin:.5f}, A/kappa={iv(2,k)/iv(1,k)/k if k else .25:.5f}",flush=True)
