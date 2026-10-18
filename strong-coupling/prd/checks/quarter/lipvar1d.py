import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
# functions of x2 = sin(al) on uniform S^3: weight cos^2 al, metric 2 sin(|da|/2)
n=300; al=np.linspace(-np.pi/2,np.pi/2,n); w=np.cos(al)**2; w/=w.sum()
D=2*np.sin(np.abs(al[:,None]-al[None])/2)
I,J=np.where(~np.eye(n,dtype=bool))
m=len(I);A=sp.csr_matrix((np.r_[np.ones(m),-np.ones(m)],(np.r_[np.arange(m),np.arange(m)],np.r_[I,J])),shape=(m,n))
best=0
rng=np.random.default_rng(0)
for s in range(8):
    g=np.sin(al) if s==0 else rng.normal(size=n)
    for it in range(15):
        a=w*(g-w@g); r=linprog(-a,A_ub=A,b_ub=D[I,J],bounds=[(-2,2)]*n,method='highs'); g=r.x
    v=w@(g-w@g)**2; best=max(best,v); print(s,v)
print("max var over functions of x2:",best," linear:",w@np.sin(al)**2)
