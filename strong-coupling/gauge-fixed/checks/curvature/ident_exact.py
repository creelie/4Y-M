"""Exact compressed curl-curl operator at the identity field on tori of side 2 and 4 (forest gauge)."""
import sympy as sp, itertools, numpy as np
SIG=[0,0,2,3,3,1,2,1,1,2,1,3,3,2,0,0]; d=4
def op(L):
    sites=list(itertools.product(range(L),repeat=d)); idx={s:i for i,s in enumerate(sites)}
    lk=lambda x,mu: idx[tuple(x)]*d+mu
    sh=lambda x,mu: [ (x[k]+(k==mu))%L for k in range(d)]
    n=len(sites)*d; A=np.zeros((n,n),dtype=np.int64)
    for x in sites:
        for mu in range(d):
            for nu in range(mu+1,d):
                row=np.zeros(n,dtype=np.int64)
                for (l,v) in [(lk(x,mu),1),(lk(sh(x,mu),nu),1),(lk(sh(x,nu),mu),-1),(lk(x,nu),-1)]: row[l]+=v
                A+=np.outer(row,row)
    act=[lk(x,m) for x in sites for m in range(d) if m!=SIG[(x[0]%2)+2*(x[1]%2)+4*(x[2]%2)+8*(x[3]%2)]]
    return A[np.ix_(act,act)]
lam=sp.symbols('lam')
for L in (2,):
    Ac=sp.Matrix(op(L))
    P=Ac.charpoly(lam).as_expr()
    print('L=2 char poly factors:',sp.factor_list(P))
Ac4=op(4); w=np.linalg.eigvalsh(Ac4.astype(float)); print('L=4 top',w[-1],' 9+sqrt17 =',9+17**0.5)
