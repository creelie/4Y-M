import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, itertools, sys
SIG=[0,0,2,3,3,1,2,1,1,2,1,3,3,2,0,0]
d=4
def run(L):
    sites=list(itertools.product(range(L),repeat=d))
    idx={s:i for i,s in enumerate(sites)}
    def lk(x,mu): return idx[tuple(x)]*d+mu
    def sh(x,mu): y=list(x); y[mu]=(y[mu]+1)%L; return y
    rows=[];cols=[];vals=[];r=0
    for x in sites:
        for mu in range(d):
            for nu in range(mu+1,d):
                for (l,v) in [(lk(x,mu),1),(lk(sh(x,mu),nu),1),(lk(sh(x,nu),mu),-1),(lk(x,nu),-1)]:
                    rows.append(r);cols.append(l);vals.append(v)
                r+=1
    B=sp.csr_matrix((vals,(rows,cols)),shape=(r,len(sites)*d))
    A=(B.T@B).tocsr()
    F=set()
    for x in sites:
        s=SIG[(x[0]%2)+2*(x[1]%2)+4*(x[2]%2)+8*(x[3]%2)]
        F.add(lk(x,s))
    act=[l for l in range(len(sites)*d) if l not in F]
    Ac=A[act][:,act]
    full=sla.eigsh(A,k=1,which='LA')[0][0]
    comp=sla.eigsh(Ac,k=1,which='LA')[0][0]
    # temporal gauge for comparison: fix direction 0 everywhere
    actT=[l for l in range(len(sites)*d) if l%d!=0]
    compT=sla.eigsh(A[actT][:,actT],k=1,which='LA')[0][0]
    print(L,"full",full,"forest",comp,"temporal",compT)
for L in map(int,sys.argv[1:]): run(L)
