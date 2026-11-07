"""Flat constrained gap in d dims via Bloch decomposition over the coarse lattice (abelian; flat SU(2) = 3 copies).
Cell: b^d sites, d*b^d links. Constraint: corner map (sum of b links on the straight line from block corner, each direction).
Quadratic form: Wilson at flat = sum_p |(da)_p|^2. Gauge removed exactly as in gaugegap. Output: min gap*b^2 over a k-grid.
Heuristic float numerics; k-grid sampling, no certification."""
import sys, itertools, json, numpy as np, scipy.linalg as sl
d=int(sys.argv[1]); b=int(sys.argv[2]); K=int(sys.argv[3]); mode=sys.argv[4] if len(sys.argv)>4 else 'corner'
cells=list(itertools.product(range(b),repeat=d)); idx={c:i for i,c in enumerate(cells)}; ns=len(cells)
def site(x,k):  # returns (index, phase) for fine site x possibly outside cell
    q=[xi//b for xi in x]; r=tuple(xi%b for xi in x); return idx[r], np.exp(1j*np.dot(k,q))
def link(x,mu,k): i,p=site(x,k); return mu*ns+i,p
nl=d*ns
def mats(k):
    # d operator: sites->links
    D=np.zeros((nl,ns),complex)
    for x in cells:
        for mu in range(d):
            y=list(x); y[mu]+=1; l,_=link(x,mu,k); i,pi=site(x,k); j,pj=site(y,k)
            D[l,j]+=pj; D[l,i]-=pi
    rows=[]
    for x in cells:
        for mu in range(d):
            for nu in range(mu+1,d):
                r=np.zeros(nl,complex); xm=list(x); xm[mu]+=1; xn=list(x); xn[nu]+=1
                for (y,m,s) in ((x,mu,1),(xm,nu,1),(xn,mu,-1),(x,nu,-1)):
                    l,p=link(y,m,k); r[l]+=s*p
                rows.append(r)
    C=np.array(rows); H=C.conj().T@C
    J=[]
    for mu in range(d):
        offs=[()] if mode=='corner' else list(itertools.product(range(b),repeat=d-1))
        for o in offs:
            r=np.zeros(nl,complex)
            for s in range(b):
                x=[0]*d; oo=list(o)
                for nu in range(d):
                    if nu!=mu: x[nu]=oo.pop(0) if oo else 0
                x[mu]=s; l,p=link(x,mu,k); r[l]+=p
            J.append(r)
    return H,np.array(J),D
best=(9e9,None)
for kk in itertools.product(range(K),repeat=d):
    k=2*np.pi*np.array(kk)/K
    H,J,D=mats(k); N=sl.null_space(J); M=J@D
    u,s,vh=np.linalg.svd(M); C=vh[(s>1e-9*max(s.max(),1)).sum():].conj().T
    G=N.conj().T@D@C
    if G.size:
        ug,sg,_=np.linalg.svd(G); kg=int((sg>1e-8).sum()); Q=ug[:,kg:]
    else: Q=np.eye(N.shape[1])
    e=np.linalg.eigvalsh(Q.conj().T@N.conj().T@H@N@Q)
    if e[0]<best[0]: best=(float(e[0]),kk)
print(json.dumps({'d':d,'b':b,'K':K,'mode':mode,'gap':best[0],'gap_b2':best[0]*b*b,'argmin_k':best[1]}),flush=True)
