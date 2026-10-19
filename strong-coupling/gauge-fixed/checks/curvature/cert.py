"""Local certificate for the forest-gauge Hessian bound.

Find coefficients alpha (per plaquette, per active position), lambda (per plaquette, per
non-F corner) and gamma (per vertex, per active end) such that
  S_p = sum_i alpha_i |Y_i|^2 - f_p - 2 sum_k lambda_k corner_k  >= 0   for every W,
  M_y (vertex matrix, diagonal gamma, perpendicular entries lambda, collinear 0) >= 0,
and every active link has  sum_{p ni e} alpha + gamma_out + gamma_in <= h.
Then +-Hess S_0 <= h |X|^2 on {X_F = 0}, for every field and every U(N).
"""
import numpy as np, itertools, sys, cvxpy as cp
SIG=[0,0,2,3,3,1,2,1,1,2,1,3,3,2,0,0]
def sig(r): return SIG[r[0]+2*r[1]+4*r[2]+8*r[3]]
def add(r,mu,s=1): t=list(r); t[mu]=(t[mu]+s)%2; return tuple(t)
res=list(itertools.product((0,1),repeat=4))
s1=np.array([[0,1],[1,0]],complex); s2=np.array([[0,-1j],[1j,0]]); s3=np.array([[1,0],[0,-1]],complex)
T=np.array([1j*s1,1j*s2,1j*s3])/np.sqrt(2)
rt=lambda A: np.real(np.trace(A))
ip=lambda A,B: -rt(A@B)
def forms(theta):
    """12x12 forms in (Y1..Y4) coordinates: f_p and the four corners (as v^T C v = 2<.,.>)."""
    W=np.diag([np.exp(1j*theta),np.exp(-1j*theta)]); Wd=W.conj().T
    f=np.zeros((4,3,4,3))
    for i in range(4):
        for a in range(3):
            for b in range(3):
                f[i,a,i,b]+=0.5*(rt(T[a]@T[b]@W)+rt(T[b]@T[a]@W))
    for (i,j) in [(0,1),(2,3)]:
        for a in range(3):
            for b in range(3):
                v=rt(T[a]@T[b]@W); f[i,a,j,b]+=v; f[j,b,i,a]+=v
    for i in (0,1):
        for j in (2,3):
            for a in range(3):
                for b in range(3):
                    v=-rt(T[a]@W@T[b]); f[i,a,j,b]+=v; f[j,b,i,a]+=v
    C={}
    for (i,j) in [(0,3),(0,1),(2,3)]:
        c=np.zeros((4,3,4,3))
        for a in range(3): c[i,a,j,a]=1; c[j,a,i,a]=1
        C[(i,j)]=c
    c=np.zeros((4,3,4,3))
    for a in range(3):
        for b in range(3):
            v=ip(T[a],W@T[b]@Wd); c[1,a,2,b]+=v; c[2,b,1,a]+=v
    C[(1,2)]=c
    return f.reshape(12,12),{k:v.reshape(12,12) for k,v in C.items()}
# plaquettes: base residue r, plane (mu<nu); links e1=(r,mu) e2=(r+mu,nu) e3=(r+nu,mu) e4=(r,nu)
plaqs=[]
for r in res:
    for mu,nu in itertools.combinations(range(4),2):
        L=[(r,mu),(add(r,mu),nu),(add(r,nu),mu),(r,nu)]
        F=[i for i,(s,m) in enumerate(L) if sig(s)==m]; assert len(F)==1
        plaqs.append(dict(r=r,mu=mu,nu=nu,links=L,i0=F[0]))
# corner -> (vertex residue, end_a, end_b); ends are (direction, 'o'|'i')
def corner_vertex(p,k):
    r,mu,nu=p['r'],p['mu'],p['nu']
    if k==(0,3): return r,(mu,'o'),(nu,'o')
    if k==(0,1): return add(r,mu),(mu,'i'),(nu,'o')
    if k==(2,3): return add(r,nu),(mu,'o'),(nu,'i')
    if k==(1,2): return add(add(r,mu),nu),(nu,'i'),(mu,'i')
def solve(nth=48,solver='CLARABEL',fixed=None):
    ths=np.linspace(0,2*np.pi,nth,endpoint=False)
    FS=[forms(t) for t in ths]
    h=cp.Variable()
    alpha={}; lam={}; gam={}
    cons=[]
    for pi,p in enumerate(plaqs):
        act=[i for i in range(4) if i!=p['i0']]
        for i in act: alpha[(pi,i)]=cp.Variable()
        for k in [(0,3),(0,1),(2,3),(1,2)]:
            if p['i0'] in k: continue
            lam[(pi,k)]=cp.Variable()
    for r in res:
        for mu in range(4):
            if mu==sig(r): continue
            for t in 'oi': gam[(r,mu,t)]=cp.Variable()
    # plaquette LMIs
    for pi,p in enumerate(plaqs):
        act=[i for i in range(4) if i!=p['i0']]
        idx=np.array([3*i+a for i in act for a in range(3)])
        for (f,C) in FS:
            S=0
            for i in act:
                D=np.zeros((12,12)); D[3*i:3*i+3,3*i:3*i+3]=np.eye(3)
                S=S+alpha[(pi,i)]*D[np.ix_(idx,idx)]
            S=S-f[np.ix_(idx,idx)]
            for k,c in C.items():
                if (pi,k) in lam: S=S-lam[(pi,k)]*c[np.ix_(idx,idx)]
            cons.append(0.5*(S+S.T)>>0)
    # vertex LMIs
    for r in res:
        ends=[(mu,t) for mu in range(4) if mu!=sig(r) for t in 'oi']
        E={e:j for j,e in enumerate(ends)}
        ent={}
        for pi,p in enumerate(plaqs):
            for k in [(0,3),(0,1),(2,3),(1,2)]:
                if (pi,k) not in lam: continue
                v,a,b=corner_vertex(p,k)
                if v!=r: continue
                ent[(E[a],E[b])]=lam[(pi,k)]
        M=[[0]*6 for _ in range(6)]
        for j,e in enumerate(ends): M[j][j]=gam[(r,e[0],e[1])]
        for (a,b),v in ent.items(): M[a][b]=v; M[b][a]=v
        Mx=cp.bmat([[M[a][b] if not isinstance(M[a][b],int) else cp.Constant(0) for b in range(6)] for a in range(6)])
        cons.append(Mx>>0)
    # link constraints
    for r in res:
        for mu in range(4):
            if mu==sig(r): continue
            tot=gam[(r,mu,'o')]+gam[(add(r,mu),mu,'i')]
            for pi,p in enumerate(plaqs):
                for i,(s,m) in enumerate(p['links']):
                    if (s,m)==(r,mu) and i!=p['i0']: tot=tot+alpha[(pi,i)]
            cons.append(tot<=h)
    prob=cp.Problem(cp.Minimize(h),cons)
    prob.solve(solver=solver,verbose=False)
    return prob.status,h.value,alpha,lam,gam
if __name__=="__main__":
    nth=int(sys.argv[1]) if len(sys.argv)>1 else 48
    st,hv,alpha,lam,gam=solve(nth)
    print('status',st,'h',hv)
    import pickle
    pickle.dump(({k:v.value for k,v in alpha.items()},{k:v.value for k,v in lam.items()},{k:v.value for k,v in gam.items()},hv),open(f'cert_{nth}.pkl','wb'))
