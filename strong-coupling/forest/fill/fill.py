"""Calm-fill test (heuristic). 3D SU(2) cube of side R (sites 0..R), boundary links fixed to
U = g_x V g_{x+mu}^-1 with g a smooth random gauge function (winding-free, S^2 boundary) and V = exp(kappa*noise).
Interior links: start Haar-random, then cool (each link -> normalised staple sum = exact local minimiser) N sweeps.
Report max interior plaquette distance |1-U_p| over plaquettes vs boundary max, as function of R."""
import sys, json, numpy as np
rng=np.random.default_rng(int(sys.argv[3]) if len(sys.argv)>3 else 0)
def q2m(q): return np.array([[q[0]+1j*q[1],q[2]+1j*q[3]],[-q[2]+1j*q[3],q[0]-1j*q[1]]])
def expsu(v):
    th=np.linalg.norm(v); n=v/th if th>0 else v
    return q2m(np.r_[np.cos(th),np.sin(th)*n])
def proj(M):  # nearest SU(2) to sum of SU(2): quaternion normalise
    q=np.array([M[0,0].real,M[0,0].imag,M[0,1].real,M[0,1].imag]); return q2m(q/np.linalg.norm(q))
def run(R,kappa,sweeps):
    n=R+1; E=[np.eye(3,dtype=int)[m] for m in range(3)]
    # smooth gauge: g(x)=exp(sum of low modes), amplitude O(1)
    a=rng.normal(size=(3,3,3))
    def g(x):
        v=np.array([sum(a[c,m,0]*np.sin(2*np.pi*(m+1)*x[m]/R+a[c,m,1]) for m in range(3)) for c in range(3)])*0.8
        return expsu(v)
    G={x:g(np.array(x)) for x in np.ndindex(n,n,n)}
    U={}
    bd=lambda x: min(x)==0 or max(x)==R
    links=[(x,m) for x in np.ndindex(n,n,n) for m in range(3) if x[m]<R]
    isb=lambda x,m: bd(x) and bd(tuple(np.add(x,E[m])))  and sum(1 for k in range(3) if k!=m and (x[k]==0 or x[k]==R))>0
    for x,m in links:
        y=tuple(np.add(x,E[m]))
        if isb(x,m):
            if MODE=='hedgehog':
                c=np.array(x,float)+0.5*E[m]-R/2
                k=[j for j in range(3) if j!=m and (x[j]==0 or x[j]==R)][0]; nrm=np.sign(c[k])*E[k]
                Am=(kappa/2)*np.cross(nrm,c)[m]   # abelian potential for field kappa along normal
                U[(x,m)]=expsu(Am*nrm.astype(float))
            else: U[(x,m)]=G[x]@expsu(kappa*rng.normal(size=3)/np.sqrt(3))@G[y].conj().T
        else: U[(x,m)]=q2m((lambda q:q/np.linalg.norm(q))(rng.normal(size=4)))
    free=[l for l in links if not isb(*l)]
    def L(x,m): return U[(x,m)]
    def staple(x,m):
        S=np.zeros((2,2),complex); y=tuple(np.add(x,E[m]))
        for nu in range(3):
            if nu==m: continue
            xn=tuple(np.add(x,E[nu]))
            if xn[nu]<=R:
                S+=L(y,nu)@L(xn,m).conj().T@L(x,nu).conj().T
            xm=tuple(np.subtract(x,E[nu]))
            if xm[nu]>=0:
                ym=tuple(np.subtract(y,E[nu]))
                S+=L(ym,nu).conj().T@L(xm,m).conj().T@L(xm,nu)
        return S
    for s in range(sweeps):
        for l in free: U[l]=proj(staple(*l).conj().T)
    def plaqs(inner):
        out=[]
        for x in np.ndindex(n,n,n):
            for m in range(3):
                for nu in range(m+1,3):
                    xm=tuple(np.add(x,E[m])); xn=tuple(np.add(x,E[nu]))
                    if xm[m]>R or xn[nu]>R: continue
                    P=L(x,m)@L(xm,nu)@L(xn,m).conj().T@L(x,nu).conj().T
                    on=all(isb(*l) for l in [(x,m),(xm,nu),(xn,m),(x,nu)])
                    if on!=inner: out.append(np.linalg.norm(np.eye(2)-P))
        return max(out) if out else 0.0
    return plaqs(False), plaqs(True)
kappa=float(sys.argv[1]); sw=int(sys.argv[2]); MODE=sys.argv[4] if len(sys.argv)>4 else 'noise'
for R in (3,4,6,8,10):
    b,i=run(R,kappa,sw); print(json.dumps({'R':R,'kappa':kappa,'bd_max':b,'int_max':i,'ratio':i/b}),flush=True)
