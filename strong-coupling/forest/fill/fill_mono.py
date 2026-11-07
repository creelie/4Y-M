"""Calm degree-N reducible collar on the boundary of a 3D cube (heuristic).
Boundary sites carry n_x in S^2 = composition of radial projection with z -> z^N (stereographic).
Link U_xy = exp(i th_xy n_y.sigma/2) Rot(n_x -> n_y) so that U (n_x.s) U^-1 = n_y.s (covariantly constant n).
Every boundary plaquette is then a rotation about n_x by angle curl(th) + Omega_p, with Omega_p the lattice solid angle.
We solve curl th = kappa - Omega_p (least squares on the closed surface; kappa fixed by the total-flux constraint mod 2pi)
-> uniform field kappa n.sigma, total flux 4 pi N. Interior: Haar start, cooling. Report max plaquette, boundary vs interior."""
import sys, json, numpy as np
rng=np.random.default_rng(5)
def q2m(q): return np.array([[q[0]+1j*q[1],q[2]+1j*q[3]],[-q[2]+1j*q[3],q[0]-1j*q[1]]])
S=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.array([[1,0],[0,-1]],complex)]
def ndot(n): return sum(a*s for a,s in zip(n,S))
def rotq(axis,ang): # SU(2) element exp(-i ang/2 axis.sigma) -> rotates by ang about axis in adjoint
    return np.cos(ang/2)*np.eye(2)-1j*np.sin(ang/2)*ndot(axis)
def rot_to(a,b):
    c=np.cross(a,b); s=np.linalg.norm(c); d=np.dot(a,b)
    if s<1e-12: return np.eye(2) if d>0 else rotq(np.array([1,0,0]) if abs(a[0])<.9 else np.array([0,1,0]),np.pi)
    return rotq(c/s,np.arctan2(s,d))
def proj(M):
    q=np.array([M[0,0].real,M[0,0].imag,M[0,1].real,M[0,1].imag]); return q2m(q/np.linalg.norm(q))
def angle_about(P,n): # P = exp(-i a/2 n.s) -> a
    return -2*np.arctan2(np.real(np.trace(P@ndot(n)))/(2j).imag*(-1) if False else -np.imag(np.trace(P@ndot(n)))/2*(-1), np.real(np.trace(P))/2)
def run(R,N,sweeps):
    n=R+1; E=[np.eye(3,dtype=int)[m] for m in range(3)]
    bd=lambda x: min(x)==0 or max(x)==R
    def nvec(x):
        c=np.array(x,float)-R/2; c/=np.linalg.norm(c)
        z=(c[0]+1j*c[1])/(1-c[2]) if c[2]<1-1e-12 else None
        if z is None: return np.array([0,0,1.])
        w=z**N; r2=abs(w)**2
        return np.array([2*w.real,2*w.imag,r2-1])/(r2+1)
    sites=[x for x in np.ndindex(n,n,n)]
    nv={x:nvec(x) for x in sites if bd(x)}
    links=[(x,m) for x in sites for m in range(3) if x[m]<R]
    def isb(x,m):
        y=tuple(np.add(x,E[m])); return bd(x) and bd(y) and any(k!=m and (x[k]==0 or x[k]==R) for k in range(3))
    blinks=[l for l in links if isb(*l)]; bidx={l:i for i,l in enumerate(blinks)}
    # boundary plaquettes
    bpl=[]
    for x in sites:
        for m in range(3):
            for nu in range(m+1,3):
                xm=tuple(np.add(x,E[m])); xn=tuple(np.add(x,E[nu]))
                if xm[m]>R or xn[nu]>R: continue
                ls=[(x,m),(xm,nu),(xn,m),(x,nu)]
                if all(isb(*l) for l in ls): bpl.append(ls)
    def base(l):
        x,m=l; y=tuple(np.add(x,E[m])); return rot_to(nv[x],nv[y])
    def plaq_mat(Uf,ls):
        (x,m),(xm,nu),(xn,_),_=ls
        return Uf(ls[0])@Uf(ls[1])@Uf(ls[2]).conj().T@Uf(ls[3]).conj().T
    # Omega_p: angle of plaquette built from base rotations, about n_x (holonomy fixes n_x)
    def ang(P,nx): return 2*np.arctan2(-np.imag(np.trace(P@ndot(nx)))/2, np.real(np.trace(P))/2)
    Om=np.array([ang(plaq_mat(base,ls).conj().T,nv[ls[0][0]]) for ls in bpl])  # convention fixed below numerically
    # curl matrix: theta enters as rotation about n_y on link; contributes +th for forward links, -th for backward (n-frame transported)
    Cm=np.zeros((len(bpl),len(blinks)))
    for i,ls in enumerate(bpl):
        for k,l in enumerate(ls): Cm[i,bidx[l]]+= 1 if k<2 else -1
    tot=Om.sum(); P=len(bpl)
    def build(th):
        def Uf(l): x,m=l; y=tuple(np.add(x,E[m])); return rotq(nv[y],th[bidx[l]])@base(l)
        return Uf
    # solve uniform target kappa; total curl = 0 on closed surface, so need sum(kappa - Om) = 0 mod 2pi
    kap=-tot/P
    th=np.linalg.lstsq(Cm,np.full(P,kap)-0*Om,rcond=None)[0]
    # iterate: correct using measured angles (handles sign/convention errors)
    for it in range(30):
        Uf=build(th); a=np.array([ang(plaq_mat(Uf,ls).conj().T,nv[ls[0][0]]) for ls in bpl])
        err=(a-np.median(a)+np.pi)%(2*np.pi)-np.pi
        if abs(err).max()<1e-9: break
        th=th-np.linalg.lstsq(Cm,err,rcond=None)[0]
    Uf=build(th); U={l:Uf(l) for l in blinks}
    for l in links:
        if l not in U: U[l]=q2m((lambda q:q/np.linalg.norm(q))(rng.normal(size=4)))
    free=[l for l in links if l not in bidx]
    def staple(x,m):
        Sm=np.zeros((2,2),complex); y=tuple(np.add(x,E[m]))
        for nu in range(3):
            if nu==m: continue
            xn=tuple(np.add(x,E[nu]))
            if xn[nu]<=R: Sm+=U[(y,nu)]@U[(xn,m)].conj().T@U[(x,nu)].conj().T
            xm=tuple(np.subtract(x,E[nu]))
            if xm[nu]>=0:
                ym=tuple(np.subtract(y,E[nu])); Sm+=U[(ym,nu)].conj().T@U[(xm,m)].conj().T@U[(xm,nu)]
        return Sm
    for s in range(sweeps):
        for l in free: U[l]=proj(staple(*l).conj().T)
    bmax=0; imax=0; prof={}
    for x in sites:
        for m in range(3):
            for nu in range(m+1,3):
                xm=tuple(np.add(x,E[m])); xn=tuple(np.add(x,E[nu]))
                if xm[m]>R or xn[nu]>R: continue
                ls=[(x,m),(xm,nu),(xn,m),(x,nu)]
                v=np.linalg.norm(np.eye(2)-plaq_mat(lambda l:U[l],ls))
                if all(l in bidx for l in ls): bmax=max(bmax,v)
                else:
                    imax=max(imax,v); d=min(min(x),R-max(xm[m],xn[nu],*x)) if True else 0
                    prof[d]=max(prof.get(d,0),v)
    return dict(R=R,N=N,kappa_angle=float(abs(np.median(a))),flux_err=float(abs(err).max()),bd_max=bmax,int_max=imax,ratio=imax/bmax,depth_profile={int(k):float(v) for k,v in sorted(prof.items())})
sw=int(sys.argv[1])
for R in (4,6,8,10,12):
    for N in (1,2,4,8):
        if 4*np.pi*N/(6*R*R)>1.2: continue
        print(json.dumps(run(R,N,sw)),flush=True)
