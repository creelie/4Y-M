"""#6 test: sign of d/dbeta log(Z-/Z) = <S>_twisted - <S>_periodic, S = sum_p z_p Re tr U_p / N.
Twist: z_p = -1 on (0,1)-plaquettes with x0=x1=L-1 (center flux through the (2,3) torus).
Groups: SU(2) (2x2), U(1) (1x1). Checkerboard Metropolis, numpy."""
import numpy as np, sys, json, itertools
from multiprocessing import Pool
D=4
def rand_near(rng, shape, N, eps):
    if N==1: return np.exp(1j*eps*rng.uniform(-1,1,shape))[...,None,None]
    a=rng.normal(size=shape+(3,))*eps; th=np.linalg.norm(a,axis=-1,keepdims=True); n=a/np.maximum(th,1e-12)
    c=np.cos(th[...,0]); s=np.sin(th[...,0])
    M=np.empty(shape+(2,2),complex)
    M[...,0,0]=c+1j*s*n[...,2]; M[...,0,1]=s*(n[...,1]+1j*n[...,0])
    M[...,1,0]=s*(-n[...,1]+1j*n[...,0]); M[...,1,1]=c-1j*s*n[...,2]
    return M
dag=lambda A: np.conj(np.swapaxes(A,-1,-2))
sh=lambda A,mu,s: np.roll(A,-s,axis=mu)   # A(x+s*mu)
def run(args):
    grp,beta,L,twist,nsw,seed=args; N=2 if grp=="SU2" else 1; rng=np.random.default_rng(seed)
    U=np.zeros((D,)+(L,)*D+(N,N),complex); U[...]=np.eye(N)
    z=np.ones((D,D)+(L,)*D)
    if twist:
        m=np.zeros((L,)*D,bool); m[L-1,L-1,:,:]=True; z[0,1][m]=-1; z[1,0][m]=-1
    par=np.indices((L,)*D).sum(0)%2
    def staple(mu):
        A=np.zeros((L,)*D+(N,N),complex)
        for nu in range(D):
            if nu==mu: continue
            f=sh(U[nu],mu,1)@dag(sh(U[mu],nu,1))@dag(U[nu])
            b=dag(sh(sh(U[nu],mu,1),nu,-1))@dag(sh(U[mu],nu,-1))@sh(U[nu],nu,-1)
            A+=z[mu,nu][...,None,None]*f+sh(z[mu,nu],nu,-1)[...,None,None]*b
        return A
    def action():
        S=0.0
        for mu in range(D):
            for nu in range(mu+1,D):
                P=U[mu]@sh(U[nu],mu,1)@dag(sh(U[mu],nu,1))@dag(U[nu])
                S+=(z[mu,nu]*np.trace(P,axis1=-2,axis2=-1).real).sum()/N
        return S
    eps=0.5 if grp=="SU2" else 1.0; meas=[]; acc=[]
    for sw in range(nsw):
        for mu in range(D):
            for p in (0,1):
                A=staple(mu); msk=(par==p)
                for h in range(3):
                    R=rand_near(rng,(L,)*D,N,eps); Un=R@U[mu]
                    dS=beta/N*(np.trace(Un@A,axis1=-2,axis2=-1).real-np.trace(U[mu]@A,axis1=-2,axis2=-1).real)
                    ok=msk&(rng.random((L,)*D)<np.exp(np.minimum(dS,0)))
                    U[mu][ok]=Un[ok]; acc.append(ok.sum()/msk.sum())
        if sw>=nsw//5: meas.append(action())
    return np.array(meas), float(np.mean(acc))
def jk(x,nb=20):
    b=np.array_split(x,nb); m=np.array([np.mean(np.concatenate(b[:i]+b[i+1:])) for i in range(nb)])
    return x.mean(), np.sqrt((nb-1)*m.var())
if __name__=="__main__":
    L=int(sys.argv[1]); nsw=int(sys.argv[2]); reps=4
    cases=[("SU2",b) for b in (1.5,2.0,2.3,2.6)]+[("U1",b) for b in (0.8,0.95,1.1,1.4)]
    jobs=[(g,b,L,t,nsw,1000*i+10*r+t) for i,(g,b) in enumerate(cases) for t in (0,1) for r in range(reps)]
    with Pool(16) as P: out=P.map(run,jobs)
    k=0
    for g,b in cases:
        res={}
        for t in (0,1):
            xs=[out[k+r][0] for r in range(reps)]; acc=np.mean([out[k+r][1] for r in range(reps)]); k+=reps
            ms=[jk(x) for x in xs]; res[t]=(np.mean([m for m,_ in ms]), np.sqrt(sum(e*e for _,e in ms))/reps, acc)
        d=res[1][0]-res[0][0]; e=np.hypot(res[1][1],res[0][1])
        print(json.dumps({"group":g,"beta":b,"L":L,"S_per":res[0][0],"S_tw":res[1][0],"Delta":d,"err":e,"z":d/e,"acc":[res[0][2],res[1][2]]}),flush=True)
