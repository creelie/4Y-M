"""Independent numerical checks of the proof of Theorem A (sharp Hessian bound 4d).
(a) local lemma: F_p <= 2 on g^4 for random SU(N) plaquettes (brute-force Hessian by exact insertion,
    corner terms built from the covariant divergence with link transports), N=2,3,4, both signs;
(b) global identity/inequality on tori: lambda_max(+-Hess S0 + D^T D) <= 4d  for random and
    centre-frustrated SU(N) fields; and lambda_max(Hess S0) = 4d at the frustrated field.
Metric <X,Y> = -Tr(XY) (SZZ)."""
import numpy as np, itertools, sys
rng=np.random.default_rng(12345)
def su_basis(N):
    B=[]
    for a in range(N):
        for b in range(a+1,N):
            E=np.zeros((N,N),complex); E[a,b]=1; E[b,a]=-1; B.append(E/np.sqrt(2))
            E=np.zeros((N,N),complex); E[a,b]=1j; E[b,a]=1j; B.append(E/np.sqrt(2))
    for k in range(1,N):
        D=np.zeros(N); D[:k]=1; D[k]=-k; D=D/np.sqrt(k*(k+1)); B.append(1j*np.diag(D))
    return B
def haar(N):
    A=rng.normal(size=(N,N))+1j*rng.normal(size=(N,N)); q,r=np.linalg.qr(A); q=q@np.diag(np.diag(r)/abs(np.diag(r)))
    return q/np.linalg.det(q)**(1/N)
def Admat(g,B):
    return np.array([[-np.trace(Bi@g@Bj@g.conj().T).real for Bj in B] for Bi in B])
def plaq_hess(Qs,B):
    """exact Hessian of ReTr(Q1 Q2 Q3^* Q4^*) along Q_e -> e^{tX_e} Q_e, basis B per link."""
    Q1,Q2,Q3,Q4=Qs; n=len(B)
    fac=[('v',0,1),('c',Q1),('v',1,1),('c',Q2),('c',Q3.conj().T),('v',2,-1),('c',Q4.conj().T),('v',3,-1)]
    sgn={f[1]:f[2] for f in fac if f[0]=='v'}
    def prod(ins):
        M=np.eye(Q1.shape[0],dtype=complex)
        for f in fac:
            if f[0]=='c': M=M@f[1]
            elif f[1] in ins: M=M@ins[f[1]]
        return M
    H=np.zeros((4*n,4*n))
    for i in range(4):
        for a in range(n):
            for j in range(4):
                for b in range(n):
                    if i==j: v=np.trace(prod({i:0.5*(B[a]@B[b]+B[b]@B[a])})).real
                    else: v=np.trace(prod({i:sgn[i]*B[a],j:sgn[j]*B[b]})).real
                    H[i*n+a,j*n+b]=v
    return H
def corners(Qs,B):
    """corner cross terms of sum_x |D_x X|^2, D_x = sum_out X_e - sum_in Ad(Q_e^*) X_e."""
    Q1,Q2,Q3,Q4=Qs; n=len(B); C=np.zeros((4*n,4*n)); I=np.eye(n)
    T=[Admat(Q.conj().T,B) for Q in Qs]   # matrix of X -> Q^* X Q
    # vertex x0: out 1, out 4 ; x1: in 1, out 2 ; x2: in 2, in 3 ; x3: out 3, in 4
    terms=[(0,I,3,I),(0,-T[0],1,I),(1,-T[1],2,-T[2]),(2,I,3,-T[3])]
    for (i,Ai,j,Aj) in terms:
        M=Ai.T@Aj      # 2 <Ai X_i, Aj X_j>
        C[i*n:(i+1)*n,j*n:(j+1)*n]+=M; C[j*n:(j+1)*n,i*n:(i+1)*n]+=M.T
    return C
def check_local(N,trials):
    B=su_basis(N); worst=-1e9
    for _ in range(trials):
        Qs=[haar(N) for _ in range(4)]
        H=plaq_hess(Qs,B); C=corners(Qs,B)
        worst=max(worst,np.linalg.eigvalsh(H+C)[-1],np.linalg.eigvalsh(-H+C)[-1])
    return worst
class Torus:
    def __init__(s,d,L):
        s.d=d; s.L=L; s.sites=list(itertools.product(range(L),repeat=d)); s.idx={}
        for x in s.sites:
            for mu in range(d): s.idx[(x,mu)]=len(s.idx)
        sh=lambda x,mu,k: tuple((x[i]+k*(i==mu))%L for i in range(d))
        s.sh=sh
        s.plaq=[(s.idx[(x,mu)],s.idx[(sh(x,mu,1),nu)],s.idx[(sh(x,nu,1),mu)],s.idx[(x,nu)]) for x in s.sites for mu in range(d) for nu in range(mu+1,d)]
def global_check(N,d,L,Q):
    T=Torus(d,L); B=su_basis(N); n=len(B); nE=len(T.idx)
    H=np.zeros((n*nE,n*nE))
    for p in T.plaq:
        Hp=plaq_hess([Q[e] for e in p],B); ids=np.concatenate([np.arange(e*n,e*n+n) for e in p])
        H[np.ix_(ids,ids)]+=Hp
    D=np.zeros((n*len(T.sites),n*nE))
    for k,x in enumerate(T.sites):
        for mu in range(d):
            e=T.idx[(x,mu)]; D[k*n:(k+1)*n,e*n:(e+1)*n]+=np.eye(n)
            ein=T.idx[(T.sh(x,mu,-1),mu)]; D[k*n:(k+1)*n,ein*n:(ein+1)*n]-=Admat(Q[ein].conj().T,B)
    DD=D.T@D
    return np.linalg.eigvalsh(H)[[0,-1]], np.linalg.eigvalsh(H+DD)[-1], np.linalg.eigvalsh(-H+DD)[-1]
if __name__=="__main__":
    for N,tr in [(2,400),(3,150),(4,40)]:
        print(f"(a) N={N}: max over {tr} random plaquettes of lambda_max(+-H_p + C_p) = {check_local(N,tr):.12f}  (claim <= 2)")
    for N,d,L in [(2,3,3),(3,3,3),(2,4,3)]:
        T=Torus(d,L); Q=[haar(N) for _ in T.idx]
        ev,a,b=global_check(N,d,L,Q)
        print(f"(b) random SU({N}) d={d} L={L}: Hess range {ev}, max(H+D'D)={a:.10f}, max(-H+D'D)={b:.10f}, 4d={4*d}")
    for N,d,L in [(2,4,2),(3,3,2)]:
        T=Torus(d,L); z=np.eye(N,dtype=complex); z[0,0]=z[1,1]=-1
        Q=[ (z if sum(x[:mu])%2 else np.eye(N,dtype=complex)) for (x,mu) in T.idx]
        ev,a,b=global_check(N,d,L,Q)
        print(f"(b) frustrated SU({N}) d={d} L={L}: Hess range {ev}, max(H+D'D)={a:.10f}, 4d={4*d}")
