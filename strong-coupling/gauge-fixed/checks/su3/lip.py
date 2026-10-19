"""Single-link Lipschitz rate for SU(3), nu_b ~ exp(ReTr(U b)) dU.
Cauchy-Schwarz bound: L(b) <= sqrt(lambda_max M), M_ij = <f_i, (-L)^{-1} f_j>, f = real coords of U.
Galerkin (lower approximation of M) with polynomial features of degree <= D. Importance sampling from Haar."""
import numpy as np, itertools, sys
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
N=3
def haar(n):
    Z=(rng.standard_normal((n,N,N))+1j*rng.standard_normal((n,N,N)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z)
    d=np.diagonal(R,axis1=1,axis2=2); Q=Q*(d/np.abs(d))[:,None,:]
    det=np.linalg.det(Q); Q=Q/(det**(1/3))[:,None,None]
    return Q
lam=np.zeros((8,3,3),complex)
lam[0][0,1]=lam[0][1,0]=1
lam[1][0,1]=-1j; lam[1][1,0]=1j
lam[2][0,0]=1; lam[2][1,1]=-1
lam[3][0,2]=lam[3][2,0]=1
lam[4][0,2]=-1j; lam[4][2,0]=1j
lam[5][1,2]=lam[5][2,1]=1
lam[6][1,2]=-1j; lam[6][2,1]=1j
lam[7]=np.diag([1,1,-2])/np.sqrt(3)
T=1j*lam/np.sqrt(2)   # -Tr(T_a T_b)=delta
def coords(U):  # n x 18
    return np.concatenate([U.real.reshape(len(U),-1),U.imag.reshape(len(U),-1)],1)
def features(U,D):
    X=coords(U); dX=np.stack([coords(np.einsum('aij,njk->nik',T[a][None],U) if False else T[a]@U) for a in range(8)],0) # 8 x n x 18
    F=[X]; G=[dX]
    if D>=2:
        idx=[(p,q) for p in range(18) for q in range(p,18)]
        P=np.array([i[0] for i in idx]); Qi=np.array([i[1] for i in idx])
        F.append(X[:,P]*X[:,Qi]); G.append(dX[:,:,P]*X[None][:,:,Qi]+X[None][:,:,P]*dX[:,:,Qi])
    if D>=3:
        idx=[(p,q,s) for p in range(18) for q in range(p,18) for s in range(q,18)]
        P=np.array([i[0] for i in idx]); Qi=np.array([i[1] for i in idx]); S=np.array([i[2] for i in idx])
        F.append(X[:,P]*X[:,Qi]*X[:,S]); G.append(dX[:,:,P]*(X[:,Qi]*X[:,S])[None]+dX[:,:,Qi]*(X[:,P]*X[:,S])[None]+dX[:,:,S]*(X[:,P]*X[:,Qi])[None])
    return X,np.concatenate(F,1),np.concatenate(G,2)
def Lrate(b,D=2,n=200000,chunk=20000,U0=None):
    """returns sqrt(lambda_max(M)), lambda_max(Cov f), Galerkin."""
    S0=0;S1=None;S2=None;SA=None;SXX=None;SX=None;SXF=None
    for c in range(0,n,chunk):
        U=haar(chunk) if U0 is None else U0[c:c+chunk]
        X,F,G=features(U,D)
        w=np.exp(np.einsum('nij,ji->n',U,b).real)
        S0+=w.sum()
        sF=(w[:,None]*F).sum(0); sFF=(w[:,None]*F).T@F; Gw=(G*np.sqrt(w)[None,:,None]).reshape(-1,G.shape[2]); sA=Gw.T@Gw
        sX=(w[:,None]*X).sum(0); sXX=(w[:,None]*X).T@X; sXF=(w[:,None]*X).T@F
        if S1 is None: S1,S2,SA,SX,SXX,SXF=sF,sFF,sA,sX,sXX,sXF
        else: S1+=sF;S2+=sFF;SA+=sA;SX+=sX;SXX+=sXX;SXF+=sXF
    mF=S1/S0; A=SA/S0; mX=SX/S0
    CovX=SXX/S0-np.outer(mX,mX); CovXF=SXF/S0-np.outer(mX,mF)
    # M = CovXF A^+ CovXF^T
    ev,V=np.linalg.eigh(A); keep=ev>1e-9*ev.max()
    Ainv=(V[:,keep]/ev[keep])@V[:,keep].T
    M=CovXF@Ainv@CovXF.T
    return np.sqrt(np.linalg.eigvalsh(M).max()), np.linalg.eigvalsh(CovX).max()
if __name__=='__main__':
    for D in (1,2):
        print('b=0 D',D, Lrate(np.zeros((3,3)),D=D,n=100000))
