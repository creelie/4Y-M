"""Adversarial test of the local lemma S_p >= 0 over W in U(N) (and SU(N)), for every frozen position."""
import numpy as np, sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from verify_global import basis, rt
from scipy.optimize import minimize
from scipy.linalg import expm
a,b,t=2.1,1.25,0.9
opp={0:2,1:3,2:0,3:1}
def Smat(W,B,i0,sgn):
    dg=len(B); n=4*dg; M=np.zeros((n,n)); Wd=W.conj().T
    BB=np.einsum('aij,bjk,kl->abil',B,B,W); tr_BBW=np.real(np.einsum('abii->ab',BB))
    BWB=np.real(np.einsum('aij,jk,bki->ab',B,W,B))
    ipW=-np.real(np.einsum('aij,jk,bkl,li->ab',B,W,B,Wd))
    I=np.eye(dg); blk=lambda i,j: (slice(i*dg,(i+1)*dg),slice(j*dg,(j+1)*dg))
    f=np.zeros((n,n))
    for i in range(4): f[blk(i,i)]+=0.5*(tr_BBW+tr_BBW.T)
    for (i,j) in [(0,1),(2,3)]: f[blk(i,j)]+=tr_BBW; f[blk(j,i)]+=tr_BBW.T
    for i in (0,1):
        for j in (2,3): f[blk(i,j)]+=-BWB; f[blk(j,i)]+=-BWB.T
    C=np.zeros((n,n))
    cor={(0,3):I,(0,1):-I,(2,3):-I,(1,2):ipW}
    for (i,j),v in cor.items():
        if i0 in (i,j): continue
        C[blk(i,j)]+=v; C[blk(j,i)]+=v.T
    D=np.zeros((n,n))
    for i in range(4):
        if i!=i0: D[blk(i,i)]=(b if i==opp[i0] else a)*I
    S=D-sgn*f-t*C
    act=[k for i in range(4) if i!=i0 for k in range(i*dg,(i+1)*dg)]
    return S[np.ix_(act,act)]
def lam(W,B,i0,sgn): return np.linalg.eigvalsh(Smat(W,B,i0,sgn))[0]
def hermbasis(N):
    Hs=[]
    for i in range(N):
        for j in range(i,N):
            E=np.zeros((N,N),complex); E[i,j]=1; E[j,i]=1; Hs.append(E)
            if i!=j: E=np.zeros((N,N),complex); E[i,j]=1j; E[j,i]=-1j; Hs.append(E)
    return np.array(Hs)
def run(N,unitary,starts,seed):
    rng=np.random.default_rng(seed); B=basis(N,unitary); Hb=hermbasis(N); worst=np.inf
    for i0 in range(4):
        for sgn in (1,-1):
            best=np.inf
            for k in range(starts):
                p0=rng.normal(size=len(Hb))*rng.uniform(0.3,3)
                obj=lambda p: lam(expm(1j*np.einsum('a,aij->ij',p,Hb)),B,i0,sgn)
                r=minimize(obj,p0,method='Nelder-Mead',options=dict(maxiter=4000,xatol=1e-9,fatol=1e-12))
                best=min(best,r.fun)
            print(f'N={N} U={unitary} i0={i0} sign={sgn:+d}: min lambda_min = {best:.6f}',flush=True)
            worst=min(worst,best)
    print(f'N={N} U={unitary}: worst {worst:.6f}',flush=True)
if __name__=="__main__":
    # sanity: SU(2) diagonal W reproduces the exact pair/diagonal minimum 0.014922 found before
    B=basis(2,False)
    th=np.linspace(0,np.pi,2001)
    print('SU(2) diagonal scan min:',min(lam(np.diag([np.exp(1j*x),np.exp(-1j*x)]),B,i0,s) for x in th for i0 in range(4) for s in (1,-1)),flush=True)
    for (N,U,starts) in [(2,True,6),(3,False,6),(3,True,6),(4,True,4)]:
        run(N,U,starts,N*10+U)
