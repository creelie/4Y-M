"""Rigorous corner data: tilted moments at b = diag(a1, a2, x3 + i y3) and the block-diagonal quadratic forms
needed by the box bound, each as (value, absolute error bound)."""
import numpy as np, math
from series import Ring
from rseries import moments_rig, gam, U_
R4=Ring(4); R2=Ring(2); N2=R2.N
PIDX=np.array([[R4.idx[tuple(x+y for x,y in zip(a,b))] for b in R2.mons] for a in R2.mons])
import os
_NPZ=os.path.join(os.path.dirname(os.path.abspath(__file__)),'epoly.npz')
if not os.path.exists(_NPZ):
    import epoly; epoly.save(_NPZ)
_ep=np.load(_NPZ)
PHIv,PHIe,GGv,GGe,GAMv,GAMe=[_ep[k] for k in ('PHIv','PHIe','GGv','GGe','GAMv','GAMe')]
BLOCKS=[[0,1,8,9,16,17],[2,3,6,7],[4,5,12,13],[10,11,14,15]]
def gk(k): return 1.5*gam(2*k+10)
def mm(A,eA,B,eB):
    """complex matmul with error bound"""
    aA=np.abs(A); aB=np.abs(B); k=A.shape[-1]
    V=A@B
    E=aA@eB+eA@aB+eA@eB+gk(k)*(aA@aB)
    return V,E
def blocks(V,E):
    return [(V[np.ix_(b,b)].copy(),E[np.ix_(b,b)].copy()) for b in BLOCKS]
def corner_data(b4,P=24,K=8):
    a1,a2,x3,y3=b4
    w,ew=moments_rig(R4,[a1,a2,complex(x3,y3)],P,K)
    Z0=w[0].real; eZ0=ew[0]+abs(w[0].imag)
    Zlo,Zhi=Z0-eZ0,Z0+eZ0
    wn=w/Z0
    ewn=(ew/Zlo+np.abs(w)*eZ0/(Z0*Zlo))+2*U_*np.abs(wn)
    W=wn[PIDX]; eW=ewn[PIDX]
    W0=W[0]; eW0=eW[0]
    # M_pq = Re sum_k GG_pqk W0_k
    Mv=np.einsum('pqk,k->pq',GGv,W0); Me=np.einsum('pqk,k->pq',np.abs(GGv),eW0)+np.einsum('pqk,k->pq',GGe,np.abs(W0))+np.einsum('pqk,k->pq',GGe,eW0)+gk(N2)*np.einsum('pqk,k->pq',np.abs(GGv),np.abs(W0))
    # means
    def vec(Av,Ae):
        V=Av@W0; E=np.abs(Av)@eW0+Ae@np.abs(W0)+Ae@eW0+gk(N2)*(np.abs(Av)@np.abs(W0)); return V.real,E
    muP=vec(PHIv,PHIe); muG=[vec(GAMv[t],GAMe[t]) for t in range(4)]
    # second moments S(A,B) = Re A W B^T
    def S2(Av,Ae,Bv,Be):
        T,eT=mm(W,eW,Bv.T,Be.T)
        V,E=mm(Av,Ae,T,eT); return V.real,E
    TP=mm(W,eW,PHIv.T,PHIe.T); TG=[mm(W,eW,GAMv[t].T,GAMe[t].T) for t in range(4)]
    def S2c(Av,Ae,T):
        V,E=mm(Av,Ae,T[0],T[1]); return V.real,E
    Spp=S2c(PHIv,PHIe,TP)
    SpG=[S2c(PHIv,PHIe,TG[t]) for t in range(4)]
    SGG=[[S2c(GAMv[t],GAMe[t],TG[u]) for u in range(4)] for t in range(4)]
    out=dict(b=b4,Zlo=Zlo,Zhi=Zhi,
             M=blocks(Mv.real,Me),
             muP=(muP[0][BLOCKS[0]],muP[1][BLOCKS[0]]),
             muG=[(m[0][BLOCKS[0]],m[1][BLOCKS[0]]) for m in muG],
             Spp=blocks(*Spp),SpG=[blocks(*x) for x in SpG],SGG=[[blocks(*x) for x in row] for row in SGG],
             maxerr=max(ew.max(),0))
    return out
if __name__=='__main__':
    import time
    t=time.time(); d=corner_data((1.36,1.02,1.02,0.0)); print(time.time()-t,'s', d['Zlo'],d['Zhi'], d['maxerr'])
    print('M block errs',[e.max() for v,e in d['M']])
    print('Spp errs',[e.max() for v,e in d['Spp']], 'SGG errs',max(e.max() for row in d['SGG'] for blk in row for v,e in blk))
