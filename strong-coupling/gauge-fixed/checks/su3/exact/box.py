"""Rigorous bound of the one-link quantity over a box of diagonal staple sums b = diag(a1, a2, x3 + i y3)."""
import numpy as np, itertools, flint
from rseries import gam, U_
from corner import BLOCKS
arb=flint.arb
G40=gam(60)
def ldl_pd(As,t):
    n=len(As); M=[[arb(t if i==j else 0.0)-arb(float(As[i][j])) for j in range(n)] for i in range(n)]
    for k in range(n):
        if not M[k][k]>0: return False
        for i in range(k+1,n):
            f=M[i][k]/M[k][k]
            for j in range(k+1,n): M[i][j]=M[i][j]-f*M[k][j]
    return True
def lam(blks):
    """certified upper bound of the largest eigenvalue of the true symmetric block-diagonal matrix"""
    best=-np.inf
    for V,E in blks:
        As=(V+V.T)/2; Es=(E+E.T)/2+2*U_*np.abs(As)
        l=np.linalg.eigvalsh(As).max(); t=l+1e-12*max(1,abs(l))
        while not ldl_pd(As,t): t+=1e-10*max(1,abs(l))
        best=max(best,(t+np.sqrt((Es**2).sum()))*(1+1e-14)+1e-300)
    return best
def quadform(cd,c1,eta,k):
    """blocks of the quadratic form h -> nu_b[(c1 phi_h + sum_t eta_t Gamma_{t,h} - k.h)^2], k supported on block 0"""
    out=[]
    for bi in range(4):
        V=c1*c1*cd['Spp'][bi][0]; E=c1*c1*cd['Spp'][bi][1]; A=np.abs(V).copy()
        for t in range(4):
            if eta[t]==0: continue
            X,eX=cd['SpG'][t][bi]; X2=c1*eta[t]*(X+X.T); V=V+X2; E=E+abs(c1*eta[t])*(eX+eX.T); A=A+np.abs(c1*eta[t])*(np.abs(X)+np.abs(X.T))
            for u in range(4):
                if eta[u]==0: continue
                Y,eY=cd['SGG'][t][u][bi]; V=V+eta[t]*eta[u]*Y; E=E+abs(eta[t]*eta[u])*eY; A=A+abs(eta[t]*eta[u])*np.abs(Y)
        if bi==0:
            mu=c1*cd['muP'][0]; emu=abs(c1)*cd['muP'][1]; amu=np.abs(mu)
            for t in range(4):
                if eta[t]==0: continue
                mu=mu+eta[t]*cd['muG'][t][0]; emu=emu+abs(eta[t])*cd['muG'][t][1]; amu=amu+abs(eta[t])*np.abs(cd['muG'][t][0])
            emu=emu+G40*amu
            kk=np.asarray(k,float)
            O=-np.outer(mu,kk)-np.outer(kk,mu)+np.outer(kk,kk)
            V=V+O; E=E+np.outer(emu,np.abs(kk))+np.outer(np.abs(kk),emu); A=A+2*np.abs(np.outer(mu,kk))+np.outer(np.abs(kk),np.abs(kk))
        out.append((V,E+G40*A))
    return out
def q30(x): return np.round(np.asarray(x)*2.0**30)/2.0**30
GAPC=283/128
def box_bound(lo,hi,cds,r,j=None,verbose=False):
    """cds: dict corner tuple -> corner data. Returns (2 r B upper bound as float, details)."""
    lo=np.asarray(lo,float); hi=np.asarray(hi,float); w=hi-lo; b0=(lo+hi)/2
    corners=[tuple(np.where(np.array(c)==1,hi,lo)) for c in itertools.product((0,1),repeat=4)]
    C=[cds[c] for c in corners]
    gapf=GAPC-r/2
    # float choice of t = j/1024 (s = 3t/8)
    lamMf=max(max(np.linalg.eigvalsh((V+V.T)/2).max() for V,E in cd['M']) for cd in C)
    def Bf(jj):
        t=jj/1024; s=3*t/8; c1=1-t; vmax=0
        for cd,bc in zip(C,corners):
            Q=quadform(cd,c1,[s*x for x in bc],cd['muP'][0]*c1+sum(s*bc[tt]*cd['muG'][tt][0] for tt in range(4)))
            vmax=max(vmax,max(np.linalg.eigvalsh((V+V.T)/2).max() for V,E in Q))
        return s*np.sqrt(lamMf)+np.sqrt(vmax/gapf)
    if j is None:
        cand=list(range(640,1100,24)); vals=[Bf(x) for x in cand]; j0=cand[int(np.argmin(vals))]
        cand=list(range(j0-24,j0+25,4)); vals=[Bf(x) for x in cand]; j=cand[int(np.argmin(vals))]
    t=j/1024; s=3*j/8192; c1=1-j/1024
    # affine centering fitted to corner means
    mus=np.array([cd['muP'][0]*c1+sum(s*bc[tt]*cd['muG'][tt][0] for tt in range(4)) for cd,bc in zip(C,corners)])
    Dm=np.array([[1.0]+list(np.array(bc)-b0) for bc in corners])
    coef,_,_,_=np.linalg.lstsq(Dm,mus,rcond=None)
    alpha=q30(coef[0]); gk=[q30(coef[1+k]) for k in range(4)]
    def cvec(bc): return alpha+sum((bc[k]-b0[k])*gk[k] for k in range(4))
    lamM=max(lam(cd['M']) for cd in C)
    lamV=max(lam(quadform(cd,c1,[s*x for x in bc],cvec(bc))) for cd,bc in zip(C,corners))
    Zmin=min(cd['Zlo'] for cd in C)
    # Hoeffding: Z(b) >= exp(-sig/8) sum_i w_i Z_i, the coefficients of Re Tr(Ub) in (a1, a2, x3, y3) being
    # Re U11, Re U22, Re U33, -Im U33 (the last two with squares adding to at most 1)
    sig=arb(float(w[0]))**2+arb(float(w[1]))**2+max(arb(float(w[2]))**2,arb(float(w[3]))**2)
    hf=(sig/8).exp()
    # interpolation error of the numerator: -d_k^2 [(A + sum delta_l B_l)^2 e^{<U,b>}] <= 2 B_k^2 e^{<U,b>}
    eN=arb(0)
    for k in range(4):
        worst=0.0
        for cd in C:
            Qk=quadform(cd,0.0,[s if l==k else 0.0 for l in range(4)],gk[k])
            worst=max(worst,cd['Zhi']*lam(Qk))
        eN+=arb(float(w[k]))**2/4*arb(worst)
    Zm=arb(Zmin)
    mQ=hf*arb(lamM); VQ=hf*(arb(lamV)+eN/Zm)
    B=arb(s)*mQ.sqrt()+(VQ/(arb(283)/128-arb(r)/2)).sqrt()
    val=2*arb(r)*B
    up=float(val.upper())*(1+2.0**-50)
    hoeffding=hf
    if verbose: print(dict(j=j,s=s,lamM=lamM,lamV=lamV,Zmin=Zmin,hoeffding=float(hoeffding.mid()),eN=float(eN.mid()),bound=up))
    return up,dict(j=j,s=s,lamM=lamM,lamV=lamV,eN=float(eN.mid()),hoeffding=float(hoeffding.mid()))
