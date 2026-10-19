"""Rigorous bound of 2 r sup_h L_h over a box of diagonal staple sums b = diag(a1, a2, x3 + i y3) with the
degree-2 corrector u = q_h = sum_k C_kh B_k (B_k a real basis of polynomials of degree 1, 2 in U, Ubar):
L_h <= ||grad q_h|| + sqrt(Var(phi_h + L q_h) / gap), gap = 283/128 - r/2.
C is chosen per box in floating point and rounded to 2^-30; from then on every quantity is bounded rigorously."""
import numpy as np, itertools, flint, time
import rs6z
from polys import *
from rseries import gam, U_
arb=flint.arb
G40=gam(60)
def gk(k): return 1.5*gam(2*k+10)
def mm(A,eA,B,eB):
    aA=np.abs(A); aB=np.abs(B); k=A.shape[-1]
    V=A@B; E=aA@eB+eA@aB+eA@eB+gk(k)*(aA@aB)
    return V,E
def ldl_pd(As,t):
    n=len(As); M=[[arb(t if i==j else 0.0)-arb(float(As[i][j])) for j in range(n)] for i in range(n)]
    for k in range(n):
        if not M[k][k]>0: return False
        for i in range(k+1,n):
            f=M[i][k]/M[k][k]
            for j in range(k+1,n): M[i][j]=M[i][j]-f*M[k][j]
    return True
def lam(V,E):
    """certified upper bound of the largest eigenvalue of the true symmetric matrix within E of V"""
    As=(V+V.T)/2; Es=(E+E.T)/2+2*U_*np.abs(As)
    l=np.linalg.eigvalsh(As).max(); t=l+1e-12*max(1,abs(l))
    while not ldl_pd(As,t): t+=1e-10*max(1,abs(l))
    return (t+np.sqrt((Es**2).sum()))*(1+1e-14)+1e-300
def q30(x): return np.round(np.asarray(x)*2.0**30)/2.0**30
def corner2(b4):
    a1,a2,x3,y3=b4
    w,ew=rs6z.moments_rig([a1,a2,complex(x3,y3)])
    Z0=w[0].real; eZ0=ew[0]+abs(w[0].imag)
    Zlo,Zhi=Z0-eZ0,Z0+eZ0
    wn=w/Z0; ewn=ew/Zlo+np.abs(w)*eZ0/(Z0*Zlo)+2*U_*np.abs(wn)
    F=SR.full
    return dict(b=b4,Zlo=Zlo,Zhi=Zhi,w=wn[F],e=ewn[F],maxerr=float(ew.max()))
def full(cd):
    w=np.zeros(R.N,complex); e=np.zeros(R.N); w[SR.full]=cd['w']; e[SR.full]=cd['e']; return w,e
GAPC=283/128
def choose_C(wbar,b0,r,kaps=np.geomspace(0.05,80,48)):
    gap=GAPC-r/2
    W3=wbar[P3]; W2=W3[:N2,:N2]
    Dir=sum((w/6)*(Db[:,:N2]@W2@Db[:,:N2].T).real for Db,w in zip(DB,W6))
    LB=LAPB+sum(b0[t]*GTB[t] for t in range(4))
    F=np.vstack([PHI2/2,LB]); G=(F@W3@F.T).real; mu=(F@wbar[:N3]).real; Cov=G-np.outer(mu,mu)
    Cpp,CpL,CLL=Cov[:18,:18],Cov[:18,18:],Cov[18:,18:]
    best=(9,None,None)
    for kap in kaps:
        A=Dir+kap*CLL/gap
        C=-np.linalg.lstsq(A+1e-11*np.trace(A)/len(A)*np.eye(len(A)),kap*CpL.T/gap,rcond=None)[0]
        Aq=C.T@Dir@C; Bq=Cpp+CpL@C+C.T@CpL.T+C.T@CLL@C
        v=np.sqrt(max(np.linalg.eigvalsh((Aq+Aq.T)/2).max(),0))+np.sqrt(max(np.linalg.eigvalsh((Bq+Bq.T)/2).max(),0)/gap)
        if v<best[0]: best=(v,kap,C)
    return best[2],best[1],2*r*best[0]
def box_bound2(lo,hi,cds,r,verbose=False):
    lo=np.asarray(lo,float); hi=np.asarray(hi,float); wd=hi-lo; b0=(lo+hi)/2
    corners=[tuple(np.where(np.array(c)==1,hi,lo)) for c in itertools.product((0,1),repeat=4)]
    CD=[cds[c] for c in corners]
    WE=[full(cd) for cd in CD]
    wbar=np.mean([w for w,e in WE],0)
    Cf,kap,fest=choose_C(wbar,b0,r)
    Cint=np.round(Cf*2.0**30); P=scaled(Cint); s=P['s']
    Gs=[g/s for g in P['G']]
    for g in Gs: assert not np.any(g[:,N2:])
    Gs=[g[:,:N2] for g in Gs]
    d12=12*s
    Ft=[f/d12 for f in P['Ft']]; eFt=[gam(2)*np.abs(f) for f in Ft]
    lamA=0.0; lamV=0.0; mus=[]; data=[]
    for (w,e),bc in zip(WE,corners):
        W3=w[P3]; E3=e[P3]; W2=W3[:N2,:N2]; E2=E3[:N2,:N2]
        SV=0; SE=0; SA=0
        for g,w6 in zip(Gs,W6):
            T,eT=mm(g,0*g.real,W2,E2); V,E=mm(T,eT,g.T,0*g.real.T)
            SV=SV+w6*V.real; SE=SE+w6*E; SA=SA+w6*np.abs(V)
        A=SV/6; eA=(SE+gam(24)*SA)/6
        lamA=max(lamA,lam(A,eA))
        Num=P['F0']+sum(bc[t]*P['Ft'][t] for t in range(4)); Na=np.abs(P['F0'])+sum(abs(bc[t])*np.abs(P['Ft'][t]) for t in range(4))
        Fc=Num/d12; eF=gam(12)*Na/d12
        w3=w[:N3]; e3=e[:N3]
        mu=(Fc@w3).real; emu=np.abs(Fc)@e3+eF@np.abs(w3)+eF@e3+gk(N3)*(np.abs(Fc)@np.abs(w3))
        T,eT=mm(W3,E3,Fc.T,eF.T); Sg,eSg=mm(Fc,eF,T,eT)
        mus.append(mu); data.append((W3,E3,w3,e3,Sg.real,eSg,mu,emu))
    Dm=np.array([[1.0]+list(np.array(bc)-b0) for bc in corners])
    coef=np.linalg.lstsq(Dm,np.array(mus),rcond=None)[0]
    alpha=q30(coef[0]); gsl=[q30(coef[1+k]) for k in range(4)]
    def centered(S,eS,mu,emu,kk):
        V=S-np.outer(mu,kk)-np.outer(kk,mu)+np.outer(kk,kk)
        E=eS+np.outer(emu,np.abs(kk))+np.outer(np.abs(kk),emu)+G40*(np.abs(S)+2*np.abs(np.outer(mu,kk))+np.outer(np.abs(kk),np.abs(kk)))
        return V,E
    for (W3,E3,w3,e3,Sg,eSg,mu,emu),bc in zip(data,corners):
        kk=alpha+sum((bc[k]-b0[k])*gsl[k] for k in range(4))
        lamV=max(lamV,lam(*centered(Sg,eSg,mu,emu,kk)))
    eN=arb(0); worst_t=[]
    for t in range(4):
        worst=0.0
        for (W3,E3,w3,e3,_,_,_,_),cd in zip(data,CD):
            f=Ft[t]; ef=eFt[t]
            mu=(f@w3).real; emu=np.abs(f)@e3+ef@np.abs(w3)+ef@e3+gk(N3)*(np.abs(f)@np.abs(w3))
            T,eT=mm(W3,E3,f.T,ef.T); Sg,eSg=mm(f,ef,T,eT)
            worst=max(worst,cd['Zhi']*lam(*centered(Sg.real,eSg,mu,emu,gsl[t])))
        worst_t.append(worst); eN+=arb(float(wd[t]))**2/4*arb(worst)
    Zmin=min(cd['Zlo'] for cd in CD)
    sig=arb(float(wd[0]))**2+arb(float(wd[1]))**2+max(arb(float(wd[2]))**2,arb(float(wd[3]))**2)
    hf=(sig/8).exp()
    mQ=hf*arb(lamA); VQ=hf*(arb(lamV)+eN/arb(Zmin))
    B=mQ.sqrt()+(VQ/(arb(283)/128-arb(r)/2)).sqrt()
    val=2*arb(r)*B
    up=float(val.upper())*(1+2.0**-50)
    det=dict(kap=kap,fest=fest,lamA=lamA,lamV=lamV,eN=float(eN.mid()),hf=float(hf.mid()),bound=up,maxerr=max(cd['maxerr'] for cd in CD),exact_bound=P['bound'])
    if verbose: print(det,flush=True)
    return up,det
if __name__=='__main__':
    import sys
    r=1.5
    for lo,wdt in [((1.453125,1.078125,1.078125,0.0),0.046875),((1.40625,1.03125,1.03125,0.0),0.09375)]:
        lo=np.array(lo); hi=lo+wdt
        corners=[tuple(np.where(np.array(c)==1,hi,lo)) for c in itertools.product((0,1),repeat=4)]
        t=time.time(); cds={c:corner2(c) for c in corners}; t1=time.time()
        up,d=box_bound2(lo,hi,cds,r,verbose=True); print('box',lo,wdt,'corners',f'{t1-t:.1f}s','bound',f'{time.time()-t1:.1f}s',up,flush=True)
