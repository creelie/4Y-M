"""Fast exact evaluation of the one-link bound at diagonal b = diag(a1, a2, x3 + i y3)."""
import numpy as np
from series import Ring, Zseries, moments
KAP=8/3
R4=Ring(4); R2=Ring(2); N2=R2.N
assert R4.mons[:N2]==R2.mons
PIDX=np.array([[R4.idx[tuple(x+y for x,y in zip(a,b))] for b in R2.mons] for a in R2.mons])
E=[]
for i in range(3):
    for j in range(3):
        M=np.zeros((3,3),complex); M[i,j]=1; E.append(M); M=np.zeros((3,3),complex); M[i,j]=1j; E.append(M)
FB=[np.diag([1,0,0]).astype(complex),np.diag([0,1,0]).astype(complex),np.diag([0,0,1]).astype(complex),np.diag([0,0,1j])]
R=R2
Um=[[R.var(3*i+j) for j in range(3)] for i in range(3)]
Ubm=[[R.var(9+3*i+j) for j in range(3)] for i in range(3)]
def lin(c,Mx): return [[sum(c[a,k]*Mx[k][b] for k in range(3)) for b in range(3)] for a in range(3)]
def tr(A): return A[0][0]+A[1][1]+A[2][2]
def trU(h): return tr(lin(h,Um))
def trUb(h): return tr(lin(np.conj(h),Ubm))
def phi(h): return (trU(h)+trUb(h))/2
def G(c,h):
    A=lin(c,Um);B=lin(h,Um);T=sum(R.mul(A[a][k],B[k][a]) for a in range(3) for k in range(3))
    A=lin(np.conj(c),Ubm);B=lin(np.conj(h),Ubm);Tb=sum(R.mul(A[a][k],B[k][a]) for a in range(3) for k in range(3))
    imc=(trU(c)-trUb(c))/2j; imh=(trU(h)-trUb(h))/2j
    return R.const(0.5*np.trace(c@h.conj().T).real)-(T+Tb)/4-R.mul(imc,imh)/3
PHI=np.array([phi(h) for h in E])                       # 18 x 190
GG=np.array([[G(E[p],E[q]) for q in range(18)] for p in range(18)])   # 18 x 18 x 190
GAM=np.array([[G(f,h) for h in E] for f in FB])         # 4 x 18 x 190
def evaluate(b4,P=18,K=8):
    """b4 = (a1,a2,x3,y3). Returns dict of exact (float) quantities."""
    a1,a2,x3,y3=b4; b=np.array([a1,a2,x3+1j*y3])
    Z=Zseries(R4,b,P=P,KMAX=K); Z0=Z[0].real; w=moments(R4,Z)/Z0
    W=w[PIDX]; W0=W[0]
    M=np.einsum('pqk,k->pq',GG,W0).real
    muP=(PHI@W0).real; muG=np.einsum('tpk,k->tp',GAM,W0).real
    Cpp=(PHI@W@PHI.T).real-np.outer(muP,muP)
    CpG=np.einsum('pk,kl,tql->tpq',PHI,W,GAM).real-np.einsum('p,tq->tpq',muP,muG)
    CGG=np.einsum('tpk,kl,uql->tupq',GAM,W,GAM).real-np.einsum('tp,uq->tupq',muG,muG)
    return dict(Z0=Z0,M=M,Cpp=Cpp,CpG=CpG,CGG=CGG,beta=np.array(b4,float))
def covX(d,s,beta=None):
    beta=d['beta'] if beta is None else beta
    c1=1-KAP*s
    CpGb=np.einsum('t,tpq->pq',beta,d['CpG']); CGGb=np.einsum('t,u,tupq->pq',beta,beta,d['CGG'])
    return c1*c1*d['Cpp']+c1*s*(CpGb+CpGb.T)+s*s*CGGb
def Bval(d,gap,svals=np.linspace(0.15,0.45,61)):
    m=np.linalg.eigvalsh(d['M']).max(); best=(9,)
    for s in svals:
        v=np.linalg.eigvalsh(covX(d,s)).max(); val=s*np.sqrt(m)+np.sqrt(v/gap)
        if val<best[0]: best=(val,s,m,v)
    return best
if __name__=='__main__':
    r=1.14; gap=283/128-r/2
    w3=np.exp(1j*np.pi/3)
    for b4 in [(r,r,r,0),(0,0,r,0),(r/2*np.cos(np.pi/3),)*0 or (0.5*r,0.5*r,0.5*r*np.cos(np.pi),0.5*r*np.sin(np.pi))]:
        import time; t=time.time(); d=evaluate(b4); print(b4,Bval(d,gap),time.time()-t)
