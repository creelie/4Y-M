"""Floating-point checks of Lemma D' (su3.md, Section 4); not used by the proof.

1. For W in SU(3), max over X in su(3), -Tr X^2 = 1, of Re Tr(W X^2) + (1/2) Re Tr W is at most 1/2,
   with equality at W = I. The maximum over X is the largest eigenvalue of an 8 x 8 matrix; for
   diagonal W it equals max(A, B) with A = (sqrt(e1^2 - 3 e2) - e1)/3, B = max_{j<k} -(c_j + c_k)/2,
   c_k = cos(theta_k), e1, e2 the elementary symmetric functions of the c_k.
2. Galerkin upper bounds for the gap of nu_b at r = 1.14 (polynomials of degree two, 10^5 Haar points),
   against gap2(r) = 283/128 - r/2 of Lemma D' and gap(r) = (sqrt(7/6) - g r/2)^2 + 3/2 - r of Lemma D."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from lip import haar, features
rng=np.random.default_rng(7)
# orthonormal basis of su(3) for <X,Y> = -Tr(XY)
E=[]
for i in range(3):
    for j in range(i+1,3):
        M=np.zeros((3,3),complex); M[i,j]=1; M[j,i]=-1; E.append(M/np.sqrt(2))
        M=np.zeros((3,3),complex); M[i,j]=1j; M[j,i]=1j; E.append(M/np.sqrt(2))
E.append(1j*np.diag([1,-1,0])/np.sqrt(2)); E.append(1j*np.diag([1,1,-2])/np.sqrt(6))
for X in E: assert abs(-np.trace(X@X)-1)<1e-12
def lhs(W):
    A=np.array([[np.trace(W@(Ea@Eb+Eb@Ea)/2).real for Eb in E] for Ea in E])
    return np.linalg.eigvalsh(A).max()+0.5*np.trace(W).real
def closed(th):
    c=np.cos(th); e1=c.sum(); e2=c[0]*c[1]+c[0]*c[2]+c[1]*c[2]
    A=(np.sqrt(max(e1*e1-3*e2,0))-e1)/3; B=max(-(c[j]+c[k])/2 for j in range(3) for k in range(j+1,3))
    return max(A,B)+e1/2
worst=-9; dev=0
for t1 in np.linspace(-np.pi,np.pi,121):
    for t2 in np.linspace(-np.pi,np.pi,121):
        th=np.array([t1,t2,-t1-t2]); v=lhs(np.diag(np.exp(1j*th)))
        dev=max(dev,abs(v-closed(th))); worst=max(worst,v)
for _ in range(2000):
    Z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)); Qm,Rm=np.linalg.qr(Z); Qm=Qm@np.diag(np.diag(Rm)/abs(np.diag(Rm)))
    Qm/=np.linalg.det(Qm)**(1/3); worst=max(worst,lhs(Qm))
print(f"1. max of Re Tr(W X^2) + Re Tr(W)/2 over the torus grid and 2000 random W: {worst:.12f} (bound 1/2)")
print(f"   largest deviation from the closed form max(A, B) on the torus: {dev:.2e}; value at W = I: {lhs(np.eye(3)):.12f}")
r=1.14; g=np.sqrt(2.06548)
gap2=283/128-r/2; gap1=(np.sqrt(7/6)-g*r/2)**2+1.5-r
U0=haar(100000)
def galerkin_gap(b,D=2):
    X,F,Gr=features(U0,D)
    w=np.exp(np.einsum('nij,ji->n',U0,b).real); w/=w.sum()
    mF=w@F; C=(F-mF).T@((F-mF)*w[:,None]); Gw=(Gr*np.sqrt(w)[None,:,None]).reshape(-1,Gr.shape[2]); A=Gw.T@Gw
    e,V=np.linalg.eigh(C); keep=e>1e-8*e.max(); P=V[:,keep]/np.sqrt(e[keep])
    return np.linalg.eigvalsh(P.T@A@P).min()
for name,b in [('b = 0',0*np.eye(3)),('b = rI',r*np.eye(3)),('b = r diag(1,1,0)',np.diag([r,r,0]))]:
    print(f"2. r = {r}, {name}: Galerkin gap {galerkin_gap(b):.4f}   (Lemma D': {gap2:.4f}; Lemma D: {gap1:.4f})")
