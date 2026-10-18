"""Check the explicit block formulas used in the proof of Lemma 2 against the brute-force local form."""
import numpy as np, sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from verify_lemma import plaq_hess, corners
from pairblock import pair, cartan
rng=np.random.default_rng(7)
# off-diagonal real basis of u(2) (orthonormal for -Tr(XY)) and the Cartan basis
Boff=[np.array([[0,1],[-1,0]],complex)/np.sqrt(2), np.array([[0,1j],[1j,0]])/np.sqrt(2)]
Bc=[np.array([[1j,0],[0,0]]), np.array([[0,0],[0,1j]])]
worst=0
for _ in range(300):
    al,ga=rng.uniform(0,2*np.pi,2); u,v=np.exp(1j*al),np.exp(1j*ga)
    W=np.diag([u,v]); I=np.eye(2,dtype=complex)
    Qs=[I,I,I,W.conj().T]                       # gauge: Q1=Q2=Q3=I, holonomy Q_p = W
    for B,ref in [(Boff, np.sort(np.repeat(np.linalg.eigvalsh(4*np.eye(4)-pair(u,v))/2,2))),
                  (Bc,  np.sort(np.concatenate([np.linalg.eigvalsh(2*np.eye(4)-cartan(u)),np.linalg.eigvalsh(2*np.eye(4)-cartan(v))])))]:
        Gb=2*np.eye(4*len(B))-(plaq_hess(Qs,B)+corners(Qs,B))
        worst=max(worst,np.max(abs(np.sort(np.linalg.eigvalsh(Gb))-ref)))
print("max |brute-force spectrum - block-formula spectrum| =",worst)
# m=1 closed form: eigenvalues 4c^2,4c^2,4(1+c)^2,4(1-c)^2, c=cos(theta/2)
w2=0
for th in rng.uniform(0,2*np.pi,200):
    c=np.cos(th/2); ev=np.sort(np.linalg.eigvalsh(4*np.eye(4)-pair(np.exp(1j*th),np.exp(-1j*th))))
    w2=max(w2,np.max(abs(ev-np.sort([4*c*c,4*c*c,4*(1+c)**2,4*(1-c)**2]))))
print("max deviation from closed-form m=1 eigenvalues =",w2)
# affine dependence on m: G(chi,theta) = (1+m)/2 G(0,theta) + (1-m)/2 G(pi,theta), m=cos chi
w3=0
for _ in range(200):
    chi,th=rng.uniform(0,2*np.pi,2); m=np.cos(chi)
    Gm=lambda ch: 4*np.eye(4)-pair(np.exp(1j*(ch+th)),np.exp(1j*(ch-th)))
    w3=max(w3,np.max(abs(Gm(chi)-((1+m)/2*Gm(0)+(1-m)/2*Gm(np.pi)))))
print("affine-in-m defect =",w3)
# Cartan closed form: 2|y|^2 - q = (y1+y3)^2+(y2+y4)^2+(1+cos phi)(y1+y2-y3-y4)^2
w4=0
for ph in rng.uniform(0,2*np.pi,200):
    a=np.array([1,0,1,0.]); b=np.array([0,1,0,1.]); z=np.array([1,1,-1,-1.])
    M=np.outer(a,a)+np.outer(b,b)+(1+np.cos(ph))*np.outer(z,z)
    w4=max(w4,np.max(abs(M-(2*np.eye(4)-cartan(np.exp(1j*ph))))))
print("Cartan SOS defect =",w4)
