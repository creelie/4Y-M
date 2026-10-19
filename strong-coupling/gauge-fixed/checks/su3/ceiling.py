"""Galerkin values of the quantity bounded in su3.md (Sections 2-4); not used by the proof.

For b in rK and |h| = 1 the route bounds sqrt(nu_b|grad psi|^2) with -L psi = phi_h - nu_b phi_h, the
H^{-1}(nu_b) norm of phi_h - nu_b phi_h. Restricting the dual problem to polynomials of degree two gives a
lower bound for that norm (lip.Lrate, 10^5 Haar points). Every matrix below lies in rK = conv(r SU(3)):
rho I with rho <= r, r diag(0,0,1) = r (I + diag(-1,-1,1))/2, (r/2) e^{i pi/3} I = r (I + omega I)/2, and
r (I + diag(e^{i pi/3}, e^{i pi/3}, e^{-2i pi/3}))/2. Output in ceiling.txt."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from lip import *
import numpy as np
U0=haar(100000)
w=np.exp(1j*np.pi/3)
def mats(r):
    return [('1.1 I',1.1*np.eye(3)),('r I',r*np.eye(3)),('r diag(0,0,1)',r*np.diag([0,0,1.])),
            ('(r/2) e^{i pi/3} I',r/2*w*np.eye(3)),('r (I + diag(w,w,w^-2))/2',r*(np.eye(3)+np.diag([w,w,w**-2]))/2)]
for r in (1.25,1.4,1.5,1.55):
    best=0
    for name,b in mats(r):
        if name=='1.1 I' and r<1.1: continue
        L,V=Lrate(b,D=2,n=100000,U0=U0); best=max(best,L)
        print(f"r={r} b={name}: L={L:.4f} 2rL={2*r*L:.3f}",flush=True)
    print(f"r={r}: largest 2rL over these matrices {2*r*best:.3f}",flush=True)
