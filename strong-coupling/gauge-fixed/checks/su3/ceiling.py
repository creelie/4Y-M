"""Galerkin values of the quantity bounded in su3.md (Sections 2-4); not used by the proof.

For b in rK and |h| = 1 the route bounds sqrt(nu_b|grad psi|^2) with -L psi = phi_h - nu_b phi_h, the
H^{-1}(nu_b) norm of phi_h - nu_b phi_h. Restricting the dual problem to polynomials of degree two gives a
lower bound for that norm (lip.Lrate, 10^5 Haar points); this script prints it and 2r times it at two
matrices b for r = 1.12, 1.25, 1.4, 1.5 (output in ceiling.txt)."""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from lip import *
import numpy as np
U0=haar(100000)
for r in (1.12,1.25,1.4,1.5):
    for name,b in (('rI',r*np.eye(3)),('r diag(1,1,-1)',r*np.diag([1,1,-1]))):
        L,V=Lrate(b,D=2,n=100000,U0=U0)
        print(f"r={r} {name}: L={L:.4f} 2rL={2*r*L:.3f}",flush=True)
