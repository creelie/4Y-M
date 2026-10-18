"""Symbolic/numeric check of the two explicit sum-of-squares identities in Lemma 2."""
import sympy as sp, numpy as np, sys, os
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from pairblock import pair
y1,y2,y3,y4,ph=sp.symbols('y1 y2 y3 y4 phi',real=True); c=sp.cos(ph)
f=-c*(y1+y2-y3-y4)**2; corner=2*y1*y4-2*y1*y2+2*y2*y3-2*y3*y4
G=2*(y1**2+y2**2+y3**2+y4**2)-f-corner
sos=(y1+y3)**2+(y2+y4)**2+(1+c)*(y1+y2-y3-y4)**2
print("Cartan identity:", sp.simplify(sp.expand(G-sos))==0)
rng=np.random.default_rng(5); worst=0
for _ in range(2000):
    th=rng.uniform(0,2*np.pi); e=np.exp(1j*th/2); cc=np.cos(th/2)
    w=rng.normal(size=4)+1j*rng.normal(size=4); z=np.array([w[0],e*w[1],np.conj(e)**2*w[2],np.conj(e)*w[3]])
    Gz=np.real(np.conj(z)@(4*np.eye(4)-pair(np.exp(1j*th),np.exp(-1j*th)))@z)
    S=2*cc**2*(abs(w[0]+w[2])**2+abs(w[1]+w[3])**2)+(1+cc)**2*abs(w[0]+w[1]-w[2]-w[3])**2+(1-cc)**2*abs(w[0]-w[1]-w[2]+w[3])**2
    worst=max(worst,abs(Gz-S))
print("pair identity (m=1) max defect:",worst)
