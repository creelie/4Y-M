"""Monte Carlo comparison of the exact tilted moments of `tilted_moments.py` (not used by the proof).

For random W in SU(3) and random unit h it compares, at r = 0.98,
  E_H[ |grad phi_h|^2 e^{r Re Tr U} ]            and
  E_H[ (G_{W,h} - (4/9) Re Tr(W h*))^2 e^{r Re Tr U} ]
computed exactly (Schur-Weyl, Weyl formula, Bessel sums) with a Haar Monte Carlo estimate.

Usage: python3 tilted_mc.py [samples] [seed]"""
import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from tilted_moments import *
import numpy as np
n=int(sys.argv[1]) if len(sys.argv)>1 else 200000
rng=np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 1)
r=0.98
def haar(m):
    Z=(rng.standard_normal((m,3,3))+1j*rng.standard_normal((m,3,3)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z)
    d=np.diagonal(R,axis1=1,axis2=2); Q=Q*(d/abs(d))[:,None,:]
    return Q/np.linalg.det(Q)[:,None,None]**(1/3)
U=haar(n); wt=np.exp(r*np.trace(U,axis1=1,axis2=2).real)
def Gmc(W,h):
    WU=W@U; hU=h@U
    return (0.5*np.trace(W@h.conj().T).real-0.5*np.trace(WU@hU,axis1=1,axis2=2).real
            -np.trace(WU,axis1=1,axis2=2).imag*np.trace(hU,axis1=1,axis2=2).imag/3)
def todict(M): return {(i,j):gq(Fr(M[i,j].real).limit_denominator(10**6),Fr(M[i,j].imag).limit_denominator(10**6)) for i in range(3) for j in range(3) if M[i,j]!=0}
def rat(M):
    d=todict(M); out=np.zeros((3,3),complex)
    for (i,j),v in d.items(): out[i,j]=complex(float(v[0]),float(v[1]))
    return out
ra=flint.arb(str(r)); worst=0
for trial in range(4):
    W=haar(1)[0]; h=rng.standard_normal((3,3))+1j*rng.standard_normal((3,3)); h/=np.linalg.norm(h)
    W=rat(W); h=rat(h)   # rational approximations; the comparison uses the same matrices on both sides
    Wd,hd=todict(W),todict(h)
    for name,poly,mc in (("gradient moment",G_poly(hd,hd),Gmc(h,h)),
                         ("tilted variance",pmul(g_poly(Wd,hd),g_poly(Wd,hd)),(Gmc(W,h)-4/9*np.trace(W@h.conj().T).real)**2)):
        ex=float(ev(expect(poly),ra)[0].mid()); est=(mc*wt).mean(); se=(mc*wt).std()/np.sqrt(n)
        worst=max(worst,abs(ex-est)/se)
        print(f"trial {trial} {name}: exact {ex:.6f}  Monte Carlo {est:.6f} +/- {se:.6f}")
print(f"largest deviation {worst:.2f} standard errors")
