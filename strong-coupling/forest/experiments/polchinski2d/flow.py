"""Falsifier for direction #5 (Polchinski-flow / multiscale Bakry-Emery) on SU(2)=S^3.

2D lattice YM in axial gauge = independent plaquettes with density e^{beta cos th}
on unit S^3 (th = distance to identity). Heuristic proxy for the BBD criterion:
  V_t = -log P_t e^{-V},  e^{-V} = e^{beta cos th},  P_t heat semigroup on S^3,
  Lam(t) = Ric(=2) + inf_th min(V_t'', V_t' cot th)   (radial / tangential Hessian),
  B = 1 / int_0^inf exp(-int_0^t Lam ds) dt   (multiscale curvature rate).
Compared with the exact Poincare gap of e^{-V} (radial ODE, angular momentum l=0,1).
If B/gap -> 0 exponentially in beta, the antipodal obstruction survives the flow.
"""
import json, numpy as np
from scipy.special import ive
from scipy.linalg import eigh_tridiagonal

TH = np.linspace(1e-4, np.pi-1e-4, 4001)
def Vt(beta, t, J=400):
    j = np.arange(J)
    a = (j+1)*ive(j+1, beta)*np.exp(-t*j*(j+2))            # scaled by e^{-beta}
    chi = np.sin(np.outer(j+1, TH))/np.sin(TH)
    return -np.log(np.maximum(a@chi, 1e-300))
def Lam(beta, t):
    V = Vt(beta, t); h = TH[1]-TH[0]
    d1 = np.gradient(V, h); d2 = np.gradient(d1, h)
    s = slice(5, -5)
    return 2 + min(d2[s].min(), (d1*np.cos(TH)/np.sin(TH))[s].min())
def gap(beta, N=3000):
    # weighted operator on th in (0,pi), weight w = sin^2 th e^{beta cos th}; symmetric FD
    x = np.linspace(0, np.pi, N+2); h = x[1]-x[0]; xm = (x[:-1]+x[1:])/2
    w = lambda y: np.sin(y)**2*np.exp(beta*(np.cos(y)-1))
    best = np.inf
    for l in (0, 1):
        wi = w(x[1:-1]); wm = w(xm)
        diag = (wm[:-1]+wm[1:])/h**2/wi + l*(l+1)/np.sin(x[1:-1])**2
        off = -wm[1:-1]/h**2/np.sqrt(wi[:-1]*wi[1:])
        ev = eigh_tridiagonal(diag, off, select='i', select_range=(0, 1))[0]
        if l == 0 and False: pass
        best = min(best, ev[1] if l == 0 else ev[0])   # l=0 ground state is constant
    return best
out = []
ts = np.concatenate([[0], np.geomspace(1e-5, 20, 400)])
for beta in [0.5, 1, 2, 4, 8, 16, 32]:
    L = np.array([Lam(beta, t) for t in ts])
    cum = np.concatenate([[0], np.cumsum((L[1:]+L[:-1])/2*np.diff(ts))])
    integ = np.trapz(np.exp(-cum), ts) + np.exp(-cum[-1])/max(L[-1], 1e-9)
    g = gap(beta)
    r = dict(beta=beta, Lam0=float(L[0]), tstar=float(ts[np.argmax(L > 0)]) if (L > 0).any() else None,
             minLam=float(L.min()), B=float(1/integ), gap=float(g), ratio=float(1/integ/g))
    out.append(r); print(json.dumps(r), flush=True)
