"""Finite diagnostics only: a scalar fold and a convex scalar collar, not SU(2) transport."""
import json
import math
from pathlib import Path
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

started = time.monotonic()
def grid(n, lo, hi):
    x, w = leggauss(n)
    return lo+(x+1)*(hi-lo)/2, w*(hi-lo)/2

def f(q, t):
    chi = 2-4*np.sqrt(2)*np.sin((q-t)/2)
    return q+chi*(t+2.94-7*q)/8

def ft(q, t):
    chi = 2-4*np.sqrt(2)*np.sin((q-t)/2)
    return (chi+2*np.sqrt(2)*np.cos((q-t)/2)*(t+2.94-7*q))/8

def fq(q, t):
    chi = 2-4*np.sqrt(2)*np.sin((q-t)/2)
    return 1-7*chi/8-2*np.sqrt(2)*np.cos((q-t)/2)*(t+2.94-7*q)/8

def bump(x):
    x=np.asarray(x); a=1-x*x
    return np.exp(-1/np.maximum(a,1e-300))*(a>0)

qstar=brentq(lambda q: fq(q,-.5), .04,.18)
critical=float(f(qstar,-.5))
def density(n, gamma):
    q,w=grid(n,.04,.18)
    lo=np.full(n,-.54);hi=np.full(n,-.46)
    valid=(f(q,lo)<gamma)&(gamma<f(q,hi))
    for _ in range(55):
        mid=(lo+hi)/2
        below=f(q,mid)<gamma
        lo=np.where(below,mid,lo);hi=np.where(below,hi,mid)
    t=(lo+hi)/2
    weight=bump((q-.11)/.07)*bump((t+.5)/.04)*valid/ft(q,t)
    return np.array([np.dot(w,weight),np.dot(w,weight*q)])

fold=[]
for n in [128,256,512]:
    row={'n':n,'critical':critical,'qstar':qstar,'density':float(density(n,critical)[0]),'responses':[]}
    for eps in [.001,.0005,.00025]:
        minus=density(n,critical-eps);plus=density(n,critical+eps)
        row['responses'].append({'epsilon':eps,'density_slope':float((plus[0]-minus[0])/(2*eps)),
            'mean_q_slope':float((plus[1]/plus[0]-minus[1]/minus[0])/(2*eps))})
    fold.append(row)
assert abs(fq(qstar,-.5))<1e-12
assert abs(fold[-1]['density']-fold[-2]['density'])<1e-7
assert abs(fold[-1]['responses'][-1]['mean_q_slope']-fold[-2]['responses'][-1]['mean_q_slope'])<1e-4

# Independent mollified delta integral: convergence in width and grid size.
mollified=[]
for n in [256,512]:
    q,wq=grid(n,.04,.18);t,wt=grid(n,-.54,-.46)
    val=f(q[:,None],t[None,:]); weight=bump((q[:,None]-.11)/.07)*bump((t[None,:]+.5)/.04)*wq[:,None]*wt[None,:]
    for eps in [.001,.0005,.00025]:
        den=float(np.sum(weight*np.exp(-.5*((val-critical)/eps)**2))/(np.sqrt(2*np.pi)*eps))
        mollified.append({'n':n,'width':eps,'density':den,'error_to_coarea':den-fold[-1]['density']})
assert abs(mollified[-1]['error_to_coarea'])<abs(mollified[-3]['error_to_coarea'])

def collar(nq):
    x,w=grid(nq,-1,1);J=.5
    d=w*np.exp(-x*x-.1*x**4)
    T=np.exp(J*x[:,None]*x[None,:])
    E=T.copy();P=np.diag(d)
    # Separate island weights include hard domains and an internal two-peak mixture.
    a0=w*(x<0)*(np.exp(-15*(x+.8)**2)+.2*np.exp(-40*(x+.25)**2))
    a1=w*(x>.1)*np.exp(-10*(x-.7)**2)
    b0=w*(x<-.1)*np.exp(-12*(x+.6)**2)
    b1=w*(x>0)*(np.exp(-20*(x-.8)**2)+.4*np.exp(-30*(x-.3)**2))
    out=[]
    for length in range(1,10):
        K=E@P@E.T
        lk=np.log(K)
        delta=max(float(np.max(np.ptp(lk[i][None,:]-lk,axis=1))) for i in range(nq))
        Z=np.array([[a@K@b for b in [b0,b1]] for a in [a0,a1]])
        logR=float(np.log(Z[1,1])+np.log(Z[0,0])-np.log(Z[1,0])-np.log(Z[0,1]))
        assert abs(logR)<=delta+1e-12
        # Exact unbounded Gaussian collar control, separately from the bounded quartic model.
        A=2*np.eye(length)-J*(np.eye(length,k=1)+np.eye(length,k=-1))
        inv=np.linalg.inv(A);gaussian_delta=4*J*J*inv[0,-1]
        out.append({'length':length,'kernel_log_cross_ratio_bound':delta,'island_log_odds':logR,'gaussian_control_delta':float(gaussian_delta)})
        P=(P@T)*d[None,:]
    return out

collars=[{'nodes':n,'rows':collar(n)} for n in [48,96]]
# Hard label jumps limit ordinary Gauss quadrature accuracy; compare kernel separately.
assert max(abs(a['kernel_log_cross_ratio_bound']-b['kernel_log_cross_ratio_bound']) for a,b in zip(collars[0]['rows'],collars[1]['rows']))<.002

# Positive-time spectral criterion: independent finite positive measures.
spectral=[]
s=.4
for a in [.1,.05,.025]:
    energies=np.array([1.,3.,1/a]); masses=np.array([.4,.3,.3])
    C=lambda t:float(np.dot(masses,np.exp(-t*energies)))
    ratio=C(s)/C(s/2)
    moment=float(np.dot(masses,energies*np.exp(-s*energies))/C(s))
    bound=2/(math.e*s*ratio)
    assert moment<=bound
    assert C(s+.7)/C(s)>=math.exp(-bound*.7)
    spectral.append({'a':a,'physical_time_ratio':ratio,'normalised_energy':moment,'moment_bound':bound,'escaping_only_ratio':math.exp(-s/(2*a))})

report={'classification':'deterministic floating diagnostics; scalar examples only, no gauge-theory bound',
        'fold':fold,'mollified_delta':mollified,'collars':collars,'spectral':spectral,'seconds':time.monotonic()-started}
Path('result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,sort_keys=True),flush=True)
