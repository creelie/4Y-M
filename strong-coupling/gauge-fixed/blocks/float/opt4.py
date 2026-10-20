import numpy as np, sys, time
from scipy.optimize import minimize
from tri import Grid
from cert2 import K_cert
from red import fields, LO, HI
beta=float(sys.argv[1]); n=int(sys.argv[2]); seed=int(sys.argv[3])
g=Grid(12,5); rng=np.random.default_rng(seed)
def val(p):
    p=np.clip(p,LO,HI); b1,b2,b3,h=fields(p,beta)
    return K_cert(b1,b2,b3,h,beta,g,thetas=np.linspace(0.3,0.7,9))[0]
t=time.time(); P=[]
for i in range(n):
    p=LO+(HI-LO)*rng.uniform(size=6)
    if i%2==0: p[:2]=1-0.1*rng.uniform(size=2)
    P.append((val(p),p))
P.sort(key=lambda a:-a[0]); print("random best",[round(a[0],5) for a in P[:5]],"t",round(time.time()-t),flush=True)
for v0,p0 in P[:3]:
    r=minimize(lambda p:-val(p),p0,method='Nelder-Mead',options={'maxfev':300,'xatol':1e-3,'fatol':1e-6})
    p=np.clip(r.x,LO,HI); print("opt",round(-r.fun,6),np.round(p,4),"t",round(time.time()-t),flush=True)
    np.save(f"cert_{beta}_{seed}_{round(-r.fun,6)}.npy",p)
