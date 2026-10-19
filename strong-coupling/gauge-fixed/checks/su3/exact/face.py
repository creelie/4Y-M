import numpy as np, sys
from multiprocessing import Pool
def work(args):
    r,rho=args
    from fast import evaluate,Bval
    d=evaluate((r,rho,rho,0.0)); B=Bval(d,283/128-r/2,svals=np.linspace(0.2,0.42,111))
    return r,rho,2*r*B[0],B[1]
jobs=[(r,f*r) for r in (1.36,1.40,1.44) for f in (0.6,0.65,0.7,0.75,0.8,0.85,0.9,1.0)]
with Pool(4) as p:
    for r,rho,v,s in p.map(work,jobs,chunksize=1): print(f"r={r} rho/r={rho/r:.2f}: 2rB={v:.4f} s={s:.3f}",flush=True)
