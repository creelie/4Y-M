"""varpi per observable (jackknife, 20 blocks). usage: analyze.py files..."""
import sys, json, numpy as np
NOBS=9; names=[f'W{w}_{t}' for w in (1,2,3) for t in ('k0all','kx1','kx3')]
for fn in sys.argv[1:]:
    hdr=dict(kv.split('=') for kv in open(fn).readline()[2:].split()); L=int(hdr['L']); n=int(hdr['n']); V=L**3
    A=np.array([[float(x) for x in l.split()[1:]] for l in open(fn) if l.startswith('C')])
    if len(A)<40: continue
    nb=20; m=len(A)//nb*nb; A=A[:m]
    def g(X):
        nm=X[:,:NOBS].mean(0); f=X[:,NOBS::2]+1j*X[:,NOBS+1::2]
        return nm/(n**3*(np.abs(f-f.mean(0))**2).mean(0)/V)
    full=g(A); d=A.reshape(nb,-1,A.shape[1]); reps=np.array([g(np.concatenate([d[j] for j in range(nb) if j!=i])) for i in range(nb)])
    err=np.sqrt((nb-1)*((reps-reps.mean(0))**2).sum(0))
    print(json.dumps(dict(L=L,beta=float(hdr['beta']),n=n,s=int(hdr['s']),nconf=len(A),ellg2=4*n/float(hdr['beta']),
        varpi={k:[float(a),float(b)] for k,a,b in zip(names,full,err)})))
