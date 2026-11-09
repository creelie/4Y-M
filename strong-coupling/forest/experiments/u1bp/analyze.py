"""varpi per observable with jackknife (20 blocks) + effective masses from slice correlators. usage: analyze.py files..."""
import sys, json, numpy as np
NOBS=5; names=['sinB12_k0','sinB12_kx1','sinB12_kx3','cosP_all_k0','cos12_k0']
def jk(f, data, nb=20):
    m=len(data)//nb*nb; d=data[:m].reshape(nb,-1,*data.shape[1:]); full=f(data[:m])
    reps=np.array([f(np.concatenate([d[j] for j in range(nb) if j!=i])) for i in range(nb)])
    return full, np.sqrt((nb-1)*((reps-reps.mean(0))**2).sum(0))
for fn in sys.argv[1:]:
    hdr=dict(kv.split('=') for kv in open(fn).readline()[2:].split())
    L=int(hdr['L']); n=int(hdr['n']); V=L**3
    A=np.array([[float(x) for x in l.split()[1:]] for l in open(fn) if l.startswith('C')])
    if len(A)<40: continue
    num=A[:,:NOBS]; F=A[:,NOBS:3*NOBS].reshape(-1,NOBS,2); F=F[...,0]+1j*F[...,1]; sl=A[:,3*NOBS:].reshape(-1,2,L)
    def varpi(idx):
        D=np.concatenate([num, F.real, F.imag],1)
        def g(X):
            nm=X[:,:NOBS].mean(0); f=X[:,NOBS:2*NOBS]+1j*X[:,2*NOBS:]
            var=(np.abs(f-f.mean(0))**2).mean(0)/V
            return nm/(n**3*var)
        return jk(g,D)
    vp,ve=varpi(None)
    def meff(X):
        out=[]
        for ch in range(2):
            s=X[:,ch,:]-X[:,ch,:].mean(); C=np.array([(s*np.roll(s,-t,1)).mean() for t in range(4)])
            out+= [np.log(C[t]/C[t+1]) if C[t+1]>0 and C[t]>0 else np.nan for t in range(1,3)]
        return np.array(out)
    mp,me=jk(meff,sl)
    print(json.dumps(dict(L=L,beta=float(hdr['beta']),n=n,s=int(hdr['s']),R=int(hdr['R']),nconf=len(A),
        varpi={k:[float(a),float(b)] for k,a,b in zip(names,vp,ve)},
        meff={k:[float(a),float(b)] for k,a,b in zip(['cosP_t1','cosP_t2','sinB12_t1','sinB12_t2'],mp,me)})))
