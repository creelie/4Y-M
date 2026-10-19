import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from fast import *
import glob
lat=Lat(4,'forest'); f=Fast(lat)
best=None
for fn in sorted(glob.glob('fast_4_*.pt')):
    Q=torch.load(fn)
    for s in (1,-1):
        w,v,_=f.top(Q,s)
        if best is None or w>best[0]: best=(w,fn,s,Q,v)
w,fn,s,Q,v=best
print('best',w,fn,s)
P=lat.P; H_=lambda A: A.conj().transpose(-1,-2)
W=Q[P[:,0]]@Q[P[:,1]]@H_(Q[P[:,2]])@H_(Q[P[:,3]])
tr=torch.diagonal(W,dim1=-2,dim2=-1).sum(-1).real/2
import numpy as np
tr=tr.numpy()
print('plaquette half-traces: min %.6f max %.6f'%(tr.min(),tr.max()))
print('distinct (rounded 4):',sorted(set(np.round(tr,4)))[:20])
import itertools, collections
sites=list(itertools.product(range(4),repeat=4))
SIGt=SIG
cnt=collections.Counter(); planes=collections.defaultdict(list)
k=0
for x in sites:
    for mu in range(4):
        for nu in range(mu+1,4):
            xm=list(x);xm[mu]=(xm[mu]+1)%4; xn=list(x);xn[nu]=(xn[nu]+1)%4
            L=[(tuple(x),mu),(tuple(xm),nu),(tuple(xn),mu),(tuple(x),nu)]
            sg=lambda y: SIGt[(y[0]%2)+2*(y[1]%2)+4*(y[2]%2)+8*(y[3]%2)]
            i0=[i for i,(y,m) in enumerate(L) if sg(y)==m][0]
            cnt[(mu,nu,i0,round(float(tr[k]),3))]+=1
            k+=1
for key in sorted(cnt): print(key,cnt[key])
# eigenvector weight per link direction relative to F
vv=v.reshape(-1,3); wt=(vv**2).sum(1)
print('eigvec weight per direction:',[float(wt[[i*4+mu for i in range(len(sites))]].sum()) for mu in range(4)])
