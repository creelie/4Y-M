"""Fast polynomial ring in the 18 entries of U, Ubar up to degree DEG, compatible with series.Ring ordering."""
import sys, os, numpy as np, itertools, math, time
EX=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','exact')
sys.path.insert(0,EX)
import series
from series import NV
class FRing(series.Ring):
    def __init__(s,deg):
        s.deg=deg; mons=[]
        for d in range(deg+1):
            for c in itertools.combinations_with_replacement(range(NV),d):
                e=[0]*NV
                for v in c: e[v]+=1
                mons.append(tuple(e))
        s.mons=mons; s.N=len(mons); s.idx={m:i for i,m in enumerate(mons)}
        s.degs=np.array([sum(m) for m in mons])
        start=[int(np.searchsorted(s.degs,d)) for d in range(deg+2)]
        I=[];J=[];K=[]
        idx=s.idx
        for i,a in enumerate(mons):
            da=s.degs[i]; jmax=start[deg-da+1]
            for j in range(jmax):
                b=mons[j]
                I.append(i);J.append(j);K.append(idx[tuple(x+y for x,y in zip(a,b))])
        s.I=np.array(I,np.int32);s.J=np.array(J,np.int32);s.K=np.array(K,np.int32)
        s.fact=np.array([float(np.prod([math.factorial(x) for x in m])) for m in mons])
if __name__=='__main__':
    t=time.time(); R=FRing(int(sys.argv[1])); print(R.N,len(R.I),time.time()-t)
from cmul import cmul as _cmul
def _fmul(s,a,b): return _cmul(s,a,b)
FRing.mul=_fmul   # series.Ring.mul replaced by the C kernel (same products, same order of summation per coefficient)
