"""Adaptive dyadic cover of the representatives b = diag(a1, a2, x3 + i y3), a1 >= a2 >= |z3|, y3 >= 0,
b in r conv T (su3.md, Lemma K), r = 45/32, by boxes on which the bound 2 r B_Q of Section 7 is below 1.

Usage: python3 cover.py [--r 45/32] [--L0 3] [--LMAX 6] [--out DIR]   (boxes of side r/2^L0 to start; four
worker processes). Writes cover_log.txt and cover_boxes.txt (every certified box with its j and bound) in
DIR (default: this directory); the corner data are cached in cover_cache.pkl there (not tracked)."""
import numpy as np, itertools, pickle, os, sys, time
from multiprocessing import Pool
from scipy.optimize import linprog
from corner import corner_data
from box import box_bound
import argparse
from fractions import Fraction
_ap=argparse.ArgumentParser(); _ap.add_argument('--r',default='45/32'); _ap.add_argument('--L0',type=int,default=3)
_ap.add_argument('--LMAX',type=int,default=6); _ap.add_argument('--out',default=os.path.dirname(os.path.abspath(__file__)))
_ap.add_argument('--table',action='store_true')
ARGS=_ap.parse_args() if __name__=='__main__' else _ap.parse_args([])
R=float(Fraction(ARGS.r))
OUT=os.path.abspath(ARGS.out); os.makedirs(OUT,exist_ok=True)
CACHE=os.path.join(OUT,'cover_cache.pkl')
LOG=os.path.join(OUT,'cover_log.txt')
def log(*a):
    s=' '.join(str(x) for x in a); print(s,flush=True)
    with open(LOG,'a') as f: f.write(s+'\n')
# ---- torus support function and separation ----
def torus_grid(n):
    t=np.linspace(0,2*np.pi,n,endpoint=False); T1,T2=np.meshgrid(t,t); T1=T1.ravel(); T2=T2.ravel()
    return np.stack([np.exp(1j*T1),np.exp(1j*T2),np.exp(-1j*(T1+T2))],1)
TGL=torus_grid(72); AEQ=np.vstack([TGL.real.T,TGL.imag.T])
NF=512; TGF=torus_grid(NF)
def hT_upper(c):
    """upper bound of max_{t in T} Re <c, t> = Re sum conj(c_j) t_j: grid maximum plus (1/2) H delta^2 / 2"""
    f=(np.conj(c)[None,:]*TGF).sum(1).real.max()
    H=abs(c[0])+abs(c[1])+2*abs(c[2]); d=2*np.pi/NF
    return f+0.5*H*d*d/2+1e-12*(1+np.abs(c).sum())
def separating(b0):
    z=np.array([b0[0],b0[1],complex(b0[2],b0[3])]); beq=np.concatenate([z.real,z.imag])
    res=linprog(np.ones(len(TGL)),A_eq=AEQ,b_eq=beq,bounds=(0,None),method='highs')
    if res.status!=0: return None
    y=res.eqlin.marginals; return y[:3]+1j*y[3:]
def separated(lo,hi):
    """True if the box misses r conv T (certified by a linear functional)"""
    c=separating((lo+hi)/2)
    if c is None: return False
    coef=np.array([c[0].real,c[1].real,c[2].real,c[2].imag])
    mn=np.where(coef>0,coef*lo,coef*hi).sum()
    return mn>R*hT_upper(c)*(1+1e-12)+1e-12
def misses_order(lo,hi):
    """True if the box misses a1 >= a2 >= |x3 + i y3|"""
    dx=max(lo[2],0,-hi[2]); dy=max(lo[3],0,-hi[3]); mz=np.hypot(dx,dy)
    return min(hi[1],hi[0])<max(lo[1],mz)
def corners_of(lo,hi): return [tuple(np.where(np.array(c)==1,hi,lo)) for c in itertools.product((0,1),repeat=4)]
def evalbox(args):
    lo,hi,cds=args
    try:
        b,d=box_bound(lo,hi,cds,R)
    except Exception as e:
        return lo,hi,9.0,dict(err=repr(e))
    return lo,hi,b,d
def table(passed):
    with open(os.path.join(OUT,'cover_boxes.txt'),'w') as f:
        f.write('# certified boxes: lower corner, upper corner in (a1, a2, x3, y3), j (s = 3j/8192), bound on 2 r sup_h L_h\n')
        for L,lo,hi,b,j in sorted(passed,key=lambda x:-x[3]):
            f.write(' '.join(f'{x:.9f}' for x in lo)+'   '+' '.join(f'{x:.9f}' for x in hi)+f'   {j}   {b:.6f}\n')
def main(L0=3,LMAX=6):
    cache={}
    if os.path.exists(CACHE):
        with open(CACHE,'rb') as f: cache=pickle.load(f)
    w=R/2**L0
    queue=[]
    for i1 in range(2**L0):
        for i2 in range(2**L0):
            for i3 in range(2**(L0+1)):
                for i4 in range(2**L0):
                    lo=np.array([i1*w,i2*w,-R+i3*w,i4*w]); queue.append((L0,lo,lo+w))
    passed=[]; discarded_order=0; discarded_sep=0; worst=0
    with Pool(4) as pool:
        level=L0
        while queue:
            t0=time.time()
            keep=[]
            for L,lo,hi in queue:
                if misses_order(lo,hi): discarded_order+=1; continue
                if separated(lo,hi): discarded_sep+=1; continue
                keep.append((L,lo,hi))
            need=sorted({c for L,lo,hi in keep for c in corners_of(lo,hi)}-set(cache))
            log(f'level {level}: {len(queue)} boxes, {len(keep)} kept, {len(need)} new corners (discarded so far: order {discarded_order}, separation {discarded_sep})')
            done=0
            for c,d in zip(need,pool.imap(corner_data,need,chunksize=2)):
                cache[c]=d; done+=1
                if done%200==0:
                    log(f'  corners {done}/{len(need)} {time.time()-t0:.0f}s')
                    with open(CACHE,'wb') as f: pickle.dump(cache,f)
            with open(CACHE,'wb') as f: pickle.dump(cache,f)
            tasks=[(lo,hi,{c:cache[c] for c in corners_of(lo,hi)}) for L,lo,hi in keep]
            nxt=[]; failed=0
            for (L,_,_),(lo,hi,b,d) in zip(keep,pool.imap(evalbox,tasks,chunksize=4)):
                if b<1: passed.append((L,lo,hi,b,d.get('j'))); worst=max(worst,b)
                else:
                    failed+=1; h=(hi-lo)/2
                    for c in itertools.product((0,1),repeat=4):
                        clo=lo+np.array(c)*h; nxt.append((L+1,clo,clo+h))
            log(f'level {level}: passed {len(keep)-failed}, failed {failed}, largest passing bound {worst:.5f}, {time.time()-t0:.0f}s')
            with open(os.path.join(OUT,'cover_passed.pkl'),'wb') as f: pickle.dump(passed,f)
            queue=nxt; level+=1
            if level>LMAX and queue:
                log(f'stopped at level {level}: {len(queue)} boxes left'); return False
    table(passed)
    log(f'COVER COMPLETE: {len(passed)} boxes, largest bound {worst:.6f}; discarded order {discarded_order}, separation {discarded_sep}')
    return True
if __name__=='__main__':
    if ARGS.table:
        with open(os.path.join(OUT,'cover_passed.pkl'),'rb') as f: table(pickle.load(f))
        sys.exit()
    log(f'r = {ARGS.r}, beta_W = r/2 = {R/2}')
    main(ARGS.L0,ARGS.LMAX)
