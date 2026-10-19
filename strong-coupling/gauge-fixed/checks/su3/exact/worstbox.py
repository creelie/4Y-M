"""Floating-point values of 2r times the bound of Lemma F' (exact moments, s optimised) on a 5^4 grid of the
box with the largest certified bound at r = 45/32, and which grid points lie in rK; not used by the proof."""
import numpy as np, itertools
from multiprocessing import Pool
r=45/32
def work(b4):
    from fast import evaluate,Bval
    from gauge import gauge
    a1,a2,x3,y3=b4; g=gauge([a1,a2,x3+1j*y3])
    return b4,g,2*r*Bval(evaluate(b4),283/128-r/2)[0]
if __name__=='__main__':
    lo=np.array([1.0546875,0.703125,0.3515625,0.703125]); hi=np.array([1.23046875,0.87890625,0.52734375,0.87890625])
    pts=[tuple(lo+np.array(c)*(hi-lo)/4) for c in itertools.product(range(5),repeat=4)]
    with Pool(int(__import__('sys').argv[1]) if len(__import__('sys').argv)>1 else 4) as p: res=p.map(work,pts)
    ins=sorted([x for x in res if x[1]<=r*(1+1e-9)],key=lambda x:-x[2])
    print(f'{len(ins)} of {len(res)} grid points in rK (gauge at most r)')
    for b4,g,v in ins: print('  ',np.round(b4,5),f'gauge {g:.4f}  value {v:.4f}')
    b4,g,v=max(res,key=lambda x:x[2]); print('largest value over the grid:',np.round(b4,5),f'gauge {g:.4f}  value {v:.4f}')
