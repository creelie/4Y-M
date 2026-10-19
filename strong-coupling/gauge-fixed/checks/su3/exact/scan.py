import numpy as np, sys, itertools, time
from multiprocessing import Pool
r=float(sys.argv[1]); n=int(sys.argv[2]) if len(sys.argv)>2 else 5
def work(b4):
    from fast import evaluate,Bval
    from gauge import gauge
    a1,a2,x3,y3=b4; g=gauge([a1,a2,x3+1j*y3])
    if g>r*(1+1e-9): return (b4,g,None,None)
    d=evaluate(b4)
    Bu=Bval(d,283/128-r/2); Bg=Bval(d,283/128-g/2)
    return (b4,g,Bu,Bg)
pts=[]
L=np.linspace(0,r,n)
for a1 in L[1:]:
    for a2 in L:
        if a2>a1+1e-12: continue
        for rho in L:
            if rho>a2+1e-12: continue
            for psi in (np.linspace(0,np.pi,5) if rho>0 else [0]):
                pts.append((a1,a2,rho*np.cos(psi),rho*np.sin(psi)))
print(len(pts),'points',flush=True)
t=time.time()
with Pool(4) as p:
    res=p.map(work,pts,chunksize=1)
inside=[x for x in res if x[2] is not None]
print(len(inside),'inside rK',time.time()-t)
inside.sort(key=lambda x:-x[2][0])
for b4,g,Bu,Bg in inside[:15]:
    print(np.round(b4,3),f"gauge {g:.3f}  2rB(uniform gap) {2*r*Bu[0]:.4f} s={Bu[1]:.3f}  2rB(gauge gap) {2*r*Bg[0]:.4f}")
inside.sort(key=lambda x:-x[3][0])
print('by gauge gap:')
for b4,g,Bu,Bg in inside[:8]:
    print(np.round(b4,3),f"gauge {g:.3f}  2rB(uniform gap) {2*r*Bu[0]:.4f}  2rB(gauge gap) {2*r*Bg[0]:.4f} s={Bg[1]:.3f}")
