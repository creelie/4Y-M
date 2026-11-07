"""shelltv: Gaussian part of Lemma D (POLYMER_ARCHITECTURE 10.4). Abelian 2D, cylinder periodic in y (period P),
x in [-b, W). Links in x<0 strip = fixed outer data. Free links: x>=0. Corner-map constraint on every block fully in
free region. Gauge fixing kappa*D Pi D^T on free non-corner sites. Output: ||Sigma_y^{-1/2} M_y||_op for the links on
column x=t*b (inner boundary), as function of t. Heuristic numerics, beta=1."""
import sys, json, numpy as np, scipy.linalg as sl
b=int(sys.argv[1]); P=int(sys.argv[2]); T=[int(t) for t in sys.argv[3].split(',')]; kappa=1.0; INNER=4*b
WP=[float(w) for w in sys.argv[4].split(',')] if len(sys.argv)>4 else [1.0]   # weak-patch plaquette weight (cos theta_p of a large-field patch)
PATCH=int(sys.argv[5]) if len(sys.argv)>5 else 2*b   # patch spans x in [X0P,X0P+PATCH), y in [0,P//4)
AT_SHELL=len(sys.argv)>6 and sys.argv[6]=='shell'   # if set, patch centred on measured column x=t*b
def run(t,wp=1.0):
    W=t*b+INNER; X0=-b
    xs=range(X0,W); lid={}
    for x in xs:
        for y in range(P):
            for mu in (0,1):
                if mu==0 and x==W-1: continue
                lid[(x,y,mu)]=len(lid)
    n=len(lid)
    def L(x,y,mu): return lid.get((x,y%P,mu))
    # curvature d: plaquettes with all links present
    rows=[]; wts=[]
    for x in xs:
        for y in range(P):
            ls=[L(x,y,0),L(x+1,y,1),L(x,y+1,0),L(x,y,1)]
            if None in ls: continue
            r=np.zeros(n); r[ls[0]]+=1; r[ls[1]]+=1; r[ls[2]]-=1; r[ls[3]]-=1; rows.append(r); x0p=(t*b-PATCH//2) if AT_SHELL else 0; wts.append(wp if (x0p<=x<x0p+PATCH and 0<=y<P//4) else 1.0)
    d=np.array(rows); H=d.T@(np.array(wts)[:,None]*d)
    fixed=[i for (x,y,mu),i in lid.items() if x<0]; free=[i for (x,y,mu),i in lid.items() if x>=0]
    # gauge: sites with x>=1 (all their links free), excluding block corners
    G=[]
    for x in range(1,W):
        for y in range(P):
            if x%b==0 and y%b==0: continue
            v=np.zeros(n)
            for l in (L(x,y,0),L(x,y,1)):
                if l is not None: v[l]+=1
            for l in (L(x-1,y,0),L(x,y-1,1)):
                if l is not None: v[l]-=1
            G.append(v)
    G=np.array(G).T; Q=H+kappa*G@G.T
    # constraint rows: corner lines of blocks with bx>=0
    J=[]
    for bx in range(0,W-b+1,b):
        for by in range(0,P,b):
            for mu in (0,1):
                r=np.zeros(n)
                for k in range(b):
                    l=L(bx+k,by,0) if mu==0 else L(bx,by+k,1)
                    if l is None: r=None;break
                    r[l]=1
                if r is not None: J.append(r[free])
    J=np.array(J); N=sl.null_space(J)            # free-link coordinates z = N w
    Qff=N.T@Q[np.ix_(free,free)]@N; Qfx=N.T@Q[np.ix_(free,fixed)]
    Ci=np.linalg.inv(Qff); Mw=-Ci@Qfx
    yl=[free.index(lid[(t*b,y,1)]) for y in range(P)]   # vertical links on inner-boundary column
    E=N[yl,:]; My=E@Mw; Sy=E@Ci@E.T
    w,V=np.linalg.eigh(Sy); k=w>1e-9*w.max(); V=V[:,k]; Sm=np.diag(w[k]**-0.5)@V.T   # restrict to support of Sy (constrained dirs removed)
    return {'wp':wp,'t':t,'n':n,'rank':int(k.sum()),'leak':float(np.linalg.norm(My-V@(V.T@My),2)),'minSy':float(w[k].min()),'minSy_all':float(w.min()),'minQ':float(np.linalg.eigvalsh(Qff).min()),'normM':float(np.linalg.norm(My,2)),'tv_op':float(np.linalg.norm(Sm@My,2))}
for wp in WP:
    for t in T:
        try: print(json.dumps(run(t,wp)),flush=True)
        except Exception as e: print(json.dumps({'wp':wp,'t':t,'err':str(e)}),flush=True)
