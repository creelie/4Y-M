"""gaugegap: as gapscan.py, but gauge directions inside ker J are removed exactly (D = linearised gauge action), not by a 1e-6 cut.\nNielsen-Olesen check for (H1) sup-over-V (2D slice of 4D SU(2) Wilson action, exact sub-block).
Background: uniform abelian plaquette angle phi = 2 pi n / L^2 on an L x L torus (Landau gauge with seam).
Fluctuations U_l -> exp(sum_a e_a T_a) U_l, T_a = i sigma_a / 2... we use T_a = i sigma_a (HS-normalised up to const).
H = Hessian of S = sum_p (1 - Re tr U_p / 2) (exact per-plaquette finite differences, sparse assembly).
Constraint: block map (a) 'corner': coarse link = product of b fine links on the straight line from the block corner;
(b) 'all': all b parallel straight lines of the block (stronger). J = linearisation (right-trivialised).
Output: min eigenvalue of H on ker J, number of negative eigenvalues, with and without constraint.
Floating point: heuristic numerics only (no certification)."""
import sys, json, numpy as np, scipy.linalg as sl
L = int(sys.argv[1]); B = [int(sys.argv[2])]; NS = [int(n) for n in sys.argv[3].split(",")]
s = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]]), np.array([[1,0],[0,-1]],complex)]
T = [1j*a/2 for a in s]
def ex(v):
    th = np.linalg.norm(v)/2; A = sum(c*t for c,t in zip(v,T))
    return np.eye(2)*np.cos(th) + (A*np.sin(th)/th if th>1e-14 else A)
def d3(a): return np.diag([np.exp(1j*a), np.exp(-1j*a)])
def lid(x,y,mu): return ((x%L)*L + (y%L))*2 + mu
def background(n):
    phi = 2*np.pi*n/L**2; U = {}
    for x in range(L):
        for y in range(L):
            U[lid(x,y,1)] = d3(phi*x/2*2/2) if False else d3(phi*x)   # U2(x,y) = e^{i phi x sigma3}
            U[lid(x,y,0)] = d3(-phi*L*y) if x==L-1 else np.eye(2)
    return U, phi
def plaq(x,y): return [(lid(x,y,0),0),(lid(x+1,y,1),0),(lid(x,y+1,0),1),(lid(x,y,1),1)]
def pval(U,ls,eps):
    M = np.eye(2)
    for k,(l,dag) in enumerate(ls):
        V = ex(eps[3*k:3*k+3]) @ U[l]; M = M @ (V.conj().T if dag else V)
    return 1 - np.real(np.trace(M))/2
def hessian(U):
    n = 2*L*L*3; H = np.zeros((n,n)); h = 1e-4; g = np.zeros(n)
    for x in range(L):
        for y in range(L):
            ls = plaq(x,y); idx = [3*l+a for l,_ in ls for a in range(3)]
            f = lambda e: pval(U,ls,e); E = np.eye(12)*h; f0 = f(np.zeros(12))
            for i in range(12):
                g[idx[i]] += (f(E[i])-f(-E[i]))/(2*h)
                for j in range(i,12):
                    v = (f(E[i]+E[j])-f(E[i]-E[j])-f(-E[i]+E[j])+f(-E[i]-E[j]))/(4*h*h)
                    H[idx[i],idx[j]] += v
                    if j!=i: H[idx[j],idx[i]] += v
    return H, g
def coords(A): return np.array([np.real(np.trace(A@t.conj().T))*2 for t in T])  # inverse of sum c T
def jac(U,b,mode):
    rows = []; h = 1e-6; n = 2*L*L*3
    for bx in range(0,L,b):
        for by in range(0,L,b):
            for mu in (0,1):
                offs = [0] if mode=='corner' else range(b)
                for o in offs:
                    line = [lid(bx+(k if mu==0 else o), by+(o if mu==0 else k), mu) for k in range(b)]
                    def prod(e,line=line):
                        M = np.eye(2)
                        for k,l in enumerate(line): M = M @ (ex(e[3*k:3*k+3]) @ U[l])
                        return M
                    P0 = prod(np.zeros(3*b))
                    R = np.zeros((3,n))
                    for k,l in enumerate(line):
                        for a in range(3):
                            e = np.zeros(3*b); e[3*k+a] = h
                            R[:,3*l+a] = coords((prod(e)-prod(-e))/(2*h) @ P0.conj().T)
                    rows.append(R)
    return np.vstack(rows)
def gauge(U):
    """Linearised gauge action: site x, generator a: links leaving x get +e_a, links entering x get -Ad_U(T_a)."""
    n=2*L*L*3; cols=[]
    for x in range(L):
        for y in range(L):
            for a in range(3):
                v=np.zeros(n)
                for mu in (0,1): v[3*lid(x,y,mu)+a]+=1.0
                for l in (lid(x-1,y,0),lid(x,y-1,1)):
                    v[3*l:3*l+3]-=coords(U[l]@T[a]@U[l].conj().T)
                cols.append(v)
    return np.array(cols).T
# gfdecay: decay of the constrained (ker J, gauge-fixed) fluctuation propagator. Heuristic numerics.
b=B[0]
for n in NS:
    U,phi=background(n); H,g=hessian(U); D=gauge(U)
    for mode in ('corner','all'):
        J=jac(U,b,mode); N=sl.null_space(J); M=J@D
        u_,s_,vt=np.linalg.svd(M); tol=1e-4*s_[0]; C=vt[(s_>tol).sum():].T
        G=N.T@D@C; ug,sg,_=np.linalg.svd(G,full_matrices=True); kg=int((sg>1e-6*max(sg[0],1)).sum())
        P=N@ug[:,kg:]; Hp=P.T@H@P; e=np.linalg.eigvalsh(Hp)
        if e[0]<=1e-8: print(json.dumps({'n':n,'mode':mode,'min':float(e[0]),'skip':'not positive'}),flush=True); continue
        Gf=P@np.linalg.solve(Hp,P.T)
        l0=lid(1,1,0); prof={}
        for x in range(L):
            for y in range(L):
                for mu in (0,1):
                    d=min(abs(x-1),L-abs(x-1))+min(abs(y-1),L-abs(y-1))
                    v=float(np.linalg.norm(Gf[3*l0:3*l0+3,3*lid(x,y,mu):3*lid(x,y,mu)+3]))
                    prof[d]=max(prof.get(d,0),v)
        print(json.dumps({'n':n,'phi':phi,'mode':mode,'b':b,'min':float(e[0]),'prof':[prof[d] for d in sorted(prof)]}),flush=True)
