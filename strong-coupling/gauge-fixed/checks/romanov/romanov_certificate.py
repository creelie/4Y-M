"""Check of Romanov's site-plus-plaquette certificate for the forest-gauge Hessian, h_F <= 14.420001.

The certificate (cert4d_h14.42.json, data of github.com/romanov360/yang-mills-4d-forest-certificate, MIT licence,
see NOTICE) is read as data. Everything else is computed here:
  1. the finite part in exact rational arithmetic: the sixteen vertex forms G_y are positive definite with
     G_y[o_mu, i_mu] = 0, the plaquette weights are nonnegative, the largest per-link total is 14420001/10^6;
  2. on the torus of side two with random U(2) and U(3) fields: the Hessian of sum_p Re Tr Q_p equals the sum of
     the plaquette forms f_p(Y, W) of Q_p(t) = e^{tY1} e^{tY2} W e^{-tY3} e^{-tY4}, and the corner terms of the
     plaquettes add up to the vertex forms minus their diagonals;
  3. every plaquette matrix M_p(phi) is positive definite for all phases: least eigenvalue at the centre of each
     cell of an adaptive grid, against a bound for |M_p(phi) - M_p(centre)| over the cell.
Usage: python3 romanov_certificate.py cert4d_h14.42.json [first last]   (part 3 for plaquettes first..last-1)
"""
import json, itertools, time, numpy as np
from fractions import Fraction as Fr
from scipy.linalg import expm
rng=np.random.default_rng(1)
d=4; sites=sorted(itertools.product((0,1),repeat=4))
PAT="0212331001332120"; mu_of={s:int(PAT[i]) for i,s in enumerate(sites)}
sh=lambda x,a:tuple((x[i]+(i==a))%2 for i in range(4))
edges=[(x,a) for x in sites for a in range(4)]; eidx={e:i for i,e in enumerate(edges)}
fixed={eidx[(y,a)] for (y,a) in edges if mu_of[sh(y,a)]==a}
PL=[]
for x in sites:
    for m in range(4):
        for n in range(m+1,4):
            p=(eidx[(x,m)],eidx[(sh(x,m),n)],eidx[(sh(x,n),m)],eidx[(x,n)])
            corners=((x,m,n),(sh(x,m),4+m,n),(sh(sh(x,m),n),4+n,4+m),(sh(x,n),m,4+n))
            PL.append((p,corners,[i for i in range(4) if p[i] not in fixed]))
assert len(fixed)==16 and all(len(f)==3 for _,_,f in PL)
import sys
cert=json.load(open(sys.argv[1]))
G={tuple(json.loads(k)):[[Fr(v) for v in r] for r in g] for k,g in cert['G'].items()}
A=[{int(i):Fr(v) for i,v in a.items()} for a in cert['A']]
# ---- 1. finite part, exact ----
def pd(M):
    M=[r[:] for r in M]; n=len(M)
    for j in range(n):
        if M[j][j]<=0: return False
        for i in range(j+1,n):
            f=M[i][j]/M[j][j]
            for k in range(j,n): M[i][k]-=f*M[j][k]
    return True
for y in sites:
    g=G[y]; assert all(g[i][j]==g[j][i] for i in range(8) for j in range(8)); assert all(g[m][4+m]==0 for m in range(4)); assert pd(g)
assert all(v>=0 for a in A for v in a.values()) and all(set(a)==set(f) for a,(_,_,f) in zip(A,PL))
tot={e:Fr(0) for e in range(64) if e not in fixed}
for (x,m) in edges:
    e=eidx[(x,m)]
    if e not in fixed: tot[e]+=G[x][m][m]+G[sh(x,m)][4+m][4+m]
for (p,_,f),a in zip(PL,A):
    for i in f: tot[p[i]]+=a[i]
h=max(tot.values()); print('finite part: 16 G_y PD, weights>=0, h =',h,float(h)); assert h==Fr(cert['h'])
# ---- 2. the plaquette form, derived here: Q_p(t)=e^{tY1}e^{tY2} W e^{-tY3} e^{-tY4} ----
def f2(Y,W):
    # t^2 coefficient times 2 of the product of the four exponentials around W
    fac=[(Y[0],1),(Y[1],1),(None,0),(Y[2],-1),(Y[3],-1)]
    I=np.eye(len(W)); ser=[I,0*I,0*I]
    for Yk,s in fac:
        if Yk is None: ser=[S@W for S in ser]; continue
        T=[I,s*Yk,Yk@Yk/2]
        ser=[sum(ser[i]@T[k-i] for i in range(k+1)) for k in range(3)]
    return 2*np.trace(ser[2]).real
ip=lambda X,Y:-np.trace(X@Y).real
def randu(N):
    Z=rng.normal(size=(N,N))+1j*rng.normal(size=(N,N)); q,r=np.linalg.qr(Z); return q*(np.diag(r)/abs(np.diag(r)))
def randalg(N):
    Z=rng.normal(size=(N,N))+1j*rng.normal(size=(N,N)); return (Z-Z.conj().T)/2
s_=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.array([[1,0],[0,-1]],complex)]
BAS=[1j*np.eye(2)/np.sqrt(2)]+[1j*m/np.sqrt(2) for m in s_]
def qmat(form):
    vec=lambda v:[sum(v[4*i+a]*BAS[a] for a in range(4)) for i in range(4)]
    E=np.eye(16); dd=np.array([form(vec(E[i])) for i in range(16)]); M=np.diag(dd)
    for i in range(16):
        for j in range(i+1,16): M[i,j]=M[j,i]=(form(vec(E[i]+E[j]))-dd[i]-dd[j])/2
    return M
D_=lambda a,b:np.diag([a,b]).astype(complex)
FB=[qmat(lambda Y,W=W:f2(Y,W)) for W in (D_(1,0),D_(1j,0),D_(0,1),D_(0,1j))]
c23=lambda W:qmat(lambda Y:ip(Y[1],W@Y[2]@W.conj().T))
C0,Cpi,Ch=c23(D_(1,1)),c23(D_(1,-1)),c23(D_(1j,1))
CB=[(C0+Cpi)/2,(C0-Cpi)/2,Ch-(C0+Cpi)/2]
# check decompositions at random phases
for _ in range(5):
    a,b=rng.uniform(0,2*np.pi,2); W=D_(np.exp(1j*a),np.exp(1j*b))
    assert np.abs(qmat(lambda Y:f2(Y,W))-(np.cos(a)*FB[0]+np.sin(a)*FB[1]+np.cos(b)*FB[2]+np.sin(b)*FB[3])).max()<1e-12
    assert np.abs(c23(W)-(CB[0]+np.cos(a-b)*CB[1]+np.sin(a-b)*CB[2])).max()<1e-12
def cm(i,j):
    M=np.zeros((16,16)); M[4*i:4*i+4,4*j:4*j+4]=0.5*np.eye(4); return M+M.T
C14,C12,C34=cm(0,3),cm(0,1),cm(2,3)
nrm=lambda M:np.linalg.norm(M,2)
def run(j,m=64,maxdepth=12):
    p,c,f=PL[j]; k=[2*float(G[s][r][q]) for s,r,q in c]
    idx=[4*i+q for i in f for q in range(4)]; R=lambda M:M[np.ix_(idx,idx)]
    D=np.zeros((16,16))
    for i in f: D[4*i:4*i+4,4*i:4*i+4]=float(A[j][i])*np.eye(4)
    K0=R(D-k[0]*C14-k[1]*C12-k[3]*C34-k[2]*CB[0]); F=[R(x) for x in FB]; C1,C2=R(CB[1]),R(CB[2])
    nF=[nrm(x) for x in F]; nC=[nrm(C1),nrm(C2)]; kk=abs(k[2])
    cells=np.array([(a,b) for a in range(m) for b in range(m)],float); n=m; total=0; lmin=np.inf
    for depth in range(maxdepth+1):
        r=np.pi/n; lam_all=[]; bad_all=[]
        for c0 in range(0,len(cells),100000):
            cc=cells[c0:c0+100000]
            ph=(2*cc+1)*np.pi/n; p1,p2=ph[:,0],ph[:,1]; dl=p1-p2
            M=(K0[None]-np.cos(p1)[:,None,None]*F[0]-np.sin(p1)[:,None,None]*F[1]-np.cos(p2)[:,None,None]*F[2]-np.sin(p2)[:,None,None]*F[3]
               -k[2]*(np.cos(dl)[:,None,None]*C1+np.sin(dl)[:,None,None]*C2))
            lam=np.linalg.eigvalsh(M)[:,0]
            # |cos(c+t)-cos c|<=|t||sin c|+t^2/2 etc., |t|<=r for phi1, phi2 and |t|<=2r for phi1-phi2
            q=lambda c,rr:(rr*np.abs(np.sin(c))+rr*rr/2, rr*np.abs(np.cos(c))+rr*rr/2)
            a1,b1=q(p1,r); a2,b2=q(p2,r); a3,b3=q(dl,2*r)
            bound=a1*nF[0]+b1*nF[1]+a2*nF[2]+b2*nF[3]+kk*(a3*nC[0]+b3*nC[1])
            lam_all.append(lam); bad_all.append(lam<=bound+1e-9)
        lam=np.concatenate(lam_all); bad=np.concatenate(bad_all); lmin=min(lmin,lam.min()); total+=len(cells)
        if lam.min()<=0: return False,depth,total,lmin
        if not bad.any(): return True,depth,total,lmin
        cells=np.concatenate([2*cells[bad]+o for o in ((0,0),(0,1),(1,0),(1,1))]); n*=2
    return False,depth,total,lmin

if __name__=='__main__':
    if len(sys.argv)<=2:
        # global identity test on the L=2 torus with random U(N) field: Hess(X,X) = sum_p f2(Y_p,W_p); corner sums
        for N in (2,3):
            Q=[randu(N) for _ in edges]; X=[randalg(N) for _ in edges]
            S=lambda t:sum(np.trace(expm(t*X[p[0]])@Q[p[0]]@expm(t*X[p[1]])@Q[p[1]]@(expm(t*X[p[2]])@Q[p[2]]).conj().T@(expm(t*X[p[3]])@Q[p[3]]).conj().T).real for p,_,_ in PL)
            eps=1e-3; fd=(S(eps)-2*S(0)+S(-eps))/eps**2
            tot2=0; corner=0
            for p,c,_ in PL:
                Q1,Q2,Q3,Q4=(Q[e] for e in p); W=Q1@Q2@Q3.conj().T@Q4.conj().T
                Y=[X[p[0]],Q1@X[p[1]]@Q1.conj().T,Q4@X[p[2]]@Q4.conj().T,X[p[3]]]
                tot2+=f2(Y,W)
            print(f'U({N}): finite-difference Hessian {fd:.6f}, sum of plaquette forms {tot2:.6f}')
            # end vectors: o at tail = X, i at head = Ad(Q^*)X
            end={}
            for (x,m) in edges:
                e=eidx[(x,m)]; end[(x,m)]=X[e]; end[(sh(x,m),4+m)]=Q[e].conj().T@X[e]@Q[e]
            lhs=0
            for (p,c,_) in PL:
                Q1,Q2,Q3,Q4=(Q[e] for e in p); W=Q1@Q2@Q3.conj().T@Q4.conj().T
                Y=[X[p[0]],Q1@X[p[1]]@Q1.conj().T,Q4@X[p[2]]@Q4.conj().T,X[p[3]]]
                vals=[ip(Y[0],Y[3]),ip(Y[0],Y[1]),ip(Y[1],W@Y[2]@W.conj().T),ip(Y[2],Y[3])]
                for (s,r,q),v in zip(c,vals): lhs+=2*float(G[s][r][q])*v
            rhs=sum(float(G[y][r][q])*ip(end[(y,r)],end[(y,q)]) for y in sites for r in range(8) for q in range(8) if r!=q)
            print(f'U({N}): sum of corner terms {lhs:.10f}, site forms minus diagonals {rhs:.10f}')
    a,b=(int(sys.argv[2]),int(sys.argv[3])) if len(sys.argv)>3 else (0,96); t0=time.time(); allok=True
    for j in range(a,b):
        ok,dep,tot,lm=run(j); allok&=ok
        print(f'plaquette {j:2d}: {"positive definite" if ok else "FAILED"}, depth {dep}, cells {tot}, least eigenvalue at centres {lm:.6f}, {time.time()-t0:.0f}s',flush=True)
    print(f'plaquettes {a}..{b-1} positive definite' if allok else 'FAILURE')
