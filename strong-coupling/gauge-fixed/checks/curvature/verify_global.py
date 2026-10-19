"""Global check of the gauge-fixed certificate on a torus of side L, group SU(N) or U(N).

For a random field Q and random tangent X with X_F = 0 it evaluates
  sum_p S_p(X) + sum_y V_y(X)   and   h|X|^2 -/+ Hess S_0(X),
which must agree, and checks S_p >= 0 (as a quadratic form, on every plaquette) and the link counts
(each active link adjacent to the frozen link in 4 plaquettes, opposite it in 2).
"""
import numpy as np, itertools, sys
from fractions import Fraction as Fr
SIG=[0,0,2,3,3,1,2,1,1,2,1,3,3,2,0,0]
d=4
a,b,t=1.05*2,1.25,0.9; gam=2*t; h=4*a+2*b+2*gam
def basis(N,unitary):
    B=[]
    for i in range(N):
        for j in range(i+1,N):
            E=np.zeros((N,N),complex); E[i,j]=1; E[j,i]=-1; B.append(E/np.sqrt(2))
            E=np.zeros((N,N),complex); E[i,j]=1j; E[j,i]=1j; B.append(E/np.sqrt(2))
    for k in range(1,N):
        D=np.zeros(N,complex); D[:k]=1; D[k]=-k; D=D/np.sqrt(k*(k+1)); B.append(np.diag(1j*D))
    if unitary: B.append(1j*np.eye(N)/np.sqrt(N))
    return np.array(B)
rt=lambda A: np.real(np.trace(A))
ipr=lambda A,B: -rt(A@B)
def rand_group(N,unitary,rng):
    Z=(rng.normal(size=(N,N))+1j*rng.normal(size=(N,N)))/np.sqrt(2)
    Q,R=np.linalg.qr(Z); Q=Q@np.diag(np.diag(R)/abs(np.diag(R)))
    if not unitary: Q=Q/np.linalg.det(Q)**(1/N)
    return Q
def fp(Y,W):
    v=sum(rt(Y[i]@Y[i]@W) for i in range(4))+2*rt(Y[0]@Y[1]@W)+2*rt(Y[2]@Y[3]@W)
    v-=2*sum(rt(Y[i]@W@Y[j]) for i in (0,1) for j in (2,3))
    return v
def cp(Y,W,i0):
    c={(0,3):ipr(Y[0],Y[3]),(0,1):-ipr(Y[0],Y[1]),(2,3):-ipr(Y[2],Y[3]),(1,2):ipr(Y[1],W@Y[2]@W.conj().T)}
    return 2*sum(v for k,v in c.items() if i0 not in k)
opp={0:2,1:3,2:0,3:1}
def Sp(Y,W,i0,sgn):
    return sum((b if i==opp[i0] else a)*ipr(Y[i],Y[i]) for i in range(4) if i!=i0)-sgn*fp(Y,W)-t*cp(Y,W,i0)
def main(L,N,unitary,seed):
    rng=np.random.default_rng(seed); B=basis(N,unitary); dg=len(B)
    sites=list(itertools.product(range(L),repeat=d)); idx={x:i for i,x in enumerate(sites)}
    sh=lambda x,mu,s=1: tuple((x[k]+(s if k==mu else 0))%L for k in range(d))
    sg=lambda x: SIG[(x[0]%2)+2*(x[1]%2)+4*(x[2]%2)+8*(x[3]%2)]
    Q={(x,mu):rand_group(N,unitary,rng) for x in sites for mu in range(d)}
    X={}
    for x in sites:
        for mu in range(d):
            X[(x,mu)]=np.zeros((N,N),complex) if mu==sg(x) else np.einsum('a,aij->ij',rng.normal(size=dg),B)
    H=lambda A: A.conj().T
    norm2=sum(ipr(v,v) for v in X.values())
    tot={1:0.0,-1:0.0}; hess=0.0; minS=np.inf; counts={}
    for x in sites:
        for mu,nu in itertools.combinations(range(d),2):
            Lk=[(x,mu),(sh(x,mu),nu),(sh(x,nu),mu),(x,nu)]
            F=[i for i,(s,m) in enumerate(Lk) if sg(s)==m]; assert len(F)==1; i0=F[0]
            for i,l in enumerate(Lk):
                if i==i0: continue
                key='opp' if i==opp[i0] else 'adj'; counts.setdefault(l,{'adj':0,'opp':0})[key]+=1
            Q1,Q2,Q3,Q4=[Q[l] for l in Lk]; W=Q1@Q2@H(Q3)@H(Q4)
            Y=[X[Lk[0]],Q1@X[Lk[1]]@H(Q1),Q4@X[Lk[2]]@H(Q4),X[Lk[3]]]
            hess+=fp(Y,W)
            for s in (1,-1): tot[s]+=Sp(Y,W,i0,s)
            # quadratic-form check of S_p on this plaquette
            act=[i for i in range(4) if i!=i0]; n=3*dg; M={1:np.zeros((n,n)),-1:np.zeros((n,n))}
            for s in (1,-1):
                for r in range(n):
                    for c_ in range(r,n):
                        def ev(vec):
                            Yv=[np.zeros((N,N),complex)]*4; Yv=list(Yv)
                            for k,i in enumerate(act): Yv[i]=np.einsum('a,aij->ij',vec[k*dg:(k+1)*dg],B)
                            return Sp(Yv,W,i0,s)
                        e1=np.zeros(n); e1[r]=1; e2=np.zeros(n); e2[c_]=1
                        if r==c_: M[s][r,r]=ev(e1)
                        else: M[s][r,c_]=M[s][c_,r]=(ev(e1+e2)-ev(e1)-ev(e2))/2
                minS=min(minS,np.linalg.eigvalsh(M[s])[0])
    V=0.0
    for y in sites:
        ends=[]
        for mu in range(d):
            if mu==sg(y): continue
            ends.append((mu,1,X[(y,mu)]))
            q=Q[(sh(y,mu,-1),mu)]; ends.append((mu,-1,H(q)@X[(sh(y,mu,-1),mu)]@q))
        for (m1,s1,A) in ends: V+=gam*ipr(A,A)
        for (m1,s1,A),(m2,s2,Bm) in itertools.permutations(ends,2):
            if m1!=m2: V+=t*s1*s2*ipr(A,Bm)
    cnt_ok=all(v=={'adj':4,'opp':2} for v in counts.values()) and len(counts)==len(sites)*(d-1)
    for s in (1,-1):
        lhs=tot[s]+V; rhs=h*norm2-s*hess
        print(f'L={L} N={N} U={unitary} sign={s:+d}: lhs={lhs:.10f} rhs={rhs:.10f} diff={lhs-rhs:.2e}')
    print('  h =',h,' min eigenvalue of S_p over plaquettes:',minS,' link counts 4 adjacent / 2 opposite:',cnt_ok,flush=True)
if __name__=="__main__":
    for (L,N,U,seed) in [(4,2,False,1),(4,2,False,2),(4,3,False,3),(4,3,True,4)]:
        main(L,N,U,seed)
