"""Exact re-check, in rational arithmetic, of the identity behind
Lemma lem:covgaffney:  sum_mu E^V_mu(v) = T1 + T2 - 2X,
for G = SU(2), Ad(q) of an integer quaternion q divided by |q|^2 (exactly
rational), on the window W=[1,n]^4, Lambda'=[0,n+1]^4, for random
integer backgrounds and random rational fields v on the inner links.
E^V_mu: covariant Dirichlet energies in nu != mu and Neumann energy in mu.
T1 = sum over inner plaquettes |Z_P(v)|^2, T2 = |d_V^* v|^2 on W°,
X  = sum_{x in W°, mu<nu} <(Ad(V_{x,nu}^{-1}V_{x,mu}) - Ad(V_{x+e_nu,mu}V_{x+e_mu,nu}^{-1})) v_{x,mu}, v_{x,nu}>.
The proof does not use this script."""
import itertools, random
from fractions import Fraction as Fr
random.seed(7)
def Ad(q):
    a,b,c,d=q; s=a*a+b*b+c*c+d*d
    M=[[a*a+b*b-c*c-d*d,2*(b*c-a*d),2*(b*d+a*c)],[2*(b*c+a*d),a*a-b*b+c*c-d*d,2*(c*d-a*b)],[2*(b*d-a*c),2*(c*d+a*b),a*a-b*b-c*c+d*d]]
    return [[Fr(x,s) for x in r] for r in M]
def T(M): return [list(r) for r in zip(*M)]          # inverse of a rotation
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def mv(A,v): return [sum(A[i][k]*v[k] for k in range(3)) for i in range(3)]
def sub(u,v): return [x-y for x,y in zip(u,v)]
def add(u,v): return [x+y for x,y in zip(u,v)]
def dot(u,v): return sum(x*y for x,y in zip(u,v))
Z3=[Fr(0)]*3
def check(n):
    e=[tuple(1 if i==m else 0 for i in range(4)) for m in range(4)]
    sh=lambda y,m,s=1: tuple(y[i]+s*e[m][i] for i in range(4))
    inL=lambda y: all(0<=c<=n+1 for c in y)
    A={}
    for y in itertools.product(range(n+2),repeat=4):
        for m in range(4):
            if inL(sh(y,m)): A[(y,m)]=Ad([random.randint(-3,3) or 1 for _ in range(4)])
    inner=lambda y,m: 1<=y[m]<=n and all(2<=y[r]<=n for r in range(4) if r!=m)
    Wo=lambda y: all(2<=c<=n for c in y)
    v={}
    for (y,m) in A:
        if inner(y,m): v[(y,m)]=[Fr(random.randint(-5,5),random.randint(1,4)) for _ in range(3)]
    g=lambda y,m: v.get((y,m),Z3)
    # covariant difference D_nu v_mu(x) = v_{x+e_nu,mu} - Ad(V_{x+e_mu,nu})^{-1} v_{x,mu}
    def Dv(x,m,nu):
        xm=sh(x,m); key=(xm,nu)
        R=T(A[key]) if key in A else None
        w=g(x,m)
        return sub(g(sh(x,nu),m), mv(R,w) if R is not None else Z3)
    E=Fr(0)
    for m in range(4):
        for x in itertools.product(range(n+2),repeat=4):
            for nu in range(4):
                if nu!=m:   # Dirichlet edges <-> inner plaquettes at x in the plane {mu,nu}
                    if 1<=x[m]<=n and 1<=x[nu]<=n and all(2<=x[r]<=n for r in range(4) if r not in (m,nu)):
                        d=Dv(x,m,nu); E+=dot(d,d)
                else:       # Neumann edges: both links inner
                    if inner(x,m) and inner(sh(x,m),m):
                        d=Dv(x,m,m); E+=dot(d,d)
    T1=Fr(0); X=Fr(0)
    for x in itertools.product(range(n+2),repeat=4):
        for m in range(4):
            for nu in range(m+1,4):
                if 1<=x[m]<=n and 1<=x[nu]<=n and all(2<=x[r]<=n for r in range(4) if r not in (m,nu)):
                    d=sub(Dv(x,m,nu),Dv(x,nu,m)); T1+=dot(d,d)
                if Wo(x):
                    gm=mm(T(A[(x,nu)]),A[(x,m)])
                    hm=mm(A[(sh(x,nu),m)],T(A[(sh(x,m),nu)]))
                    X+=dot(sub(mv(gm,g(x,m)),mv(hm,g(x,m))),g(x,nu))
    T2=Fr(0)
    for y in itertools.product(range(n+2),repeat=4):
        if Wo(y):
            s=Z3
            for m in range(4):
                s=add(s,sub(mv(A[(y,m)],g(y,m)),g(sh(y,m,-1),m)))
            T2+=dot(s,s)
    print(f"n={n}: E = {float(E):.6f}, T1+T2-2X = {float(T1+T2-2*X):.6f}, exact equality: {E==T1+T2-2*X}")
    return E==T1+T2-2*X
ok=all(check(n) for n in (2,3))
print("ALL EXACT" if ok else "FAILED")
