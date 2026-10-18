"""Exact (sympy) verification of the pair-block steps of Lemma 2.
H(u,v) = Hq(u,v) + Hq(v,u)^T is the Hermitian matrix of [f_p + c_p] on an off-diagonal pair.
(1) with u=e^{i(chi+t)}, v=e^{i(chi-t)}, 4I - H is affine in m = cos(chi) (entries depend on chi only via cos chi
    once the common phase is removed by a diagonal unitary);  (2) at chi=0 the conjugated matrix equals the
    Gram matrix of the explicit sum of squares;  (3) the chi=pi case equals the chi=0 case at t+pi."""
import sympy as sp
I=sp.I
def Hq(u,v):
    H=sp.zeros(4,4)
    def off(i,j,c):
        H[j,i]+=c; H[i,j]+=sp.conjugate(c)
    off(0,1,-(1+u)); off(2,3,-(1+u)); off(0,2,v); off(0,3,1+v); off(1,2,v*(1+sp.conjugate(u))); off(1,3,v)
    return H-sp.re(u)*sp.eye(4)
def pair(u,v): return Hq(u,v)+Hq(v,u).T
t,chi=sp.symbols('t chi',real=True)
def M(chi_,t_):
    u=sp.exp(I*(chi_+t_)); v=sp.exp(I*(chi_-t_))
    return 4*sp.eye(4)-pair(u,v)
# (2) chi=0: conjugate by D=diag(1,e^{it/2},e^{-it},e^{-it/2})  (z = D w)
D=sp.diag(1,sp.exp(I*t/2),sp.exp(-I*t),sp.exp(-I*t/2))
G=(D.H*M(0,t)*D).applyfunc(lambda e: sp.simplify(sp.expand_complex(sp.expand(e))))
c=sp.cos(t/2)
w=sp.Matrix(sp.symbols('w1:5'))
def gram(vec): return vec*vec.T   # real coefficient vectors -> rank-one Gram
S=2*c**2*(gram(sp.Matrix([1,0,1,0]))+gram(sp.Matrix([0,1,0,1])))+(1+c)**2*gram(sp.Matrix([1,1,-1,-1]))+(1-c)**2*gram(sp.Matrix([1,-1,-1,1]))
diff=(G-S).applyfunc(lambda e: sp.simplify(sp.expand_trig(sp.expand_complex(e)).rewrite(sp.cos)))
print("(2) chi=0 conjugated matrix minus SOS Gram matrix is zero:", all(sp.simplify(x)==0 for x in diff))
# (1) affine in cos(chi): remove the common phase by z_i -> conjugation with diag(1,1,1,1)? check entries of M directly
Mc=M(chi,t).applyfunc(lambda e: sp.simplify(sp.expand_complex(sp.expand(e))))
m=sp.symbols('m',real=True)
ok=True
for e in Mc:
    e2=sp.expand(sp.expand_trig(e))
    # entries are polynomials in cos(chi), sin(chi); check they are affine in (cos chi, sin chi) jointly with sin chi coefficient zero
    p=sp.Poly(e2.subs({sp.sin(chi):sp.Symbol('S'),sp.cos(chi):sp.Symbol('C')}),sp.Symbol('S'),sp.Symbol('C'))
    for (a,b),coef in p.terms():
        if a>0 or b>1:
            if sp.simplify(coef)!=0: ok=False
print("(1) every entry of 4I-H is affine in cos(chi) with no sin(chi):", ok)
# (3)
d3=(M(sp.pi,t)-M(0,t+sp.pi)).applyfunc(lambda e: sp.simplify(sp.expand_complex(sp.expand(e))))
print("(3) M(pi,t) == M(0,t+pi):", all(x==0 for x in d3))
