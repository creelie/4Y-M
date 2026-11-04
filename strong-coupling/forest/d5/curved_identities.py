"""Exact symbolic verification of the plaquette identities used in
Proposition prop:curvedsharp (SU(2), N=2).

Background links A = ca + i sigma_1 sa (direction 0), B = cb + i sigma_2 sb
(direction 1), with ca^2+sa^2 = cb^2+sb^2 = 1.  Test field values on the four
links of a (0,1)-plaquette: f X_1, gp X_2, fp X_1 (reversed), g X_2 (reversed),
X_j = i sigma_j/sqrt2.  Every identity is checked as a polynomial identity in
(ca,sa,cb,sb,f,fp,g,gp) modulo ca^2+sa^2-1, cb^2+sb^2-1 (exact arithmetic).
"""
import sympy as sp
ca,sa,cb,sb,f,fp,g,gp,t=sp.symbols('ca sa cb sb f fp g gp t',real=True)
I=sp.I; r2=sp.sqrt(2)
s1=sp.Matrix([[0,1],[1,0]]); s2=sp.Matrix([[0,-I],[I,0]]); s3=sp.Matrix([[1,0],[0,-1]])
X1,X2,X3=I*s1/r2,I*s2/r2,I*s3/r2
Id=sp.eye(2)
A=ca*Id+I*s1*sa; B=cb*Id+I*s2*sb
def adj(M): return M.T.conjugate()
def ip(X,Y): return sp.expand(-(X*Y).trace())
REL=[ca**2+sa**2-1, cb**2+sb**2-1]
def red(e):
    e=sp.expand(e)
    q,r=sp.reduced(e,REL,ca,cb,sa,sb,f,fp,g,gp)
    return sp.expand(r)
def zero(e,name):
    r=red(e)
    print(f"{name}: {'OK' if r==0 else 'FAIL '+str(r)}")
    assert r==0
def Ad(h,X): return h*X*adj(h)
def coords(X): return [sp.expand(ip(X,Xj)) for Xj in (X1,X2,X3)]
s2a,c2a=2*sa*ca,ca**2-sa**2
s2b,c2b=2*sb*cb,cb**2-sb**2
# rotation facts
for (lhs,rhs,name) in [
    (Ad(A,X2), c2a*X2-s2a*X3, 'Ad(A)X2'),
    (Ad(A,X3), c2a*X3+s2a*X2, 'Ad(A)X3'),
    (Ad(B,X1), c2b*X1+s2b*X3, 'Ad(B)X1'),
    (Ad(B,X3), c2b*X3-s2b*X1, 'Ad(B)X3'),
    (Ad(A,X1), X1, 'Ad(A)X1'), (Ad(B,X2), X2, 'Ad(B)X2')]:
    for e in (lhs-rhs): zero(e,name)
# Y_i of the (0,1)-plaquette, left-trivialised
Y1=f*Ad(B*A*adj(B),X1); Y2=gp*Ad(B*A,X2); Y3=-fp*Ad(B,X1); Y4=-g*X2
Z=Y1+Y2+Y3+Y4
print('coordinates of Y1..Y4 in X1,X2,X3:')
for Y in (Y1,Y2,Y3,Y4): print('  ',[red(c) for c in coords(Y)])
Z2=ip(Z,Z)
S=sa**2*sb**2
zero(Z2-(f**2+fp**2+g**2+gp**2-2*c2b*f*fp-2*c2a*g*gp+2*s2a*s2b*f*g),'|Z|^2 formula')
# holonomy and its su(2) part
c=A*B*adj(A)*adj(B)
zero(sp.expand(c.trace()/2)-(1-2*S),'Re part of c = 1-2S')
Chat=c-sp.expand(c.trace()/2)*Id
print('coordinates of Chat:',[sp.factor(red(x)) for x in coords(Chat)])
zero(2*S - (1-sp.expand(c.trace()/2)), 'S')
# operator norm of I-c: c has eigenvalues e^{+-i gamma}, cos gamma = 1-2S, |I-c|^2 = 2-2cos gamma = 4S
# Hessian from the exact second derivative of s(V_p(t)), via the product of exponentials
def ex(X,u):
    # exp(u X) for X in su(2) with -Tr X^2 = 1: X^2 = -1/2, exp(uX)=cos(u/sqrt2)+sqrt2 X sin(u/sqrt2)
    return sp.cos(u/r2)*Id+r2*X*sp.sin(u/r2)
Vp=A*ex(X1,t*f)*B*ex(X2,t*gp)*adj(A*ex(X1,t*fp))*adj(B*ex(X2,t*g))
s=1-(Vp.trace())/2
H=sp.diff(s,t,2).subs(t,0)
H=sp.expand(sp.simplify(H))
F,G=f+fp,g+gp
zero(sp.re(H)-(sp.Rational(1,2)*Z2-S*(F**2+G**2)+s2a*s2b/2*(F*G-2*f*g)),'Hess s(V_p) formula')
# decomposition H = (1/2 - S)|Z|^2 + 1/2 <Chat, sum_{i<i'} [Y_i,Y_i']>
Ys=[Y1,Y2,Y3,Y4]
comm=sp.zeros(2,2)
for i in range(4):
    for j in range(i+1,4):
        comm+=Ys[i]*Ys[j]-Ys[j]*Ys[i]
zero(sp.re(H)-((sp.Rational(1,2)-S)*Z2+sp.Rational(1,2)*ip(Chat,comm)),"H=(1/2-S)|Z|^2+<Chat,sum[Y,Y']>/2")
# general second-derivative formula (i) of lem:curvedloop at this plaquette
E2=Z*Z+comm
zero(sp.re(H)-sp.expand(-(c*E2).trace()/2),"Hess = -(1/N) ReTr(V_C(V0)(Z^2+sum[Y,Y']))")
# flat plaquettes in planes (0,nu), nu=2,3 and (1,nu)
Vq=A*ex(X1,t*f)*adj(A*ex(X1,t*fp))
Hq=sp.simplify(sp.diff(1-Vq.trace()/2,t,2).subs(t,0))
zero(Hq-(f-fp)**2/2,'flat (0,nu) plaquette')
Vr=B*ex(X2,t*gp)*adj(B*ex(X2,t*g))
Hr=sp.simplify(sp.diff(1-Vr.trace()/2,t,2).subs(t,0))
zero(Hr-(gp-g)**2/2,'flat (1,nu) plaquette')
print('all identities verified')
