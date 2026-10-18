import sympy as sp
t=sp.symbols('t',real=True)
def expmY(Y,s):  # exact for small: use series to 2nd order
    I=sp.eye(Y.shape[0]); return I+s*Y+s**2*Y*Y/2
def ip(X,Y): return sp.re(sp.expand(-(X*Y).trace()))  # <X,Y>=-Tr(XY) real for skew-Hermitian
def fc(Ys,W):
    Y1,Y2,Y3,Y4=Ys
    Qp=expmY(Y1,t)*expmY(Y2,t)*W*expmY(-Y3,t)*expmY(-Y4,t)
    f=sp.diff(sp.expand(Qp.trace()),t,2).subs(t,0)
    f=sp.re(sp.expand(f))
    c=2*ip(Y1,Y4)-2*ip(Y1,Y2)+2*ip(Y2,W*Y3*W.H)-2*ip(Y3,Y4)
    return sp.expand(f),sp.expand(c)
# Cartan block (N=1)
y=sp.symbols('y1:5',real=True); ph=sp.symbols('phi',real=True)
W=sp.Matrix([[sp.exp(sp.I*ph)]]); Ys=[sp.Matrix([[sp.I*yi]]) for yi in y]
f,c=fc(Ys,W)
F=sp.expand(2*sum(yi**2 for yi in y)-f-c)
S=(y[0]+y[2])**2+(y[1]+y[3])**2+(1+sp.cos(ph))*(y[0]+y[1]-y[2]-y[3])**2
print("Cartan identity:",sp.simplify(sp.expand_complex(F-S))==0)
# pair block N=2
x=sp.symbols('x1:5',real=True); w=sp.symbols('w1:5',real=True)
z=[x[i]+sp.I*w[i] for i in range(4)]
ch,th=sp.symbols('chi theta',real=True)
u=sp.exp(sp.I*(ch+th)); v=sp.exp(sp.I*(ch-th))
W=sp.diag(u,v); Ys=[sp.Matrix([[0,zi],[-sp.conjugate(zi),0]]) for zi in z]
f,c=fc(Ys,W)
F=sp.expand_complex(sp.expand(2*sum(2*(x[i]**2+w[i]**2) for i in range(4))-f-c))
F0=F.subs(ch,0); Fpi=F.subs(ch,0).subs(th,th+sp.pi)
print("convex combination:",sp.simplify(sp.expand_trig(F-((1+sp.cos(ch))/2*F0+(1-sp.cos(ch))/2*Fpi)))==0)
# SOS at chi=0 with substitution
a=[sp.symbols('a%d'%i,real=True) for i in range(4)]; b=[sp.symbols('b%d'%i,real=True) for i in range(4)]
W_=[a[i]+sp.I*b[i] for i in range(4)]
phase=[1,sp.exp(sp.I*th/2),sp.exp(-sp.I*th),sp.exp(-sp.I*th/2)]
sub={}
for i in range(4):
    zi=phase[i]*W_[i]; sub[x[i]]=sp.re(sp.expand_complex(zi)); sub[w[i]]=sp.im(sp.expand_complex(zi))
G=sp.expand_complex(F0.subs(sub,simultaneous=True))
cc=sp.cos(th/2)
ab=lambda e: sp.expand_complex(e*sp.conjugate(e))
S=2*cc**2*(ab(W_[0]+W_[2])+ab(W_[1]+W_[3]))+(1+cc)**2*ab(W_[0]+W_[1]-W_[2]-W_[3])+(1-cc)**2*ab(W_[0]-W_[1]-W_[2]+W_[3])
print("pair SOS identity:",sp.simplify(sp.expand(sp.expand_trig(sp.expand_complex(G-S)).rewrite(sp.cos)))==0)
