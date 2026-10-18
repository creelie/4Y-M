# Derive the pair block M of Lemma 3.2 from the four position rules and compare.
import sympy as sp
chi,th=sp.symbols('chi theta',real=True)
wa=sp.exp(sp.I*(chi+th)); wb=sp.exp(sp.I*(chi-th))
x=sp.symbols('x1:5',real=True); y=sp.symbols('y1:5',real=True)
z=[x[i]+sp.I*y[i] for i in range(4)]
cj=sp.conjugate
def YYW(i,j): return -sp.re(z[i]*cj(z[j])*wa+cj(z[i])*z[j]*wb)
def YWY(i,j): return -sp.re(z[i]*cj(z[j])*wb+cj(z[i])*z[j]*wa)
def ip(i,j): return 2*sp.re(cj(z[i])*z[j])
f=sum(YYW(i,i) for i in range(4))+2*YYW(0,1)+2*YYW(2,3)-2*sum(YWY(i,j) for i in (0,1) for j in (2,3))
c=2*ip(0,3)-2*ip(0,1)+2*sp.re(2*cj(z[1])*z[2]*wa*cj(wb))-2*ip(2,3)
F=sp.expand(2*sum(2*(x[i]**2+y[i]**2) for i in range(4))-f-c)
m=sp.cos(chi); e=sp.exp
M=sp.zeros(4,4)
for i in range(4): M[i,i]=4+2*m*sp.cos(th)
M[0,1]=2+2*m*e(-sp.I*th); M[2,3]=M[0,1]
M[0,2]=-2*m*e(sp.I*th); M[1,3]=M[0,2]
M[0,3]=-2-2*m*e(sp.I*th); M[1,2]=-2*e(2*sp.I*th)-2*m*e(sp.I*th)
for i in range(4):
    for j in range(i):
        M[i,j]=cj(M[j,i])
zz=sp.Matrix(z)
G=sp.expand((zz.H*M*zz)[0])
print('F - z*Mz =', sp.simplify(sp.expand_complex(F-G)))
