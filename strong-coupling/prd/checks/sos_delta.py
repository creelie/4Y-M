# Check M(-1,th)=M(1,th+pi) and the sum-of-squares form of z*M(1,th)z after z = Delta u.
import sympy as sp
th,m=sp.symbols('theta m',real=True)
e=sp.exp; I=sp.I
def Mat(m,th):
    M=sp.zeros(4,4)
    for i in range(4): M[i,i]=4+2*m*sp.cos(th)
    M[0,1]=2+2*m*e(-I*th); M[2,3]=M[0,1]
    M[0,2]=-2*m*e(I*th); M[1,3]=M[0,2]
    M[0,3]=-2-2*m*e(I*th); M[1,2]=-2*e(2*I*th)-2*m*e(I*th)
    for i in range(4):
        for j in range(i): M[i,j]=sp.conjugate(M[j,i])
    return M
print('reflection:',sp.simplify((Mat(-1,th)-Mat(1,th+sp.pi)).applyfunc(lambda t: sp.simplify(sp.expand_complex(t)))))
D=sp.diag(1,e(I*th/2),e(-I*th),e(-I*th/2))
N=(D.H*Mat(1,th)*D).applyfunc(lambda t: sp.simplify(sp.expand_complex(t)))
c=sp.cos(th/2)
v=lambda a: sp.Matrix(a)
S=2*c**2*(v([1,0,1,0])*v([1,0,1,0]).T+v([0,1,0,1])*v([0,1,0,1]).T)+(1+c)**2*v([1,1,-1,-1])*v([1,1,-1,-1]).T+(1-c)**2*v([1,-1,-1,1])*v([1,-1,-1,1]).T
print('SOS:',(N-S).applyfunc(lambda t: sp.simplify(sp.expand_complex(sp.expand_trig(t)))))
print(N)
