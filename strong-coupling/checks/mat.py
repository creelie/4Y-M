import sympy as sp
exec(open('sos.py').read().split('# SOS at chi=0')[0].replace('print(','(lambda *a:None)('))
# Hermitian matrix: F = sum conj(z_i) M_ij z_j ; recover via polarization
M=sp.zeros(4,4)
zs=[x[i]+sp.I*w[i] for i in range(4)]
def Fz(vec):
    sub={}
    for i in range(4): sub[x[i]]=sp.re(vec[i]); sub[w[i]]=sp.im(vec[i])
    return sp.simplify(F.subs(sub,simultaneous=True))
for i in range(4):
    e=[0]*4; e[i]=1; M[i,i]=sp.simplify(Fz(e))
for i in range(4):
  for j in range(i+1,4):
    e=[0]*4; e[i]=1; e[j]=1; r=Fz(e)-M[i,i]-M[j,j]   # = 2 Re M_ij
    e=[0]*4; e[i]=1; e[j]=sp.I; s=Fz(e)-M[i,i]-M[j,j] # = 2 Re(i M_ij) = -2 Im M_ij
    Mij=sp.simplify((r - sp.I*s)/2)
    M[i,j]=sp.simplify(sp.expand_complex(Mij).rewrite(sp.exp)); M[j,i]=sp.conjugate(M[i,j])
for i in range(4):
  for j in range(4):
    print(i+1,j+1, sp.simplify(sp.expand_complex(M[i,j]).rewrite(sp.cos)))
