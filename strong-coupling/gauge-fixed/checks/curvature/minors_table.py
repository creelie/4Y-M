import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from exact4 import *
u=sp.symbols('u')
for i0 in range(4):
    M=pair_matrix(i0,1,1)
    Ms=sp.simplify(M.subs(S,sp.sqrt(1-C**2)))
    lead=[M.extract(list(range(k)),list(range(k))).det() for k in (1,2,3)]
    out=[]
    for mi in lead:
        v=redS(mi); re_=redS(sp.re(v)); im=redS(sp.im(v))
        out.append(sp.factor(sp.expand(re_).subs(C**2,u)) if not re_.has(S) else ('S!',re_))
        assert im==0
    print('position',i0+1,'pair block (m=1), leading minors in u=cos^2(theta/2):',out)
    print('   matrix:',[[sp.simplify(sp.expand(Ms[r,c])) for c in range(3)] for r in range(3)])
    for cphi in (1,-1):
        D=diag_matrix(i0,1,cphi)
        print('   diag block cos phi=%d:'%cphi,D.tolist(),'leading minors',[D.extract(list(range(k)),list(range(k))).det() for k in (1,2,3)])
