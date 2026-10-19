"""Exact check of the gauge-fixed local lemma for each position of the frozen link.

S_p = sum_i alpha_i |Y_i|^2 -/+ f_p - t c_p, restricted to the three active links;
alpha = a on the two links adjacent to the frozen one, b on the opposite link.
Pair block (z_i = (Y_i)_{ab}, w_a = m q^2, w_b = m qbar^2, q = e^{i theta/2}, m = cos chi = +-1)
and diagonal block (phase cos phi = +-1); both are affine in m and cos phi.
Rules (Appendix A of the paper):
  ReTr(Y_iY_jW) -> -Re(z_i zb_j w_a + zb_i z_j w_b),  ReTr(Y_iWY_j) -> -Re(z_i zb_j w_b + zb_i z_j w_a),
  <Y_i,Y_j> -> 2Re(zb_i z_j),  <Y_2,WY_3W*> -> 2Re(zb_2 z_3 w_a wb_b),  |Y_i|^2 -> 2|z_i|^2.
"""
import sympy as sp, itertools, sys, math
a,b,t=[sp.Rational(x) for x in (sys.argv[1:4] if len(sys.argv)>3 else ('21/10','5/4','9/10'))]
C,S=sp.symbols('C S',real=True)
q=C+sp.I*S; qb=C-sp.I*S
z=sp.symbols('z1:5'); zb=sp.symbols('zb1:5')
y=sp.symbols('y1:5',real=True)
def Re2(e,ec): return (e+ec)   # 2 Re(e), ec = conjugate expression
def build_pair(m):
    wa=m*q**2; wb=m*qb**2; cwa=m*qb**2; cwb=m*q**2
    YYW=lambda i,j: -sp.Rational(1,2)*Re2(z[i]*zb[j]*wa+zb[i]*z[j]*wb, zb[i]*z[j]*cwa+z[i]*zb[j]*cwb)
    YWY=lambda i,j: -sp.Rational(1,2)*Re2(z[i]*zb[j]*wb+zb[i]*z[j]*wa, zb[i]*z[j]*cwb+z[i]*zb[j]*cwa)
    ip=lambda i,j: Re2(zb[i]*z[j], z[i]*zb[j])                      # <Y_i,Y_j>
    ip23=Re2(zb[1]*z[2]*wa*cwb, z[1]*zb[2]*cwa*wb)                   # <Y_2,WY_3W*>
    f=sum(YYW(i,i) for i in range(4))+2*YYW(0,1)+2*YYW(2,3)-2*sum(YWY(i,j) for i in (0,1) for j in (2,3))
    corner={(0,3):2*ip(0,3),(0,1):-2*ip(0,1),(1,2):2*ip23,(2,3):-2*ip(2,3)}
    return f,corner
opp={0:2,1:3,2:0,3:1}
def pair_matrix(i0,sgn,m):
    f,corner=build_pair(m)
    act=[i for i in range(4) if i!=i0]
    al={i:(b if i==opp[i0] else a) for i in act}
    e=sum(al[i]*2*z[i]*zb[i] for i in act)-sgn*f-t*sum(v for k,v in corner.items() if i0 not in k)
    e=sp.expand(e.subs({z[i0]:0,zb[i0]:0}))
    return sp.Matrix(3,3,lambda r,c: sp.expand(e.coeff(zb[act[r]]).coeff(z[act[c]])))
def diag_matrix(i0,sgn,cphi):
    act=[i for i in range(4) if i!=i0]
    al={i:(b if i==opp[i0] else a) for i in act}
    s=y[0]+y[1]-y[2]-y[3]
    cd={(0,3):2*y[0]*y[3],(0,1):-2*y[0]*y[1],(1,2):2*y[1]*y[2],(2,3):-2*y[2]*y[3]}
    e=sum(al[i]*y[i]**2 for i in act)+sgn*cphi*s**2-t*sum(v for k,v in cd.items() if i0 not in k)
    e=sp.expand(e.subs(y[i0],0))
    return sp.Matrix(3,3,lambda r,c: sp.Rational(1,2)*sp.diff(e,y[act[r]],y[act[c]]))
def minors(M):
    return [M.extract(list(Sx),list(Sx)).det() for k in (1,2,3) for Sx in itertools.combinations(range(3),k)]
def redS(v):
    v=sp.expand(v)
    for _ in range(12): v=sp.expand(v.subs(S**2,1-C**2))
    return v
if __name__=="__main__":
    ok=True; polys=set()
    for i0 in range(4):
        for sgn in (1,-1):
            for m in (1,-1):
                M=pair_matrix(i0,sgn,m)
                herm=sp.simplify(redS(M-M.H.subs({})))
                for mi in minors(M):
                    v=redS(mi)
                    im=redS(sp.im(v)); re_=redS(sp.re(v))
                    if im!=0: ok=False; print('nonreal minor',i0,sgn,m,im)
                    if re_.has(S):
                        vals=[float(re_.subs({C:math.cos(x),S:math.sin(x)})) for x in [math.pi*k/1000 for k in range(2000)]]
                        print('minor depends on S',i0,sgn,m,min(vals)); ok=ok and min(vals)>0
                        continue
                    P=sp.Poly(re_,C)
                    roots=[r for r in sp.real_roots(P) if -1<=r<=1]
                    if roots or P.eval(0)<=0: ok=False; print('FAIL',i0,sgn,m,sp.factor(re_))
                    polys.add(sp.factor(re_))
            for cphi in (1,-1):
                for mi in minors(diag_matrix(i0,sgn,cphi)):
                    if mi<=0: ok=False; print('FAIL diag',i0,sgn,cphi,mi)
        print('position',i0,'checked',flush=True)
    print('distinct pair minors:',len(polys))
    for p in sorted(polys,key=str): print('  ',p)
    print('ALL PRINCIPAL MINORS POSITIVE:',ok)
