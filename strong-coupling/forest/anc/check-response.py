"""Exact finite checks for the response appendix; standard library only.
These do not evaluate continuous Wilson integrals or prove the continuum mass gap.
"""
from fractions import Fraction as F
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parent
checks=[]
def check(name,ok):
    assert ok,name
    checks.append(name)
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def mv(A,v):return [dot(row,v) for row in A]
def mul(A,B):return [[dot(row,col) for col in zip(*B)] for row in A]
def add(A,B):return [[x+y for x,y in zip(a,b)] for a,b in zip(A,B)]
def scale(c,A):return [[c*x for x in row] for row in A]
def ident(n):return [[F(int(i==j)) for j in range(n)] for i in range(n)]
def trace(A):return sum(A[i][i] for i in range(len(A)))
def transpose(A):return [list(c) for c in zip(*A)]

# Gaussian Sylvester and noncommuting affine-transport fixture.
C=[[F(1),F(0)],[F(0),F(2)]]; Cp=[[F(0),F(1)],[F(1),F(0)]]
Sigma=[[F(1),F(0)],[F(0),F(1,2)]]
K=[[F(0),F(-1,3)],[F(-1,3),F(0)]]
H=[[F(0),F(-1,4)],[F(-1,4),F(0)]]
check('noncommuting_precision',mul(C,Cp)!=mul(Cp,C))
check('symmetric_Sylvester',add(mul(C,K),mul(K,C))==scale(-1,Cp))
check('optimal_cost',trace(mul(mul(K,Sigma),K))==F(1,6))
A0=[[F(0),F(-1,2)],[F(-1,4),F(0)]]
check('alternative_cost',trace(mul(mul(A0,Sigma),transpose(A0)))==F(3,16))
check('cost_difference',F(3,16)-F(1,6)==F(1,48))
R=add(Cp,add(mul(C,H),mul(H,C)))
check('Sylvester_residual',R==[[F(0),F(1,4)],[F(1,4),F(0)]])
Delta=add(H,scale(-1,K))
check('residual_response',trace(mul(mul(Delta,Delta),Sigma))==F(1,96))
check('residual_upper',F(1,96)<=F(1,32))

# Exact mixed Poisson/flux certificate.
L=[[F(2),F(-1)],[F(-1),F(2)]]
Li=[[F(2,3),F(1,3)],[F(1,3),F(2,3)]]
k=[F(1),F(0)];h=[F(1,2),F(0)]
r=[x-y for x,y in zip(k,mv(L,h))]
V=2*dot(k,h)-dot(h,mv(L,h));X=dot(k,mv(Li,k))
j=[F(-1,6),F(1,6),F(1,3)]
check('mixed_residual',r==[F(0),F(1,2)])
check('mixed_response',X==F(2,3))
check('mixed_identity',X==V+dot(r,mv(Li,r)))
check('exact_flux',[j[0]+j[1],-j[0]+j[2]]==r)
check('exact_flux_bound',V+dot(j,j)==X)
check('nonoptimal_flux_bound',V+F(1,4)==F(3,4)>=X)

# Genuine conditional expectations for the four-state sign fixture.
w=[F(3,8),F(1,8),F(1,8),F(3,8)]
P1=[[F(3,4),0,F(1,4),0],[0,F(1,4),0,F(3,4)],
    [F(3,4),0,F(1,4),0],[0,F(1,4),0,F(3,4)]]
P2=[[F(3,4),F(1,4),0,0],[F(3,4),F(1,4),0,0],
    [0,0,F(1,4),F(3,4)],[0,0,F(1,4),F(3,4)]]
I=ident(4); W=[[w[i] if i==j else F(0) for j in range(4)] for i in range(4)]
Q1=add(I,scale(-1,P1));Q2=add(I,scale(-1,P2))
L=add(Q1,Q2);T=scale(F(1,2),add(P1,P2))
ip=lambda a,b:sum((p*x*y for p,x,y in zip(w,a,b)),F(0))
s=[F(1),F(1),F(-1),F(-1)];t=[F(1),F(-1),F(1),F(-1)]
for index,P in enumerate([P1,P2],1):
    check('projection_'+str(index),mul(P,P)==P)
    check('weighted_symmetry_'+str(index),mul(W,P)==mul(transpose(P),W))
    check('constant_preservation_'+str(index),mv(P,[F(1)]*4)==[F(1)]*4)
vectors=[([s[i]+t[i] for i in range(4)],F(1,2)),
         ([s[i]-t[i] for i in range(4)],F(3,2)),
         ([s[i]*t[i]-F(1,2) for i in range(4)],F(2))]
for z,(v,lam) in enumerate(vectors):check('eigenpair_'+str(z),mv(L,v)==[lam*x for x in v])
check('cross_covariance',ip(s,t)==F(1,2))
h=[F(0)]*4;power=t[:];rho=F(5,6)
for m in range(1,9):
    h=[a+b/2 for a,b in zip(h,power)]
    power=mv(T,power)
    residual=[a-b for a,b in zip(t,mv(L,h))]
    check('Neumann_residual_'+str(m),residual==power)
    check('residual_rate_'+str(m),ip(residual,residual)<=rho**(2*m)*ip(t,t))
    resolved=ip(mv(Q1,s),mv(Q1,h))+ip(mv(Q2,s),mv(Q2,h))
    check('static_covariance_identity_'+str(m),ip(s,t)==resolved+ip(s,residual))
    if m==1:
        check('disjoint_support_resolved',resolved==0)
        check('first_residual_covariance',ip(s,residual)==F(1,2))

# Entropy, convexification polynomial identities and two-level arithmetic.
c0=F(487,5000)
check('entropy_curvature',F(5,16)-F(3,8)==F(-1,16))
check('indefinite_joint_Hessian',c0/F(8)-2*c0*c0==F(-169963,25000000))
check('two_level_trace',1/c0+240*(1+F(4,27)/c0)==F(2696920,4383))
check('exponential_comparison',F(65,24)+F(1,100)==F(1631,600)<F(25,9))
def pmul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out
p=[F(1),F(57),F(-15),F(1)]
a=pmul([F(49),F(-14),F(1)],[F(-1),F(1)]);a[0]+=50;a[1]-=6
b=pmul([F(64),F(-16),F(1)],[F(1),F(1)]);b[0]-=63;b[1]+=9
check('convexification_middle_interval',a==p)
check('convexification_last_interval',b==p)
for n in [2,3,4,8]:
    v=[F(1)]+[F(1,n)]*(n*n)
    check('projector_norm_'+str(n),dot(v,v)==2)
    check('projector_row_'+str(n),sum((abs(v[0]*x/2) for x in v),F(0))==F(n+1,2))
result={'status':'passed','exact_checks':len(checks),'arithmetic':'Python fractions.Fraction; no floating-point certification','checks':checks,'scope':'Finite fixtures for the written response arguments. No continuous Wilson integrals, new Lean build, uniform D estimate, continuum E construction or physical mass gap is certified by this script.'}
(ROOT/'response-checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'status':'passed','exact_checks':len(checks)}))
