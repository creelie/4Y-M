import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from cert import forms
import numpy as np, cvxpy as cp
# F at position 4 (index 3): adjacent e1,e3 (idx 0,2), opposite e2 (idx 1); corners (0,1) sign -1, (1,2) sign +1
idx=np.array([0,1,2,3,4,5,6,7,8])
def S_of(a,b,t,th,sgn=1):
    f,C=forms(th)
    D=lambda i: np.diag([1.0 if 3*i<=k<3*i+3 else 0 for k in range(12)])
    S=a*(D(0)+D(2))+b*D(1)-sgn*f-t*(-C[(0,1)]+C[(1,2)])
    return S[:9,:9]
def sdp(nth):
    a,b,t=cp.Variable(),cp.Variable(),cp.Variable()
    cons=[]
    for th in np.linspace(0,2*np.pi,nth,endpoint=False):
        f,C=forms(th)
        D=lambda i: np.diag([1.0 if 3*i<=k<3*i+3 else 0 for k in range(12)])[:9,:9]
        S=a*(D(0)+D(2))+b*D(1)-f[:9,:9]-t*(-C[(0,1)]+C[(1,2)])[:9,:9]
        cons.append(0.5*(S+S.T)>>0)
    pr=cp.Problem(cp.Minimize(4*a+2*b+4*t),cons); pr.solve(solver='CLARABEL')
    return pr.status,pr.value,a.value,b.value,t.value
for n in [90,360,1440]:
    print(n,sdp(n),flush=True)
