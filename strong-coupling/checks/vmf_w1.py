# W1 (chordal) Lipschitz rate of the von Mises-Fisher law on S^3 under a change of
# its parameter M in the direction of M itself (zonal case), by an exact LP over
# zonal 1-Lipschitz test functions.  Also Var(<x,u>) for parallel/perpendicular u.
import numpy as np
from scipy.optimize import linprog
from scipy.integrate import quad
def moments(k):
    w=lambda t: np.exp(k*np.cos(t))*np.sin(t)**2
    Z=quad(w,0,np.pi)[0]
    E1=quad(lambda t:w(t)*np.cos(t),0,np.pi)[0]/Z
    E2=quad(lambda t:w(t)*np.cos(t)**2,0,np.pi)[0]/Z
    return E1,E2
def lp_rate(k,n=240):
    th=(np.arange(n)+0.5)*np.pi/n
    w=np.exp(k*np.cos(th))*np.sin(th)**2; w/=w.sum()
    g=np.cos(th); g=g-(w*g).sum()          # derivative of density along M-hat
    c=-(w*g)                                # maximise sum w g f
    A=[];b=[]
    for i in range(n-1):                   # adjacent constraints suffice? no: use all pairs
        pass
    I,J=np.triu_indices(n,1)
    d=2*np.sin((th[J]-th[I])/2)
    m=len(I)
    from scipy.sparse import coo_matrix
    rows=np.concatenate([np.arange(m),np.arange(m),m+np.arange(m),m+np.arange(m)])
    cols=np.concatenate([I,J,I,J]); vals=np.concatenate([np.ones(m),-np.ones(m),-np.ones(m),np.ones(m)])
    A=coo_matrix((vals,(rows,cols)),shape=(2*m,n))
    res=linprog(c,A_ub=A,b_ub=np.concatenate([d,d]),bounds=[(None,None)]*n,method='highs')
    return -res.fun
for k in [0,0.25,0.5,1,2,3,5,8]:
    E1,E2=moments(k)
    var_par=E2-E1**2; var_perp=(1-E2)/3
    print(f"kappa={k:5}  Var_par={var_par:.5f}  Var_perp={var_perp:.5f}  LP_zonal_rate={lp_rate(k,160):.5f}")
