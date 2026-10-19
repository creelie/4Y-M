import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from lip import *
U0=haar(100000)
def galerkin_gap(b,D=2):
    X,F,G=features(U0,D)
    w=np.exp(np.einsum('nij,ji->n',U0,b).real); w/=w.sum()
    mF=w@F; C=(F-mF).T@((F-mF)*w[:,None]); Gw=(G*np.sqrt(w)[None,:,None]).reshape(-1,G.shape[2]); A=Gw.T@Gw
    ev,V=np.linalg.eigh(C); keep=ev>1e-8*ev.max(); P=V[:,keep]/np.sqrt(ev[keep])
    lam=np.linalg.eigvalsh(P.T@A@P); return lam.min()
def mgrad(b):
    X,F,G=features(U0,1)
    w=np.exp(np.einsum('nij,ji->n',U0,b).real); w/=w.sum()
    Gw=(G*np.sqrt(w)[None,:,None]).reshape(-1,18); A=Gw.T@Gw
    # f_h = Re Tr(U h) is linear in X with |coefficient| = |h|_HS; Re Tr(Uh)=sum Re(U_ij h_ji)
    return np.linalg.eigvalsh(A).max()
r=0.732
for name,b in [('rI',r*np.eye(3)),('r(1,1,0)',np.diag([r,r,0])),('r(1,1,0.5i)',np.diag([r,r,0.5j*r])),('r(1,.5,.5)',np.diag([r,r/2,r/2]))]:
    print(name,'Galerkin gap (upper bound) %.4f   bound used 1.0750'%galerkin_gap(b),
          '  nu|grad f|^2 max %.4f  bound 0.82'%mgrad(b))
