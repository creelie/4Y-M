# Chordal W1-Lipschitz rate of the vMF family on S^3 in direction u=cos(phi)e1+sin(phi)e2,
# mean direction e1, concentration k.  Reduced to the half 2-sphere (x1,x2,r), weight exp(k x1) r.
import numpy as np, ot, sys
def grid(n):
    # (theta,psi) grid on half-sphere: x1=cos t, x2=sin t cos p, r=sin t sin p, p in (0,pi)
    t=(np.arange(n)+0.5)*np.pi/n; p=(np.arange(n)+0.5)*np.pi/n
    T,P=np.meshgrid(t,p,indexing='ij')
    X=np.column_stack([np.cos(T).ravel(),(np.sin(T)*np.cos(P)).ravel(),(np.sin(T)*np.sin(P)).ravel()])
    dA=(np.sin(T)).ravel()
    return X,dA
def rate(k,phi,n=40):
    X,dA=grid(n)
    w=np.exp(k*(X[:,0]-1))*X[:,2]*dA; w/=w.sum()
    g=np.cos(phi)*X[:,0]+np.sin(phi)*X[:,1]; g=g-(w*g).sum()
    s=w*g; a=np.clip(s,0,None); b=np.clip(-s,0,None)
    M=a.sum(); C=ot.dist(X,X,metric='euclidean')
    return ot.emd2(a/M,b/M,C,numItermax=10**7)*M
n=int(sys.argv[1]) if len(sys.argv)>1 else 40
for k in [0,0.5,1,2,3,5]:
    print(k,[round(rate(k,phi,n),5) for phi in np.linspace(0,np.pi/2,5)],flush=True)
