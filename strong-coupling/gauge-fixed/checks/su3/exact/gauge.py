import numpy as np
from scipy.optimize import linprog
def torus_grid(n=72):
    t=np.linspace(0,2*np.pi,n,endpoint=False); T1,T2=np.meshgrid(t,t); T1=T1.ravel(); T2=T2.ravel()
    return np.stack([np.exp(1j*T1),np.exp(1j*T2),np.exp(-1j*(T1+T2))],1)
TG=torus_grid()
AEQ=np.vstack([TG.real.T,TG.imag.T])
def gauge(z):
    z=np.asarray(z,complex); beq=np.concatenate([z.real,z.imag])
    res=linprog(np.ones(len(TG)),A_eq=AEQ,b_eq=beq,bounds=(0,None),method='highs')
    return res.fun if res.status==0 else np.inf
if __name__=='__main__':
    w=np.exp(1j*np.pi/3)
    for z in [[1,1,1],[0,0,1],[w/2]*3,[0.6*w]*3,[(1+w)/2,(1+w)/2,(1+w**-2)/2],[1,0,0],[0.5,0.5,0]]:
        print(np.round(z,3),gauge(z))
