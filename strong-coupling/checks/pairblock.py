"""Unordered-pair block {a,b} (a!=b): form S(z)=q_{u,v}(z)+q_{v,u}(-conj z) on C^4, need S <= 4|z|^2.
Cartan block a=b: z=i*y, y real: need q_{u,u}(iy) <= 2|y|^2."""
import numpy as np
def Hq(u,v):
    H=np.zeros((4,4),complex)
    def off(i,j,c): H[j,i]+=c; H[i,j]+=np.conj(c)
    off(0,1,-(1+u)); off(2,3,-(1+u)); off(0,2,v); off(0,3,1+v); off(1,2,v*(1+np.conj(u))); off(1,3,v)
    return H-np.real(u)*np.eye(4)
def pair(u,v):
    # q_{v,u}(-conj z) = conj(z)^T Hq(v,u) conj(z) = z^T conj(Hq(v,u))... as Hermitian form in z: Hq(v,u)^T
    return Hq(u,v)+Hq(v,u).T
def cartan(u):
    return np.real(Hq(u,u))   # for z=iy, z^*Hz = y^T H y, real part -> Re H
if __name__=="__main__":
    rng=np.random.default_rng(0); m1=-1e9; m2=-1e9; arg=None
    for _ in range(300000):
        a,g=rng.uniform(0,2*np.pi,2); u,v=np.exp(1j*a),np.exp(1j*g)
        e=np.linalg.eigvalsh(pair(u,v))[-1]
        if e>m1: m1,arg=e,(a,g)
        m2=max(m2,np.linalg.eigvalsh(cartan(u))[-1])
    print("max eig pair block (need <=4):",m1,"at",arg)
    print("max eig cartan block (need <=2):",m2)
