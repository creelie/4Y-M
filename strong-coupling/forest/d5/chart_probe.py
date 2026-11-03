"""Nonabelian local-pair diagnostic for CANONICAL_CHART.md (not a certificate).

Finite differences independently compare the physical path-link Jacobian
with H=I+N_X. Exponential coordinates use an HS-orthonormal SU(2) basis.
"""
import numpy as np
from math import sin, cos, atan2, sqrt

def mul(a,b):
    return np.r_[a[0]*b[0]-a[1:]@b[1:],a[0]*b[1:]+b[0]*a[1:]+np.cross(a[1:],b[1:])]
def inv(a): return np.r_[a[0],-a[1:]]
def exp(x):
    r=np.linalg.norm(x)/sqrt(2)
    return np.r_[cos(r), x*(sin(r)/(sqrt(2)*r) if r else 1/sqrt(2))]
def log(q):
    s=np.linalg.norm(q[1:])
    return q[1:]*(sqrt(2)*atan2(s,q[0])/s if s else sqrt(2))

def run(scale):
    rng=np.random.default_rng(20261103)
    c=exp(np.array([.3,-.2,.1]))
    # Pair mu=0: offset zero has tree u0=w0=I, and all other u,w are active.
    u=[exp(scale*rng.normal(size=3)) for _ in range(7)]
    w=[exp(scale*rng.normal(size=3)) for _ in range(7)]
    ell=[mul(c,exp(scale*rng.normal(size=3))) for _ in range(7)]
    A=[mul(inv(c),mul(ua,la)) for ua,la in zip(u,ell)]
    A += [mul(inv(c),mul(la,wa)) for la,wa in zip(ell,w)]
    def N(X): return sum((log(mul(exp(-X),a)) for a in A),np.zeros(3))/16
    X=np.zeros(3)
    for iteration in range(1000):
        new=-N(X)
        if np.linalg.norm(new-X)<1e-16: X=new;break
        X=new
    residual=np.linalg.norm(X+N(X))
    base=mul(c,exp(X))
    def coarse_from_path(path):
        Nx=sum((log(mul(inv(path),mul(c,a))) for a in A),np.zeros(3))/16
        return mul(path,exp(Nx))
    errors=[]
    for h in [1e-4,1e-5,1e-6]:
        H=np.eye(3);J=np.zeros((3,3))
        for j in range(3):
            e=np.eye(3)[j]*h
            H[:,j]+=(N(X+e)-N(X-e))/(2*h)
            plus=coarse_from_path(mul(base,exp(e)))
            minus=coarse_from_path(mul(base,exp(-e)))
            J[:,j]=(log(mul(inv(c),plus))-log(mul(inv(c),minus)))/(2*h)
        err=abs(np.linalg.det(J)/np.linalg.det(H)-1)
        errors.append(err)
        print(f'scale={scale:g} step={h:g} detJ={np.linalg.det(J):.12g} detH={np.linalg.det(H):.12g} relerr={err:.3g}')
    assert residual<1e-13 and min(errors)<1e-7
    print(f'  fixed-point residual={residual:.3g}, iterations={iteration+1}')

for scale in [1e-3,1e-4,1e-5]: run(scale)
print('Nonabelian determinant diagnostic passed (floating computation only).')
