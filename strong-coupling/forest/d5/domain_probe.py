"""Floating diagnostic of the proved asymptotic domain obstruction.

The proof is in LABEL_DOMAIN.md. These finite-scale values are not
interval certificates, and do not establish an explicit validity radius.
"""
import numpy as np
from scipy.optimize import root
from chart_probe import mul, exp, log

def solve(z, initial):
    A=[exp(za) for za in z]
    def H(X):
        return X+sum((log(mul(exp(-X),a)) for a in A),np.zeros(3))/8
    out=root(H,initial,tol=1e-11)
    residual=np.linalg.norm(H(out.x))
    assert residual<1e-12
    return out.x,residual

print('Nonconvex physical-label pullback: floating diagnostic only')
for radius in [.1,.03,.01,.003]:
    h=.02*radius
    a=radius+radius*h*h/6
    centre=np.tile(np.array([-a/7,0.,0.]),(7,1))
    v=np.zeros_like(centre);v[0,1]=1;v[1,1]=-1
    plus,rp=solve(centre+h*v,[a,0,0])
    minus,rm=solve(centre-h*v,[a,0,0])
    midpoint,r0=solve(centre,[a,0,0])
    free_max=max(np.linalg.norm(za) for za in centre+h*v)
    assert max(np.linalg.norm(plus),np.linalg.norm(minus))<radius
    assert np.linalg.norm(midpoint)>radius
    assert free_max<radius
    print(f'r={radius:g}: endpoint path radii {np.linalg.norm(plus):.15g}, {np.linalg.norm(minus):.15g}; midpoint {np.linalg.norm(midpoint):.15g}')
    print(f'  max free radius={free_max:.9g}; max equation residual={max(rp,rm,r0):.3g}')
print('At tested radii, endpoints are good and their coordinate midpoint is bad.')
