"""Single SU(2)-variable gradient-flow calibration, not four-dimensional Wilson flow."""
from fractions import Fraction as F
import json
from pathlib import Path
import time
import numpy as np
from numpy.polynomial.legendre import leggauss

start=time.monotonic()
mom=[F(1)]
for k in range(1,5): mom.append(mom[-1]*(F(3,2)+k-1))
ew=mom[2]/24-mom[1]/6
exw=mom[3]/24-mom[2]/6
exxw=mom[4]/24-mom[3]/6
cumulant=exxw-mom[2]*ew-2*mom[1]*(exw-mom[1]*ew)
assert cumulant==F(7,16)
assert mom[3]-mom[1]*mom[2]==F(15,2)
assert F(7,16)+F(15,2)*F(1,6)==F(27,16)

def correlations(g,tau,n):
    r,w=leggauss(n);r=6*(r+1);w=6*w
    theta=g*r/np.sqrt(2)
    Q=8*np.sin(theta/2)**2
    weight=w*r*r*(np.sin(theta)/theta)**2*np.exp(-Q/g**2)
    weight/=weight.sum()
    c=np.exp(-2*tau)
    flowed=2*np.arctan(c*np.tan(theta/2))
    A=8*np.sin(flowed/2)**2/g**2
    B=Q/g**2
    cov=lambda X,Y:float(np.dot(weight,X*Y)-np.dot(weight,X)*np.dot(weight,Y))
    return cov(A,B)/cov(B,B)

rows=[]
for tau in [.1,.25,.5]:
    c2=np.exp(-4*tau);expected=5/8*c2*(1-c2)
    gs=np.array([.18,.15,.12,.09,.06])
    vals=np.array([correlations(g,tau,256) for g in gs])
    checks=np.array([correlations(g,tau,512) for g in gs])
    assert max(abs(vals-checks))<1e-11
    effective=(vals-c2)/gs**2
    extrap=float(np.polynomial.polynomial.polyfit(gs**2,effective,3)[0])
    assert abs(extrap-expected)<1e-6
    rows.append({'flow_time':tau,'R0':float(c2),'R2_full_exact_formula':float(expected),
                 'R2_linear_flow_only':float(expected/3),'R2_extrapolated':extrap,
                 'g':gs.tolist(),'effective_R2':effective.tolist(),
                 'quadrature_max_difference':float(max(abs(vals-checks)))})
report={'classification':'exact rational Wick factors plus floating compact-group quadrature calibration only',
        'measure_cumulant':str(cumulant),'observable_covariance':str(mom[3]-mom[1]*mom[2]),'rows':rows,'seconds':time.monotonic()-start}
Path('result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,sort_keys=True))
