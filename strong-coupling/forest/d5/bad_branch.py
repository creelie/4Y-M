"""Scalar diagnostic and rational sign bound for BAD_BRANCHES.md."""
from fractions import Fraction as F
from math import sin,cos,sqrt

left=F(1,2);right=F(49,100)
# Strict elementary bounds sin(1/4)<1/4, cos(1/4)>31/32,
# 1<sqrt(2)<3/2 imply this strict upper bound for f'(0).
upper=-F(3,4)+F(7,2)*F(3,2)*F(1,4)-F(61,100)*F(31,32)
assert upper == -F(91,3200) and upper<0
sine_lower=F(1,4)-F(1,4)**3/6
assert 2*F(7,5)*sine_lower>F(1,2)
assert 2*F(3,2)*F(1,4)<1

def f(q):
    differences=[-.5-q]+[.49-q]*6
    m=max(2*sqrt(2)*abs(sin(d/2)) for d in differences)
    chi=max(0.,min(1.,2-2*m))
    return q+chi*sum(differences)/8

gamma=f(0)
fp=-.75+sqrt(2)*(3.5*sin(.25)-.61*cos(.25))
print('Rational proof: cutoff is in its transition branch at q=0')
print('Rational strict upper bound for f\'(0): -91/3200 < 0')
print(f'Numerical f(0)={gamma:.12g}; f\'(0)={fp:.12g}')
# Locate the two additional scalar preimages; q=0 is the middle one.
def bisect(lo,hi):
    flo=f(lo)-gamma
    for _ in range(80):
        mid=(lo+hi)/2; fm=f(mid)-gamma
        if flo*fm<=0:hi=mid
        else:lo=mid;flo=fm
    return (lo+hi)/2
roots=[bisect(-1,-1e-3),0.,bisect(1e-3,1)]
for q in roots:
    assert abs(f(q)-gamma)<1e-12
    print(f'path angle q={q:.12g}, average angle f(q)={f(q):.12g}')
print('Three path-link solutions for the same seven non-path links and average.')
