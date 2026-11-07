from fractions import Fraction as F
from itertools import product

for k in range(12):
    a=F(1,2**k)
    a_next=a/2
    for j in range(1,101):
        assert a_next*max(a_next,F(1,j)) <= a*max(a,F(1,j))
r=F(1,2)
v=(1,-1)
masked_q=2*r*v[0]*v[1]
assert masked_q<0
print('PASS transfer comparison monotonicity at 1200 rational entries')
print('PASS cross-cell mask has a negative rational quadratic form')
