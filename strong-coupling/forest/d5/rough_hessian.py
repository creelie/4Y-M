"""Search a centre-valued canonical fibre for an exact negative action direction.

Floating eigensolver proposes an integer vector; all final constraints and
signed-curvature sums are checked with integer arithmetic. Q''=2 sum sign*(d xi)^2
for a unit Hilbert--Schmidt generator. No sampling or continuum claim.
"""
from itertools import product, combinations
import json
import numpy as np
from scipy.linalg import eigh

n = 2
fine = 2*n
sites = list(product(range(fine), repeat=4))
coarse = list(product(range(n), repeat=4))
def shift(x, mu, step=1):
    out = list(x); out[mu] = (out[mu]+step)%fine; return tuple(out)
def tree(x, mu):
    return x[mu]%2 == 0 and all(x[nu]%2 == 0 for nu in range(mu+1,4))
def path(x, mu):
    return x[mu]%2 == 1 and all(x[nu]%2 == 0 for nu in range(4) if nu != mu)
edges = [(x,mu) for x in sites for mu in range(4) if not tree(x,mu)]
ix = {e:i for i,e in enumerate(edges)}
free = [e for e in edges if not path(*e)]
fidx = {e:i for i,e in enumerate(free)}
# Integer matrix S2=2S avoids halves in the constraint section.
S2 = np.zeros((len(edges),len(free)), dtype=np.int64)
for e,j in fidx.items(): S2[ix[e],j]=2
for y in coarse:
    for mu in range(4):
        origin=tuple(2*v for v in y)
        dp=shift(origin,mu)
        row=ix[(dp,mu)]
        for o in product(range(2),repeat=4):
            x=tuple(2*y[k]+o[k] for k in range(4))
            for e in ((x,mu),(shift(x,mu),mu)):
                if e in fidx: S2[row,fidx[e]] -= 1
        # The path link occurs twice in the average, hence these coefficients.
# canonical c_y,mu=(-1)^(sum_{nu<mu} y_nu) has every coarse plaquette -I.
def sign_link(x,mu):
    if x[mu]%2 == 0: return 1
    return (-1)**sum(x[nu]//2 for nu in range(mu))
rows=[]; signs=[]; plaqs=[]
for x in sites:
    for mu,nu in combinations(range(4),2):
        es=[((x,mu),1),((shift(x,mu),nu),1),((shift(x,nu),mu),-1),((x,nu),-1)]
        row=np.zeros(len(edges),dtype=np.int64)
        sg=1
        for e,s in es:
            if e in ix: row[ix[e]]+=s
            sg*=sign_link(*e)
        rows.append(row); signs.append(sg); plaqs.append((x,mu,nu))
D=np.array(rows,dtype=np.int64); signs=np.array(signs,dtype=np.int64)
DS=D@S2
H=DS.T@(signs[:,None]*DS)
val,vec=eigh(H.astype(float),subset_by_index=[0,0])
print(f'coarse side {n}; active {len(edges)}; independent {len(free)}')
print(f'min eigenvalue of S2^T d^T sign d S2 (computation): {val[0]:.12g}')
if val[0] >= -1e-8: raise RuntimeError('No negative direction found')
z=np.rint(100*vec[:,0]/max(abs(vec[:,0]))).astype(np.int64)
xi=S2@z
# All 64 exact linear average constraints, integer numerator of 16V.
constraints=[]
for y in coarse:
    for mu in range(4):
        c=0
        for o in product(range(2),repeat=4):
            x=tuple(2*y[k]+o[k] for k in range(4))
            for e in ((x,mu),(shift(x,mu),mu)):
                if e in ix: c+=int(xi[ix[e]])
        constraints.append(c)
assert all(c==0 for c in constraints)
dxi=D@xi
positive=sum(int(t)**2 for t,s in zip(dxi,signs) if s==1)
negative=sum(int(t)**2 for t,s in zip(dxi,signs) if s==-1)
assert positive-negative<0
record={'coarse_side':n,'generator_HS_norm':1,'z':z.tolist(),
        'positive_sum':positive,'negative_sum':negative,'signed_sum':positive-negative,
        'norm_squared':sum(int(t)**2 for t in xi),
        'all_average_numerators_zero':True,
        'meaning':'xi=S2*z; Q second derivative = 2*signed_sum'}
print(json.dumps({k:v for k,v in record.items() if k!='z'},sort_keys=True))
print('Discovery only; the sparse exact witness is in rough_certificate.py')
