"""Actual finite spatial plaquette flow vector field, not the interacting Wick integral."""
import os
os.environ['JAX_PLATFORMS']='cpu'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE']='false'
import itertools as it
import math
import json
from pathlib import Path
import time
import numpy as np
import jax
jax.config.update('jax_enable_x64',True)
import jax.numpy as j

start=time.monotonic()
shape=(3,2,2,2)
vertices=list(it.product(*(range(n) for n in shape)))
edges=[(v,mu) for v in vertices for mu in range(4) if v[mu]+1<shape[mu]]
lookup={e:i for i,e in enumerate(edges)}
tree=set()
for v in vertices[1:]:
    mu=next(k for k in range(4) if v[k]);parent=list(v);parent[mu]-=1
    tree.add(lookup[(tuple(parent),mu)])
ids=[]
for v in vertices:
    for mu,nu in it.combinations(range(1,4),2):
        if v[mu]+1>=shape[mu] or v[nu]+1>=shape[nu]:continue
        vm=list(v);vm[mu]+=1;vn=list(v);vn[nu]+=1
        ids.append([lookup[(v,mu)],lookup[(tuple(vm),nu)],lookup[(tuple(vn),mu)],lookup[(v,nu)]])
ids=j.array(ids);signs=j.array([1.,1.,-1.,-1.])

def mul(a,b):
    return j.concatenate([(a[...,0]*b[...,0]-j.sum(a[...,1:]*b[...,1:],axis=-1))[...,None],
        a[...,:1]*b[...,1:]+b[...,:1]*a[...,1:]+j.cross(a[...,1:],b[...,1:])],axis=-1)

def exact_action(A,g):
    # Smooth scalar functions of r^2 at zero, to make derivatives on tree edges valid.
    X=A[ids]*signs[None,:,None];r2=j.sum(X*X,axis=-1)
    # Initial fields are bounded; Taylor-to-degree-16 is independently below roundoff here.
    theta2=g*g*r2/2
    scalar=sum((-1.)**k*theta2**k/float(math.factorial(2*k)) for k in range(9))
    sinc=sum((-1.)**k*theta2**k/float(math.factorial(2*k+1)) for k in range(9))
    qs=j.concatenate([scalar[...,None],(g/j.sqrt(2))*sinc[...,None]*X],axis=-1)
    q=qs[:,0]
    for slot in range(1,4):q=mul(q,qs[:,slot])
    return 4*j.sum(1-q[:,0])

def coefficients(A):
    X=A[ids]*signs[None,:,None];r2=j.sum(X*X,axis=-1)
    z=j.zeros_like(r2);v0=j.zeros_like(X)
    qs=j.stack([j.concatenate([j.ones_like(r2)[...,None],v0],axis=-1),
                j.concatenate([z[...,None],X/j.sqrt(2)],axis=-1),
                j.concatenate([-r2[...,None]/4,v0],axis=-1),
                j.concatenate([z[...,None],-r2[...,None]*X/(12*j.sqrt(2))],axis=-1),
                j.concatenate([r2[...,None]**2/96,v0],axis=-1)],axis=0)
    q=qs[:,:,0,:]
    for slot in range(1,4):
        q=j.stack([sum((mul(q[a],qs[k-a,:,slot,:]) for a in range(k+1)),j.zeros_like(q[0])) for k in range(5)])
    return -4*j.sum(q[2:5,:,0],axis=1)

grad_exact=jax.jit(jax.grad(exact_action,argnums=0))
grad_coeff=jax.jit(jax.jacfwd(coefficients))
rng=np.random.default_rng(20261103)
A=rng.normal(size=(len(edges),3))*.2
A[list(tree)]=0
A=j.asarray(A)
G2,G3,G4=np.asarray(grad_coeff(A))
AA=np.asarray(A);r2=np.sum(AA*AA,axis=-1)
ad2=2*(AA*np.sum(AA*G2,axis=-1)[:,None]-r2[:,None]*G2)
rows=[]
for g in [.2,.1,.05,.025]:
    grad=np.asarray(grad_exact(A,g))
    theta=g*np.sqrt(r2/2)
    factor=1/np.sinc(theta/np.pi)**2
    radial=AA*np.sum(AA*grad,axis=-1)[:,None]/np.where(r2>0,r2,1)[:,None]
    exact=-(radial+factor[:,None]*(grad-radial))/g**2
    pred=-G2-g*G3+g*g*(-G4+ad2/12)
    omitted=-G2-g*G3-g*g*G4
    rows.append({'g':g,'complete_error':float(np.linalg.norm(exact-pred)),
                 'omitted_metric_error':float(np.linalg.norm(exact-omitted)),
                 'tree_velocity_norm':float(np.linalg.norm(exact[list(tree)]))})
assert all(a['complete_error']/b['complete_error']>6 for a,b in zip(rows,rows[1:]))
assert rows[-1]['omitted_metric_error']>2*rows[-1]['complete_error']
assert rows[-1]['tree_velocity_norm']>.01
report={'classification':'finite spatial-action Taylor vector-field check, not flowed Wick covariance or B2',
        'shape':shape,'edges':len(edges),'spatial_plaquettes':len(ids),'tree_edges':len(tree),'rows':rows,
        'metric_correction_norm':float(np.linalg.norm(ad2/12)),'seconds':time.monotonic()-start}
Path('result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,sort_keys=True))
