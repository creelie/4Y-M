"""Batched SU(2) constraint Jacobians near path folds and hard-label boundaries.

Computation only. Small controls selected away from label boundaries need
not span even though all eight private links give a regular submersion.
No Wasserstein bound or divergent integrated response is inferred.
"""
import argparse
import os
p=argparse.ArgumentParser();p.add_argument('--device',choices=['cpu','gpu'],required=True)
p.add_argument('--batch',type=int,default=64);p.add_argument('--batches',type=int,default=1)
p.add_argument('--seed',type=int,default=20261103);p.add_argument('--reference',action='store_true')
args=p.parse_args()
os.environ['JAX_PLATFORMS']='cpu' if args.device=='cpu' else 'cuda'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE']='false'
import json
import math
import time
from pathlib import Path
import numpy as np
import jax
jax.config.update('jax_enable_x64',True)
import jax.numpy as j
SQ2=math.sqrt(2)

def mul(a,b):
    return j.concatenate(((a[...,0]*b[...,0]-j.sum(a[...,1:]*b[...,1:],axis=-1))[...,None],
         a[...,0,None]*b[...,1:]+b[...,0,None]*a[...,1:]+j.cross(a[...,1:],b[...,1:])),axis=-1)
def inv(q):return q*j.array([1.,-1.,-1.,-1.])
def qexp(z):
    t=j.sum(z*z,axis=-1)/2;r=j.sqrt(j.maximum(t,1e-30))
    c=j.where(t<1e-12,1-t/2+t*t/24,j.cos(r))
    s=j.where(t<1e-12,1-t/6+t*t/120,j.sin(r)/r)/SQ2
    return j.concatenate((c[...,None],s[...,None]*z),axis=-1)
def qlog(q):
    v=q[...,1:];s2=j.sum(v*v,axis=-1);s=j.sqrt(j.maximum(s2,1e-30))
    ratio=j.where(s2<1e-16,SQ2/(j.where(j.abs(q[...,0])>1e-12,q[...,0],1.)),SQ2*j.arctan2(s,q[...,0])/s)
    return ratio[...,None]*v

def value(ell,u,w):
    paths=j.concatenate((mul(u,ell[1:]),mul(ell[1:],w)),axis=0)
    relative=mul(inv(ell[0]),paths)
    distances=j.sqrt(j.maximum(2*j.sum((relative-j.array([1.,0.,0.,0.]))**2,axis=-1),1e-30))
    ordered=j.sort(distances)
    m=ordered[-1];chi=j.clip(2-2*m,0.,1.)
    N=j.sum(qlog(relative),axis=0)/16
    return mul(ell[0],qexp(chi*N)),m,ordered[-1]-ordered[-2]

def jacobian(ell,u,w):
    def shifted(delta):return value(mul(ell,qexp(delta)),u,w)[0]
    q,m,gap=value(ell,u,w)
    raw=jax.jacfwd(shifted)(j.zeros((8,3))).reshape(4,24)
    # Right-trivialise dq: sqrt(2)*vector(q^-1 dq).
    basis=j.eye(4)
    tangent=SQ2*jax.vmap(lambda v:mul(inv(q),v)[1:])(basis).T
    return tangent@raw,q,m,gap

@jax.jit
def batch_diagnostic(ells,us,ws,radii):
    J,c,m,gap=jax.vmap(jacobian)(ells,us,ws)
    relative=mul(inv(c)[:,None,:],ells[:,1:,:])
    distance=j.linalg.norm(qlog(relative),axis=-1)
    clearance=j.abs(distance-radii[:,None])
    eligible=clearance>radii[:,None]*.05
    mask=j.concatenate((j.ones((ells.shape[0],1),dtype=bool),eligible),axis=1)
    allowed=J*j.repeat(mask,3,axis=1)[:,None,:]
    full=j.linalg.svd(J,compute_uv=False)[:,-1]
    path=j.linalg.svd(J[:,:,:3],compute_uv=False)[:,-1]
    interior=j.linalg.svd(allowed,compute_uv=False)[:,-1]
    return j.stack((full,path,interior,m,gap,j.sum(eligible,axis=1)),axis=1),J,c

def nmul(a,b):
    return np.r_[a[0]*b[0]-a[1:]@b[1:],a[0]*b[1:]+b[0]*a[1:]+np.cross(a[1:],b[1:])]
def ninv(q):return q*np.array([1.,-1.,-1.,-1.])
def nexp(z):
    s=np.linalg.norm(z)/SQ2
    return np.r_[math.cos(s),z*(math.sin(s)/(SQ2*s) if s else 1/SQ2)]
def nlog(q):
    s=np.linalg.norm(q[1:]);return q[1:]*(SQ2*math.atan2(s,q[0])/s if s else SQ2)
def nvalue(ell,u,w):
    paths=[nmul(ua,la) for ua,la in zip(u,ell[1:])]+[nmul(la,wa) for la,wa in zip(ell[1:],w)]
    rel=[nmul(ninv(ell[0]),pa) for pa in paths]
    m=max(np.linalg.norm(q-np.array([1.,0,0,0]))*SQ2 for q in rel)
    chi=max(0.,min(1.,2-2*m));N=sum((nlog(q) for q in rel),np.zeros(3))/16
    return nmul(ell[0],nexp(chi*N))

def generate(batch,seed,reference=False):
    rng=np.random.default_rng(seed)
    eps=np.exp(rng.uniform(np.log(.02 if reference else 1e-7),np.log(.2),batch))
    radius=np.exp(rng.uniform(np.log(1e-5),np.log(.02),batch))
    perturb=eps*np.exp(rng.uniform(np.log(.0001),np.log(.3),batch))
    ells=[];us=[];ws=[]
    for e,r,sigma in zip(eps,radius,perturb):
        a=.5;split=.001
        q=2*math.asin((1-e)/(2*SQ2))-a-split
        chi=2*e
        N=(1-7*chi/8)/(2*SQ2*math.cos((a+split+q)/2))
        b=(8*N+a+7*q)/6
        left=[-a-split]+[b]*6;right=[-a+split]+[b]*6
        P=[nexp(np.array([SQ2*t,0,0])+sigma*rng.normal(size=3)) for t in left]
        Pp=[nexp(np.array([SQ2*t,0,0])+sigma*rng.normal(size=3)) for t in right]
        path=nexp(np.array([SQ2*q,0,0]))
        # Compute the target average before reparameterising the other private links.
        temp=np.array([path]+[np.array([1.,0,0,0])]*7)
        c=nvalue(temp,np.array(P),np.array(Pp))
        priv=[]
        for k in range(7):
            direction=np.array([1.,0,0])+(.3 if k else .001)*rng.normal(size=3)
            direction/=np.linalg.norm(direction)
            distance=r if k==0 else max(20*r,.1)
            priv.append(nmul(c,nexp(distance*direction)))
        ell=np.array([path]+priv)
        u=np.array([nmul(pa,ninv(la)) for pa,la in zip(P,priv)])
        w=np.array([nmul(ninv(la),pa) for pa,la in zip(Pp,priv)])
        ells.append(ell);us.append(u);ws.append(w)
    return np.array(ells),np.array(us),np.array(ws),radius,eps,perturb

started=time.monotonic()
devices=[str(x) for x in jax.devices()]
if args.device=='gpu':assert all(x.platform=='gpu' for x in jax.devices())
else:assert all(x.platform=='cpu' for x in jax.devices())
summary={'classification':'floating computation; not a transport theorem','device':args.device,
         'devices':devices,'jax_version':jax.__version__,'numpy_version':np.__version__,
         'seed':args.seed,'batch':args.batch,'requested_batches':args.batches,'samples':0,'candidates':[]}
compile_seconds=None;scan_seconds=0.;fd_error=None;full_min=1e99
bins={}
for index in range(args.batches):
    host=generate(args.batch,args.seed+index,args.reference)
    ell,u,w,r,e,pert=host
    before=time.monotonic()
    result,J,c=batch_diagnostic(j.asarray(ell),j.asarray(u),j.asarray(w),j.asarray(r))
    result=np.array(result);J=np.array(J);c=np.array(c)
    duration=time.monotonic()-before
    if index==0:compile_seconds=duration
    else:scan_seconds+=duration
    assert np.isfinite(result).all() and np.isfinite(J).all()
    smooth=(result[:,4]>1e-8)&(np.abs(result[:,3]-1)>1e-8)&(np.abs(result[:,3]-.5)>1e-8)
    if smooth.any():
        full_min=min(full_min,float(result[smooth,0].min()))
        assert full_min>.0309,('full-submersion bound contradicted; stop and audit',full_min)
    if args.reference and index==0:
        errs=[]
        for k in np.flatnonzero(smooth)[:8]:
            delta=np.random.default_rng(args.seed+int(k)).normal(size=(8,3))
            delta/=np.linalg.norm(delta);h=1e-6
            plus=np.array([nmul(la,nexp(h*z)) for la,z in zip(ell[k],delta)])
            minus=np.array([nmul(la,nexp(-h*z)) for la,z in zip(ell[k],delta)])
            raw=(nvalue(plus,u[k],w[k])-nvalue(minus,u[k],w[k]))/(2*h)
            fd=SQ2*nmul(ninv(c[k]),raw)[1:]
            errs.append(float(np.linalg.norm(fd-J[k]@delta.reshape(-1))))
        assert len(errs)>=4 and max(errs)<2e-6,errs
        fd_error=max(errs)
    for k in np.flatnonzero(smooth):
        key=str(int(math.floor(math.log10(e[k]))))
        row=bins.setdefault(key,{'samples':0,'min_interior_sigma':1e99,'max_interior_sigma':0.,'min_full_sigma':1e99})
        row['samples']+=1;row['min_interior_sigma']=min(row['min_interior_sigma'],float(result[k,2]))
        row['max_interior_sigma']=max(row['max_interior_sigma'],float(result[k,2]))
        row['min_full_sigma']=min(row['min_full_sigma'],float(result[k,0]))
    # Retain bounded, reproducible candidates, not full device arrays.
    for k in np.argsort(result[:,2])[:2]:
        if not smooth[k]:continue
        summary['candidates'].append({'batch_index':index,'sample_index':int(k),'epsilon':float(e[k]),
          'radius':float(r[k]),'nonabelian_perturbation':float(pert[k]),'metrics':result[k].tolist(),
          'ell':ell[k].tolist(),'u':u[k].tolist(),'w':w[k].tolist(),'c':c[k].tolist()})
    summary['candidates']=sorted(summary['candidates'],key=lambda x:x['metrics'][2])[:8]
    summary['samples']+=args.batch
    print(json.dumps({'batch_index':index,'samples':summary['samples'],'seconds':duration,
                      'min_full_sigma':full_min,'min_interior_sigma':min(x['metrics'][2] for x in summary['candidates'])}),flush=True)
    # Saturation: several independently seeded batches and all epsilon decades covered.
    if not args.reference and index>=7 and len(bins)>=6 and summary['candidates'][0]['metrics'][2]<1e-6:
        summary['early_stop']='candidate saturation after eight independent batches';break
summary.update(compile_first_batch_seconds=compile_seconds,scan_seconds=scan_seconds,
               total_seconds=time.monotonic()-started,finite_difference_error=fd_error,
               min_full_sigma=full_min,epsilon_decades=bins)
Path('result.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='candidates'},indent=2,sort_keys=True))
