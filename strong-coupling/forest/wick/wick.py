"""Finite open-box SU(2) weak-coupling coefficients by exact Gaussian Wick contractions.

Computation, not interval certification. Full maximal-tree gauge, beta_W=4/g^2,
Q/g^2 action, and the Haar Jacobian are retained. Observables are unflowed
spatial plaquette energy / g^2. This is not an infinite-volume vacuum claim.
"""
import argparse
import os
p=argparse.ArgumentParser();p.add_argument('--device',choices=['cpu','gpu'],required=True)
p.add_argument('--geometries',default='3x2');p.add_argument('--reference',action='store_true')
p.add_argument('--cpu-reference');p.add_argument('--tree-order',choices=['first','last'],default='first');p.add_argument('--compare-total-only',action='store_true');args=p.parse_args()
os.environ['JAX_PLATFORMS']='cpu' if args.device=='cpu' else 'cuda'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE']='false'
import itertools as it
import json
import math
import time
from pathlib import Path
import numpy as np
import jax
jax.config.update('jax_enable_x64',True)
import jax.numpy as j

SIGNS=np.array([1.,1.,-1.,-1.])

def geometry(lt,ls):
    shape=(lt,ls,ls,ls)
    vertices=list(it.product(*(range(s) for s in shape)))
    # Each non-origin vertex has a parent obtained by lowering its first positive coordinate.
    tree=set()
    for v in vertices[1:]:
        choices=[k for k in range(4) if v[k]]
        mu=choices[0] if args.tree_order=='first' else choices[-1]
        parent=list(v);parent[mu]-=1
        tree.add((tuple(parent),mu))
    edges=[(v,mu) for v in vertices for mu in range(4) if v[mu]+1<shape[mu] and (v,mu) not in tree]
    lookup={e:k for k,e in enumerate(edges)}
    plaquettes=[];labels=[];n=len(edges)
    for v in vertices:
        for mu,nu in it.combinations(range(4),2):
            if v[mu]+1>=shape[mu] or v[nu]+1>=shape[nu]:continue
            vm=list(v);vm[mu]+=1
            vn=list(v);vn[nu]+=1
            es=[(v,mu),(tuple(vm),nu),(tuple(vn),mu),(v,nu)]
            plaquettes.append([lookup.get(e,n) for e in es])
            labels.append(v[0] if mu>0 and nu>0 else -1)
    indices=np.array(plaquettes,dtype=np.int64);labels=np.array(labels,dtype=np.int64)
    D=np.zeros((len(indices),n))
    for row,es in enumerate(indices):
        for slot,e in enumerate(es):
            if e<n:D[row,e]+=SIGNS[slot]
    H=D.T@D
    eig=np.linalg.eigvalsh(H)
    assert eig[0]>0,('tree gauge quadratic form is not positive',eig[0])
    mats=[D[labels==t].T@D[labels==t] for t in range(lt)]
    K=np.zeros((n,n,n))
    Ks=[np.zeros_like(K) for _ in range(lt)]
    for es,t in zip(indices,labels):
        for slots in it.combinations(range(4),3):
            if any(es[s]==n for s in slots):continue
            coefficient=math.sqrt(2)/6*np.prod(SIGNS[list(slots)])
            for perm in it.permutations(range(3)):
                inversions=sum(perm[a]>perm[b] for a in range(3) for b in range(a+1,3))
                ix=tuple(es[slots[k]] for k in perm)
                val=coefficient*(-1)**inversions
                K[ix]+=val
                if t>=0:Ks[t][ix]+=val
    return {'lt':lt,'ls':ls,'n':n,'H':H,'K':K,'Ks':Ks,'M':mats,
            'indices':indices,'labels':labels,'condition':float(eig[-1]/eig[0])}


def quartic_expectation(C,indices):
    # B_ij is the covariance of oriented v_i=sign_i*A_i/sqrt(2), per colour.
    full=j.pad(C,((0,1),(0,1)))
    B=full[indices[:,:,None],indices[:,None,:]]*j.asarray(SIGNS)[None,:,None]*j.asarray(SIGNS)[None,None,:]/2
    diag=j.diagonal(B,axis1=1,axis2=2)
    scalar=j.sum(15*diag*diag/24,axis=1)
    for a,b in it.combinations(range(4),2):
        scalar+=(9*diag[:,a]*diag[:,b]+6*B[:,a,b]**2)/4
        scalar+=15*(diag[:,a]+diag[:,b])*B[:,a,b]/6
        for k in range(4):
            if k not in (a,b):scalar+=(9*B[:,a,b]*diag[:,k]+6*B[:,a,k]*B[:,b,k])/2
    def dd(a,b,c,d):return 9*B[:,a,b]*B[:,c,d]+3*(B[:,a,c]*B[:,b,d]+B[:,a,d]*B[:,b,c])
    scalar+=dd(0,1,2,3)-dd(0,2,1,3)+dd(0,3,1,2)
    return -4*scalar


def cubic_pair(K,C):
    # E[(6 K_ijk A_i^0 A_j^1 A_k^2)^2]/2 = 18 <K,C^tensor3 K>.
    T=j.einsum('ia,abc->ibc',C,K)
    T=j.einsum('jb,ibc->ijc',C,T)
    T=j.einsum('kc,ijc->ijk',C,T)
    return 18*j.sum(K*T)


def calculations(g,target,reference):
    H=j.asarray(g['H']);K=j.asarray(g['K']);Kt=j.asarray(g['Ks'][target]);Kr=j.asarray(g['Ks'][reference])
    Mt=j.asarray(g['M'][target]);Mr=j.asarray(g['M'][reference])
    ids=j.asarray(g['indices']);labels=j.asarray(g['labels'])
    def parts(st,H,K,Kt,Kr,Mt,Mr,ids,labels):
        s,t=st
        C=j.linalg.inv(2*(H-s*Mt-t*Mr))
        q4=quartic_expectation(C,ids)
        action=-j.sum(q4)-.5*j.trace(C)  # Haar log j=-g^2 sum|A|^2/6.
        observables=s*j.sum(j.where(labels==target,q4,0.))+t*j.sum(j.where(labels==reference,q4,0.))
        cubic=cubic_pair(-K+s*Kt+t*Kr,C)
        # For the measure-only part, keep the action cubic independent of source.
        measure_cubic=cubic_pair(-K,C)
        return j.array([action+measure_cubic,observables+cubic-measure_cubic,action+observables+cubic])
    def curvature(st,H,K,Kt,Kr,Mt,Mr,ids,labels):
        return jax.jacfwd(jax.jacrev(parts,argnums=0),argnums=0)(st,H,K,Kt,Kr,Mt,Mr,ids,labels)
    compute=jax.jit(curvature)
    inputs=(j.zeros(2),H,K,Kt,Kr,Mt,Mr,ids,labels)
    before=time.monotonic();second=np.array(compute(*inputs));seconds=time.monotonic()-before
    warm_seconds=None
    if args.reference:
        before=time.monotonic();np.array(compute(*inputs));warm_seconds=time.monotonic()-before
    C=np.linalg.inv(2*g['H'])
    free=6*np.trace(g['M'][target]@C@g['M'][reference]@C)
    coefficients=second[:,0,1]
    assert np.isfinite(coefficients).all()
    assert abs(coefficients[0]+coefficients[1]-coefficients[2])<1e-6*max(1,abs(coefficients[2]))
    return {'target':target,'reference':reference,'C0':float(free),
            'C2_measure':float(coefficients[0]),'C2_observable_and_cross':float(coefficients[1]),
            'C2_total':float(coefficients[2]),'seconds_including_compile':seconds,'warm_kernel_seconds':warm_seconds}


def numerical_reference(g):
    """Independent quaternion-polynomial Q3/Q4 and Gaussian moment checks."""
    rng=np.random.default_rng(20261103)
    def qm(a,b):return np.r_[a[0]*b[0]-a[1:]@b[1:],a[0]*b[1:]+b[0]*a[1:]+np.cross(a[1:],b[1:])]
    def coefficients(a):
        total=np.zeros((5,4));total[0,0]=1
        for z,sign in zip(a,SIGNS):
            norm=z@z
            q=np.zeros((5,4));q[0,0]=1;q[1,1:]=sign*z/math.sqrt(2)
            q[2,0]=-norm/4;q[3,1:]=-sign*norm*z/(12*math.sqrt(2));q[4,0]=norm*norm/96
            out=np.zeros_like(total)
            for n in range(5):
                for k in range(n+1):out[n]+=qm(total[k],q[n-k])
            total=out
        return -4*total[3,0],-4*total[4,0]
    A=rng.normal(size=(32,g['n'],3))
    q3=[];q4=[];q4_direct=[]
    for field in A:
        padded=np.vstack((field,np.zeros(3)));s3=s4=0.
        for es in g['indices']:
            c3,c4=coefficients(padded[es]);s3+=c3;s4+=c4
        q3.append(s3);q4.append(s4)
        v=padded[g['indices']]*SIGNS[None,:,None]/math.sqrt(2)
        norms=np.sum(v*v,axis=2);dots=np.einsum('pia,pja->pij',v,v)
        scalar=np.sum(norms*norms/24,axis=1)
        for i,k in it.combinations(range(4),2):
            scalar+=norms[:,i]*norms[:,k]/4+(norms[:,i]+norms[:,k])*dots[:,i,k]/6
            for z in range(4):
                if z not in (i,k):scalar+=dots[:,i,k]*norms[:,z]/2
        scalar+=dots[:,0,1]*dots[:,2,3]-dots[:,0,2]*dots[:,1,3]+dots[:,0,3]*dots[:,1,2]
        q4_direct.append(float(-4*np.sum(scalar)))
    predicted=6*np.einsum('ijk,si,sj,sk->s',g['K'],A[:,:,0],A[:,:,1],A[:,:,2],optimize=True)
    err=float(np.max(np.abs(predicted-q3)))
    assert err<1e-10,err
    err4=float(np.max(np.abs(np.array(q4_direct)-q4)))
    assert err4<1e-10,err4
    # Independent exact Gaussian Wick enumeration on one plaquette (12 variables).
    cov=rng.normal(size=(4,4));cov=cov@cov.T/10
    # Finite-difference of the quartic polynomial's expectation via exact cubature:
    # product 3-point Gauss-Hermite integrates all total-degree <=4 polynomials.
    # 3^12 is avoidable: fourth derivatives at zero are only 12^4 small entries.
    # Use a reproducible MC diagnostic with its measured standard error instead.
    size=32768
    Z=rng.normal(size=(size,4,3));L=np.linalg.cholesky(cov)
    X=np.einsum('ij,sjc->sic',L,Z)
    v=X*SIGNS[None,:,None]/math.sqrt(2)
    norm=np.sum(v*v,axis=2);dots=np.einsum('sia,sja->sij',v,v)
    scalar=np.sum(norm*norm/24,axis=1)
    for i,k in it.combinations(range(4),2):
        scalar+=norm[:,i]*norm[:,k]/4+(norm[:,i]+norm[:,k])*dots[:,i,k]/6
        for z in range(4):
            if z not in (i,k):scalar+=dots[:,i,k]*norm[:,z]/2
    scalar+=dots[:,0,1]*dots[:,2,3]-dots[:,0,2]*dots[:,1,3]+dots[:,0,3]*dots[:,1,2]
    samples=-4*scalar
    exact=float(quartic_expectation(j.asarray(cov),j.asarray([[0,1,2,3]]))[0])
    se=float(samples.std(ddof=1)/math.sqrt(size));mc=float(samples.mean())
    assert abs(mc-exact)<6*se,(mc,exact,se)
    # Cubic Wick contraction on the full small box, also an independent MC check.
    C=np.linalg.inv(2*g['H']);L=np.linalg.cholesky(C)
    Z=rng.normal(size=(16384,g['n'],3));X=np.einsum('ij,sjc->sic',L,Z)
    cubic=6*np.einsum('ijk,si,sj,sk->s',g['K'],X[:,:,0],X[:,:,1],X[:,:,2],optimize=True)
    exact3=float(2*cubic_pair(j.asarray(g['K']),j.asarray(C)))
    se3=float((cubic*cubic).std(ddof=1)/math.sqrt(len(cubic)))
    assert abs(float(np.mean(cubic*cubic))-exact3)<6*se3
    return {'quaternion_cubic_max_error':err,'quaternion_quartic_max_error':err4,'quartic_wick_exact':exact,'quartic_mc':mc,'quartic_mc_se':se,
            'cubic_square_wick_exact':exact3,'cubic_square_mc':float(np.mean(cubic*cubic)),'cubic_square_mc_se':se3,
            'reference_scope':'quaternion algebra exact to roundoff; stochastic Wick diagnostics at six-SE tolerance'}

def main():
    start=time.monotonic()
    report={'classification':'finite-volume perturbative computation, not interval-certified and not continuum vacuum',
            'device':args.device,'devices':[str(d) for d in jax.devices()],'jax_version':jax.__version__,'tree_order':args.tree_order,'device_kinds':[d.device_kind for d in jax.devices()],'geometries':[]}
    if args.device=='gpu':assert len(jax.devices())==1 and all(d.platform=='gpu' for d in jax.devices())
    reference=json.loads(Path(args.cpu_reference).read_text()) if args.cpu_reference else None
    for text in args.geometries.split(','):
        lt,ls=map(int,text.split('x'));before=time.monotonic();g=geometry(lt,ls)
        row={'lt':lt,'ls':ls,'independent_edges':g['n'],'quadratic_condition':g['condition'],
             'geometry_seconds':time.monotonic()-before}
        if args.reference:row['reference_checks']=numerical_reference(g)
        ref=lt//2
        same=calculations(g,ref,ref)
        rows=[same]
        for target in sorted(set([0,ref-1])):
            if target>=0 and target!=ref:rows.append(calculations(g,target,ref))
        for r in rows:
            r['R0']=r['C0']/same['C0']
            for part in ['measure','total']:
                r['R2_'+part]=(r['C2_'+part]-r['R0']*same['C2_'+part])/same['C0']
        row['correlations']=rows
        if reference:
            matching=[v for v in reference['geometries'] if v['lt']==lt and v['ls']==ls]
            if matching:
                errors=[]
                for actual in rows:
                    expected=next(v for v in matching[0]['correlations'] if v['target']==actual['target'])
                    for field in (['C0','C2_total'] if args.compare_total_only else ['C0','C2_measure','C2_total']):
                        errors.append(abs(actual[field]-expected[field])/max(1,abs(expected[field])))
                assert max(errors)<1e-8,('CPU/GPU reference disagreement',errors)
                row['cpu_reference_max_scaled_error']=max(errors)
        report['geometries'].append(row)
        print(json.dumps(row,sort_keys=True),flush=True)
        Path('result.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    report['total_seconds']=time.monotonic()-start
    if args.device=='gpu':report['device_memory_stats']=jax.devices()[0].memory_stats()
    Path('result.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'total_seconds':report['total_seconds'],'device':args.device}))

if __name__=='__main__':main()
