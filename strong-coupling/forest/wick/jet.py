"""CPU-only explicit mixed Taylor jets, avoiding large AD intermediates.

Same finite-box Wick formula; an independent derivative implementation.
C_s=2 C M C, C_st=4 C M C N C+4 C N C M C.
"""
from wick import args,geometry,SIGNS
import json
import time
from pathlib import Path
import numpy as np
assert args.device=='cpu'


def quartic(C,ids):
    full=np.pad(C,((0,1),(0,1)))
    B=full[ids[:,:,None],ids[:,None,:]]*SIGNS[None,:,None]*SIGNS[None,None,:]/2
    d=np.diagonal(B,axis1=1,axis2=2)
    q=np.sum(15*d*d/24,axis=1)
    for a in range(4):
        for b in range(a+1,4):
            q+=(9*d[:,a]*d[:,b]+6*B[:,a,b]**2)/4+15*(d[:,a]+d[:,b])*B[:,a,b]/6
            for k in range(4):
                if k not in (a,b):q+=(9*B[:,a,b]*d[:,k]+6*B[:,a,k]*B[:,b,k])/2
    def dd(a,b,c,e):return 9*B[:,a,b]*B[:,c,e]+3*(B[:,a,c]*B[:,b,e]+B[:,a,e]*B[:,b,c])
    return -4*(q+dd(0,1,2,3)-dd(0,2,1,3)+dd(0,3,1,2))


def transform(Y,A,B,C):
    t=np.einsum('ia,abc->ibc',A,Y,optimize=True)
    t=np.einsum('jb,ibc->ijc',B,t,optimize=True)
    return np.einsum('kc,ijc->ijk',C,t,optimize=True)


def calculate(g,target,ref):
    start=time.monotonic();C=np.linalg.inv(2*g['H']);M=g['M'][target];N=g['M'][ref]
    Cs=2*C@M@C;Ct=2*C@N@C
    Cst=4*(C@M@C@N@C+C@N@C@M@C)
    B=-g['K'];Bs=g['Ks'][target];Bt=g['Ks'][ref]
    def pairing(X,A,D,E,Y):return float(np.vdot(X,transform(Y,A,D,E)).real)
    cubic_measure=18*(3*pairing(B,Cst,C,C,B)+6*pairing(B,Cs,Ct,C,B))
    cubic_cross=18*(2*pairing(Bs,C,C,C,Bt)+6*pairing(Bs,Ct,C,C,B)+6*pairing(Bt,Cs,C,C,B))
    def bilin(X,Y):return (quartic(X+Y,g['indices'])-quartic(X,g['indices'])-quartic(Y,g['indices']))/2
    qst=2*bilin(Cst,C)+2*bilin(Cs,Ct)
    qobs=2*np.sum(bilin(Ct,C)[g['labels']==target])+2*np.sum(bilin(Cs,C)[g['labels']==ref])
    measure=-np.sum(qst)-.5*np.trace(Cst)+cubic_measure
    total=measure+qobs+cubic_cross
    return {'target':target,'reference':ref,'C0':float(6*np.trace(M@C@N@C)),
            'C2_measure':float(measure),'C2_total':float(total),
            'C2_observable_and_cross':float(qobs+cubic_cross),'seconds':time.monotonic()-start}

start=time.monotonic();reference=json.loads(Path(args.cpu_reference).read_text())
report={'classification':'independent explicit-jet floating CPU check, not interval certificate',
        'derivative_method':'analytic mixed covariance jets; no automatic Hessian','tree_order':args.tree_order,'geometries':[]}
for text in args.geometries.split(','):
    lt,ls=map(int,text.split('x'));g=geometry(lt,ls);ref=lt//2
    rows=[calculate(g,ref,ref)]
    for target in sorted(set([0,ref-1])):
        if target>=0 and target!=ref:rows.append(calculate(g,target,ref))
    expected=next(z for z in reference['geometries'] if z['lt']==lt and z['ls']==ls)
    errors=[]
    for row in rows:
        old=next(z for z in expected['correlations'] if z['target']==row['target'])
        for field in ['C0','C2_measure','C2_total']:
            errors.append(abs(row[field]-old[field])/max(1,abs(old[field])))
        row['R0']=row['C0']/rows[0]['C0']
        for part in ['measure','total']:
            row['R2_'+part]=(row['C2_'+part]-row['R0']*rows[0]['C2_'+part])/rows[0]['C0']
    assert max(errors)<1e-8,('Explicit-jet / GPU disagreement',lt,ls,errors)
    entry={'lt':lt,'ls':ls,'independent_edges':g['n'],'max_scaled_error':max(errors),'correlations':rows}
    report['geometries'].append(entry)
    Path('result.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(entry,sort_keys=True),flush=True)
report['total_seconds']=time.monotonic()-start
Path('result.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps({'total_seconds':report['total_seconds']}))
