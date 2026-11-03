"""Exact rational certificate for the smallest finite-box Wick coefficients.

No numerical libraries. K=(sqrt(2)/6)T, T integer, so the Wick prefactor
18 cancels the irrational factors exactly. This is a finite computation.
"""
from fractions import Fraction as F
from itertools import product,combinations,permutations
from functools import reduce
from math import gcd
import json
import time
from pathlib import Path
S=[1,1,-1,-1]

def mm(A,B):
    Bt=list(zip(*B))
    return [[sum(a*b for a,b in zip(row,col)) for col in Bt] for row in A]
def scale(A,k):return [[k*a for a in row] for row in A]
def plus(A,B):return [[a+b for a,b in zip(x,y)] for x,y in zip(A,B)]
def inverse(A):
    n=len(A);work=[[F(a) for a in row]+[F(i==j) for j in range(n)] for i,row in enumerate(A)]
    for k in range(n):
        pivot=next(i for i in range(k,n) if work[i][k]);work[k],work[pivot]=work[pivot],work[k]
        v=work[k][k];work[k]=[x/v for x in work[k]]
        for i in range(n):
            if i!=k and work[i][k]:
                v=work[i][k];work[i]=[x-v*y for x,y in zip(work[i],work[k])]
    return [row[n:] for row in work]
def lcm(a,b):return a//gcd(a,b)*b

def build(last=False):
    vertices=list(product(range(3),range(2),range(2),range(2)));shape=[3,2,2,2];tree=set()
    for v in vertices[1:]:
        choices=[i for i in range(4) if v[i]];mu=choices[-1] if last else choices[0]
        parent=list(v);parent[mu]-=1;tree.add((tuple(parent),mu))
    edges=[(v,m) for v in vertices for m in range(4) if v[m]+1<shape[m] and (v,m) not in tree]
    n=len(edges);lookup={e:i for i,e in enumerate(edges)};plaq=[];labels=[];D=[]
    T={};Ts=[{} for _ in range(3)]
    for v in vertices:
        for mu,nu in combinations(range(4),2):
            if v[mu]+1>=shape[mu] or v[nu]+1>=shape[nu]:continue
            a=list(v);a[mu]+=1;b=list(v);b[nu]+=1
            es=[lookup.get(e,n) for e in [(v,mu),(tuple(a),nu),(tuple(b),mu),(v,nu)]]
            label=v[0] if mu>0 else -1;plaq.append(es);labels.append(label)
            row=[0]*n
            for i,e in enumerate(es):
                if e<n:row[e]+=S[i]
            D.append(row)
            for slots in combinations(range(4),3):
                if any(es[k]==n for k in slots):continue
                sign=S[slots[0]]*S[slots[1]]*S[slots[2]]
                for p in permutations(range(3)):
                    inv=sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
                    ix=tuple(es[slots[i]] for i in p);value=sign*(-1)**inv
                    T[ix]=T.get(ix,0)+value
                    if label>=0:Ts[label][ix]=Ts[label].get(ix,0)+value
    DT=[list(x) for x in zip(*D)];H=mm(DT,D);C=inverse(scale(H,2))
    denominator=reduce(lcm,(x.denominator for row in C for x in row),1)
    A=[[int(x*denominator) for x in row] for row in C]
    M=[]
    for t in range(3):
        Dt=[row for row,label in zip(D,labels) if label==t]
        M.append(mm([list(x) for x in zip(*Dt)],Dt))
    clean=lambda d:{k:v for k,v in d.items() if v}
    return n,denominator,A,M,plaq,labels,clean(T),[clean(t) for t in Ts]


def q4(C,plaq):
    n=len(C);out=[]
    for es in plaq:
        B=[[F(S[i]*S[k],2)*C[es[i]][es[k]] if es[i]<n and es[k]<n else F(0) for k in range(4)] for i in range(4)]
        d=[B[i][i] for i in range(4)]
        q=sum((F(15,24)*x*x for x in d),F(0))
        for i,k in combinations(range(4),2):
            q+=(9*d[i]*d[k]+6*B[i][k]**2)/4+F(15,6)*(d[i]+d[k])*B[i][k]
            for z in range(4):
                if z not in (i,k):q+=(9*B[i][k]*d[z]+6*B[i][z]*B[k][z])/2
        def dd(i,k,z,w):return 9*B[i][k]*B[z][w]+3*(B[i][z]*B[k][w]+B[i][w]*B[k][z])
        out.append(-4*(q+dd(0,1,2,3)-dd(0,2,1,3)+dd(0,3,1,2)))
    return out

def pair(T,U,A,B,C):
    answer=0
    for (i,j,k),v in T.items():
        for (a,b,c),w in U.items():answer+=v*w*A[i][a]*B[j][b]*C[k][c]
    return answer

def run(last):
    n,d,A,M,plaq,labels,T,Ts=build(last);print('Exact geometry',last,n,'nonzero cubic entries',len(T),flush=True)
    C=[[F(a,d) for a in row] for row in A];ref=1;results=[]
    for target in [ref,0]:
        Ms,Mt=M[target],M[ref]
        As=scale(mm(mm(A,Ms),A),2);At=scale(mm(mm(A,Mt),A),2)
        Ast=scale(plus(mm(mm(mm(mm(A,Ms),A),Mt),A),mm(mm(mm(mm(A,Mt),A),Ms),A)),4)
        Cs=[[F(a,d*d) for a in row] for row in As];Ct=[[F(a,d*d) for a in row] for row in At]
        Cst=[[F(a,d**3) for a in row] for row in Ast]
        cm=F(3*pair(T,T,Ast,A,A)+6*pair(T,T,As,At,A),d**5)
        cc=F(2*pair(Ts[target],Ts[ref],A,A,A),d**3)-F(6*(pair(Ts[target],T,At,A,A)+pair(Ts[ref],T,As,A,A)),d**4)
        def bil(X,Y):return [(a-b-c)/2 for a,b,c in zip(q4(plus(X,Y),plaq),q4(X,plaq),q4(Y,plaq))]
        qst=[2*a+2*b for a,b in zip(bil(Cst,C),bil(Cs,Ct))]
        measure=-sum(qst)-sum(Cst[i][i] for i in range(n))/2+cm
        obs=2*sum(v for v,label in zip(bil(Ct,C),labels) if label==target)+2*sum(v for v,label in zip(bil(Cs,C),labels) if label==ref)
        total=measure+obs+cc
        mat=mm(mm(mm(Ms,A),Mt),A);free=F(6*sum(mat[i][i] for i in range(n)),d*d)
        results.append({'target':target,'C0':free,'C2_measure':measure,'C2_total':total})
    for r in results:
        r['R0']=r['C0']/results[0]['C0'];r['R2_total']=(r['C2_total']-r['R0']*results[0]['C2_total'])/results[0]['C0']
    return results

start=time.monotonic();first=run(False);last=run(True)
for a,b in zip(first,last):
    assert a['C0']==b['C0'] and a['C2_total']==b['C2_total'] and a['R2_total']==b['R2_total']
assert first[1]['R2_total']>0
report={'classification':'exact rational finite-box computation, not continuum or infinite-vacuum result',
        'first_tree':[{k:str(v) for k,v in row.items()} for row in first],
        'last_tree':[{k:str(v) for k,v in row.items()} for row in last],
        'tree_invariance':'exact equality of complete coefficients','seconds':time.monotonic()-start}
Path('result.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
