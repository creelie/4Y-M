import re,sys
from fractions import Fraction as F
L=open(sys.argv[1]).read().split('\n')
def rank(M):
    M=[r[:] for r in M];n=len(M);r=0
    for c in range(n):
        p=next((i for i in range(r,n) if M[i][c]!=0),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        for i in range(n):
            if i!=r and M[i][c]!=0:
                f=M[i][c]/M[r][c]; M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        r+=1
    return r
for i,l in enumerate(L):
    m=re.match(r'theorem (romanov_G_\d+)',l)
    if not m: continue
    s=L[i+2].split('≤',1)[1].split(':= by')[0]
    A=[[F(0)]*8 for _ in range(8)]
    for c,a,b in re.findall(r'\((-?\d+(?:/\d+)?) : ℝ\) \* x(\d) \* x(\d)',s):
        A[int(a)][int(b)]+=F(c)/2; A[int(b)][int(a)]+=F(c)/2
    v=[1,0,0,0,-1,0,0,0]
    Av=[sum(A[r][c]*v[c] for c in range(8)) for r in range(8)]
    print(m.group(1),'rank(A-1e-4 I)=',rank([[A[r][c]-(F(1,10000) if r==c else 0) for c in range(8)] for r in range(8)]),
          'rank(A)=',rank(A),'denoms:',sorted({x.denominator for row in A for x in row})[-3:])
print('--- precise lambda_min (exact rational bisection)')
def pd(A,t):
    M=[[A[i][j]-(t if i==j else 0) for j in range(8)] for i in range(8)]
    for k in range(8):
        if M[k][k]<=0: return False
        for i in range(k+1,8):
            f=M[i][k]/M[k][k]
            for j in range(k,8): M[i][j]-=f*M[k][j]
    return True
for i,l in enumerate(L):
    m=re.match(r'theorem (romanov_G_\d+)',l)
    if not m or m.group(1) not in ('romanov_G_0000','romanov_G_1011'): continue
    s=L[i+2].split('≤',1)[1].split(':= by')[0]
    A=[[F(0)]*8 for _ in range(8)]
    for c,a,b in re.findall(r'\((-?\d+(?:/\d+)?) : ℝ\) \* x(\d) \* x(\d)',s):
        A[int(a)][int(b)]+=F(c)/2; A[int(b)][int(a)]+=F(c)/2
    lo,hi=F(0),F(1,1000)
    for _ in range(40):
        mid=(lo+hi)/2
        if pd(A,mid): lo=mid
        else: hi=mid
    print(m.group(1),'lambda_min ~ %.10e'%float(lo))
