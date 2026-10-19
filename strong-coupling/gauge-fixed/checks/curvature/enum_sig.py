import itertools
sites=list(itertools.product((0,1),repeat=4))
code=lambda r: r[0]+2*r[1]+4*r[2]+8*r[3]
def add(r,mu): s=list(r); s[mu]^=1; return tuple(s)
plaq=[]
for x in sites:
    for mu,nu in itertools.combinations(range(4),2):
        plaq.append([(x,mu),(add(x,mu),nu),(add(x,nu),mu),(x,nu)])
order=sorted(sites,key=code)
sols=[]
sig={}
def ok():
    for p in plaq:
        cnt=0; und=False
        for (s,m) in p:
            if s not in sig: und=True
            elif sig[s]==m: cnt+=1
        if cnt>1 or (not und and cnt!=1): return False
    return True
def dfs(i):
    if i==16: sols.append(tuple(sig[s] for s in order)); return
    for m in range(4):
        sig[order[i]]=m
        if ok(): dfs(i+1)
        del sig[order[i]]
dfs(0)
print(len(sols))
for s in sols[:30]: print(list(s))
