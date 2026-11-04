# Brute-force cross-check of the hand counts 1441, 2310, 216 of patches/A05_counts.tex (not used by the proofs).
import itertools as it
planes=list(it.combinations(range(4),2))
def has_vertex_in_box(c,s,lo,hi):
    return any(all(lo[i]<=c[i]+(1 if i in s and i in d else 0)<=hi[i] for i in range(4)) for d in [set(),{s[0]},{s[1]},set(s)])
def NP_box(z,pi):
    lo=[2*z[i]-2 for i in range(4)]; hi=[2*z[i]+(5 if i in pi else 3) for i in range(4)]; return lo,hi
# N_ov: plaquettes with corner in {0,1}^4 (all parity classes), count P containing them
best=0
for c in it.product(range(2),repeat=4):
    for s in planes:
        n=0
        for pi in planes:
            for z in it.product(range(-3,4),repeat=4):
                lo,hi=NP_box(z,pi)
                if has_vertex_in_box(c,s,lo,hi): n+=1
        best=max(best,n)
print("N_ov",best)
# witness regions meeting Pi(0): plaquettes owned by block 0 have corner in {0,1}^4
cnt=0
for pi in planes:
    for z in it.product(range(-3,4),repeat=4):
        lo,hi=NP_box(z,pi)
        if any(has_vertex_in_box(c,s,lo,hi) for c in it.product(range(2),repeat=4) for s in planes): cnt+=1
print("regions meeting Pi(y)",cnt)
n=sum(1 for s in planes for c in it.product(range(-2,3),repeat=4) if has_vertex_in_box(c,s,[0]*4,[1]*4))
print("plaquettes with a vertex in one block",n)
