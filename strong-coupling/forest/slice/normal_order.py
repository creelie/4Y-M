# Normal-ordered form of L^3, L = -d^2 + x^2 + 4 in one variable, as
# sum c_{ab} x^a d^b (x to the left), in exact integer arithmetic, and the
# constant sum_{k>=0} (2k+5)^(-7/4) used in Lemma lem:hermite (mpmath iv).
from math import comb, factorial
from collections import defaultdict
def mul(A,B):
    C=defaultdict(int)
    for (a,b),u in A.items():
        for (c,d),v in B.items():
            for j in range(0,min(b,c)+1):
                C[(a+c-j,b-j+d)]+=u*v*comb(b,j)*factorial(c)//factorial(c-j)
    return {k:v for k,v in C.items() if v}
L={(0,2):-1,(2,0):1,(0,0):4}
L3=mul(L,mul(L,L))
print(sorted(L3.items()))
lam=sum(abs(v) for v in L3.values())
print('lambda =',lam)
print('max a, max b =',max(a for a,b in L3),max(b for a,b in L3))
from mpmath import iv, mpf
iv.dps=30
# sum_{k>=0}(2k+5)^(-7/4) <= 5^(-7/4) + int_0^inf (2x+5)^(-7/4) dx = 5^(-7/4) + (2/3) 5^(-3/4)
s=iv.mpf(5)**(iv.mpf(-7)/4)+iv.mpf(2)/3*iv.mpf(5)**(iv.mpf(-3)/4)
print('sum bound =',s)
print('pi*sum bound =',iv.pi*s)
