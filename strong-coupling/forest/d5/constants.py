"""Constants for D5.md (exact sums, mpmath for the window).

N(2a) = sum over links e of Z^4 of exp(-2 a d(e, j)), d = l-infinity distance
between link midpoints, j a fixed link.  Links at midpoint distance in [k, k+1)
are counted by shells of their base point; we use the upper bound obtained by
giving every link with base at l-inf distance k from base(j) the distance k-1
(as in the proof of Theorem constrainedcube), so the result is an upper bound.
"""
from mpmath import mp, mpf, exp, log, sqrt
mp.dps = 30
a = mpf('2.5e-4')
# links with base at l-inf distance k>=1: 4((2k+1)^4-(2k-1)^4) = 4(64k^3+16k)
def N(alpha):
    s = mpf(4)  # k = 0: four links at the base of j, distance <= 1 -> weight <= 1
    k = 1
    while True:
        t = 4*(64*k**3+16*k)*exp(-alpha*(k-1))
        s += t
        if k > 10 and t < s*mpf('1e-25'):
            break
        k += 1
    return s
N2a = N(2*a)
print('N(2 alpha_G) <=', mp.nstr(N2a, 6), ' sqrt =', mp.nstr(sqrt(N2a), 6))
pbar = mpf('4.8')
# |d/dt E F| <= (120/lam) * (2 pbar sqrt(N)) * (pbar e^a sqrt(18) cQ lam)
for cQ in [mpf(4), mpf(8)]:
    C = 120*2*pbar*sqrt(N2a)*pbar*exp(a)*sqrt(18)*cQ
    print('c_Q =', cQ, ' C_G <=', mp.nstr(C, 4), '  n >= 17 C_G/alpha: 17C_G =', mp.nstr(17*C, 4))
# window: need lam*|V|*eps_W <= 1 with eps_W from Prop calmfibre, K = 4.8e-18
K = mpf('4.8e-18')
for alpha in [mpf('0.5'), mpf('0.9')]:
    C = 120*2*pbar*sqrt(N2a)*pbar*exp(a)*sqrt(18)*8
    n = 17*C/alpha
    W = (n/2)**4
    nP = 384*W
    V = 4*n**4
    # log eps_W = 45 W log(3e8 nP^1.5/K^3) - 3 lam K^2/4 ; require <= -log(lam V)
    L = 45*W*log(mpf('3e8')*nP**mpf('1.5')/K**3)
    lam = 4*(L + 400)/(3*K**2)   # 400 >> log(lam V) here
    print('alpha =', alpha, ' n ~', mp.nstr(n, 3), ' |W| ~', mp.nstr(W, 3),
          ' beta_W = 4 lam >~', mp.nstr(4*lam, 3))
    # Prop calmrare(ii): calm exterior probability <= 7e-32 beta^1.5 |W|^-3
    print('    P(calm exterior) <=', mp.nstr(mpf('7e-32')*(4*lam)**mpf('1.5')/W**3, 3))
