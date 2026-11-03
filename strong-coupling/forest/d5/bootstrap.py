"""Window bootstrap for CALM_LINK.md section 3 (order of magnitude, SU(2)).

Window: cube of side 2R+1 blocks around the link j; m = (2R+1)^4 blocks.
  |L1| <= 81 m free links, n_P <= 576 m plaquettes (ISLANDS.md section 1).
Boundary links of the window: N_d(R) <= 8 (2R+1)^3 * 81.
Response of a window expectation to one boundary link at distance R:
  C e^{-alpha_G R}, C = max(C_*, C_b) = 9.2e13 (OpenGluon CALM_DATA.md).
Bootstrap (ANNEALED.md Lemma 3) at rate m' = alpha_G/2:
  N_d(R) C e^{-alpha_G R/2} <= 1/2.
Floor of the bootstrap: 2 eta, eta = lam |L1| (eps_loc + p_nc) (osc factors O(1) for local G).
  (i)  eps_loc small:  lam |L1| eps_loc <= 1e-3, lam = beta/4, K = 4.8e-18 (as in islands.py).
  (ii) p_nc from torus chessboard rarity (RARE.md): p_nc <= n_P q(rho_c), rho_c = K/(8 n_P).
  (iii) p_nc from the density bound localcalm (any exterior, but only a fraction of positions):
       fraction 4096 m eps, eps >= 32 (6 log beta + 31)/(beta rho^2), rho = K/(3072 m).
"""
from mpmath import mp, mpf, log, exp, e
mp.dps = 40
aG = mpf('2.5e-4'); C = mpf('9.2e13'); K = mpf('4.8e-18')

def bisect(f, lo, hi, n=300):
    # smallest x in [lo,hi] with f(x) <= 0, f decreasing
    for _ in range(n):
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return hi

R = bisect(lambda R: log(8*(2*R+1)**3*81*C) - aG*R/2 - log(mpf(1)/2), mpf(1), mpf(1e8))
R = mp.ceil(R)
m = (2*R+1)**4; L1 = 81*m; nP = 576*m
print(f'bootstrap radius R = {mp.nstr(R,4)} blocks; window m = {mp.nstr(m,3)} blocks, |L1| <= {mp.nstr(L1,3)}')
print(f'  contraction factor at R: {mp.nstr(8*(2*R+1)**3*81*C*exp(-aG*R/2),4)}')

# (i) eps_loc
A = 73*m*log(mpf('3e8')*nP**mpf('1.5')/K**3)
lb = bisect(lambda lb: A + log(exp(lb)/4*L1*1000) - mpf(3)/4*exp(lb)/4*K**2, mpf(1), mpf(2000))
print(f'(i)  lam|L1|eps_loc <= 1e-3 needs beta >= {mp.nstr(exp(lb),3)}')

# (ii) torus rarity
rc = K/(8*nP)
logq = lambda lb: -exp(lb)*rc**2/16 + 6*lb + mpf('24.7') + log(6*e*80)
lb2 = bisect(lambda lb: lb - log(4) + log(L1) + log(nP) + logq(lb) - log(mpf('1e-3')), mpf(1), mpf(2000))
print(f'(ii) torus: lam|L1| n_P q(rho_c) <= 1e-3 needs beta >= {mp.nstr(exp(lb2),3)}  (rho_c = {mp.nstr(rc,3)})')

# (iii) density bound: lam |L1| * 4096 m eps, with the smallest admissible eps
rho = K/(3072*m)
def dens(beta):
    eps = 32*(6*log(beta)+31)/(beta*rho**2)
    return beta/4*L1*4096*m*eps
for beta in [mpf('1e10'), mpf('1e100'), mpf('1e1000')]:
    print(f'(iii) density bound: lam|L1| p_nc >= {mp.nstr(dens(beta),3)} at beta = {mp.nstr(beta,3)}')
print('      (beta-independent up to log beta: the density bound can never make the floor small)')
