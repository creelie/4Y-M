"""Order-of-magnitude thresholds for ISLANDS.md (crude; see the note).

Window of side n0 lattice units (blocks of side 2): m = (n0/2)^4 blocks.
Local calm-fibre bound (Prop L): exponent 73 m, n_P <= 576 m, K = 4.8e-18.
Calm exterior threshold for the window: rho_c = K/(8 n_P)  (Cor. calmexterior).
Rarity (RARE.md): log q(rho) = -beta rho^2/16 + 6 log beta + 24.7 + log(6e*80).
Island + collar cost per site: (4.5 beta * C_G)^(4*81).
"""
from mpmath import mp, mpf, log, e
mp.dps = 40
K = mpf('4.8e-18'); CG = mpf('3e13')
for alpha in [mpf('0.9'), mpf('0.5')]:
    n0 = 17*CG/alpha
    m = (n0/2)**4
    nP = 576*m
    # (1) off-chart error in the window small: (3/4) lam K^2 >= 73 m log(3e8 nP^1.5/K^3) + log(lam*4 n0^4) + 10
    L1 = 73*m*log(mpf('3e8')*nP**mpf('1.5')/K**3)
    lam1 = 4*(L1 + 500)/(3*K**2); b1 = 4*lam1
    # (2) windows calm with tails: q(rho_c) * (4.5 beta CG)^(324) < 1/2
    rho = K/(8*nP)
    f = lambda lb: -mp.e**lb*rho**2/16 + 6*lb + mpf('24.7') + log(6*e*80) + 324*(log(mpf('4.5'))+lb+log(CG)) + log(2)
    lo, hi = mpf(1), mpf(2000)
    for _ in range(200):
        mid=(lo+hi)/2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    lb = hi
    print(f'alpha={alpha}: n0={mp.nstr(n0,3)}  m={mp.nstr(m,3)}  beta for small off-chart error >= {mp.nstr(b1,3)}'
          f'  rho_c={mp.nstr(rho,3)}  beta for calm-window islands >= {mp.nstr(mp.e**lb,3)}')
