"""Order-of-magnitude ledger for LARGE_FIELD.md (crude; see the note).

Bad set = links leaving a fixed geodesic ball of radius r1 in fibre gauge.
Sketch assumption: a bad link forces a plaquette with delta > rho = r1/4 in its block closure.
Rarity (RARE.md, torus): log q(rho) = -beta rho^2/16 + 6 log beta + 24.7 + log(6e*80).
Island + collar cost per block (ISLANDS.md 2): (4.5 beta C_G)^324.
Coarse scales by chessboard alone: a coarse block of side 2^k with holonomy defect r1 forces
a fine plaquette > r1/4^(k+1) (Stokes over <= 4^(k+1) plaquettes); usable while the rarity
condition still holds at that threshold.
"""
from mpmath import mp, mpf, log, e
mp.dps = 40
CG = mpf('3e13')

def island_ok(lb, rho):
    return -mp.e**lb*rho**2/16 + 6*lb + mpf('24.7') + log(6*e*80) + 324*(log(mpf('4.5'))+lb+log(CG)) + log(2) < 0

def beta_needed(rho):
    lo, hi = mpf(0), mpf(5000)
    for _ in range(300):
        mid = (lo+hi)/2
        lo, hi = (lo, mid) if island_ok(mid, rho) else (mid, hi)
    return mp.e**hi

print('basin radius r1 | plaquette level rho | beta for islands beaten by rarity')
for r1 in ['4.8e-18', '4.7e-16', '1e-8', '1e-4', '1e-2', '1e-1']:
    r1 = mpf(r1); rho = r1/4
    print(f'{mp.nstr(r1,3):>8} | {mp.nstr(rho,3):>8} | {mp.nstr(beta_needed(rho),3)}')

print()
print('chessboard-only coarse scales k (coarse blocks of side 2^k) at which the island condition still holds')
for r1, beta in [('4.8e-18', '1e38'), ('4.8e-18', '1e64'), ('1e-2', '1e8'), ('1e-2', '1e20')]:
    r1 = mpf(r1); beta = mpf(beta); k = 0
    while island_ok(log(beta), r1/4**(k+1)):
        k += 1
    print(f'r1={mp.nstr(r1,3)} beta={mp.nstr(beta,3)}: scales k = 0..{k-1} usable ({k} scales)')
