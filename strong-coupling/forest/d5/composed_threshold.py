"""
Composed island threshold at scale j (SU(2), Wilson measure on the torus).

Inputs are the constants of main.tex:
  C_0 = 1/2 + (2/3) log(24 pi^3)                     (Lemma wilsonZ)
  q_rho(beta) = exp(6(log beta + C_0) - beta rho^2/16) (Theorem defects)
  eta_j <= 1.023 * 4^j * eta_0, eta_0 = 1.01 rho_0    (Theorem defectscales(iv), text after it)
  spreading 3(a_j - a)                                 (Theorem defectscales(ii))
  witness box at the block scale: corners in 2y+[-5,5]^4 (proof of Theorem islandthreshold)
  Delta = (4s+3)^4 - 1, log q* = -H(log 2 + 2 log Delta) - log(e M H)

Output: beta*_j for the composed labels at scale j (blocks of side 2^{j+1} a),
with M_j, H_j computed below, and the comparison with k_* of Corollary defectscales.
"""
from mpmath import mp, mpf, log, exp, e, pi, findroot, floor
mp.dps = 40

C0 = mpf(1)/2 + mpf(2)/3*log(24*pi**3)
C1 = C0 + log(12800)/6

def witness_counts(j):
    """M_j = number of fine unit cubes that can hold the fine witness of a bad
    scale-j block, H_j = max number of scale-j blocks whose witness set contains
    a given fine cube.  Derivation: Lemma scalewitness."""
    if j == 0:
        return 11**4, 6**4
    side = 17*2**j - 4            # corners x - 2^{j+1} y in [-8 2^j + 2, 9 2^j - 3]
    lo, hi = -8*2**j + 2, 9*2**j - 3
    # number of integers y_c with 2^{j+1} y_c in [x_c - hi, x_c - lo]: brute force over x_c
    step = 2**(j+1)
    hmax = 0
    for xc in range(step):
        cnt = sum(1 for yc in range(-20, 21) if xc - hi <= step*yc <= xc - lo)
        hmax = max(hmax, cnt)
    return side**4, hmax**4

def logqstar(M, H, s):
    Delta = (4*s+3)**4 - 1
    return -H*(log(2) + 2*log(Delta)) - log(e*M*H)

def beta_star(rho0, lq):
    f = lambda b: 6*(log(b) + C0) - b*rho0**2/16 - lq
    bmax = 96/rho0**2
    assert f(bmax) > 0
    # bisection on [bmax, big]
    lo, hi = bmax, bmax
    while f(hi) > 0:
        hi *= 2
    for _ in range(200):
        mid = (lo+hi)/2
        if f(mid) > 0: lo = mid
        else: hi = mid
    return hi

def kstar(beta):
    """Corollary defectscales: largest k with 4^k * 1.01 rho_c <= 1e-4, rho_c^2 = 96(log beta + C_1)/beta."""
    rhoc = (96*(log(beta)+C1)/beta)**0.5
    k = 0
    while 4**(k+1)*mpf('1.01')*rhoc <= mpf('1e-4'):
        k += 1
    return k, rhoc

def composed(j, s, rho=mpf('1e-8')):
    M, H = witness_counts(j)
    lq = logqstar(M, H, s)
    rho0 = rho/(mpf('1.01')*mpf('1.023')*4**j) if j >= 1 else rho
    b = beta_star(rho0, lq)
    return M, H, lq, rho0, b

if __name__ == '__main__':
    print("check block scale (paper): s=1,5 at rho=1e-8")
    for s in (1,5):
        M,H,lq,rho0,b = composed(0,s)
        print(f"  s={s}: M={M} H={H} log q*={float(lq):.1f} beta*={float(b):.4e}")
    print()
    for s in (5,1):
        print(f"s = {s}, rho = 1e-8 (label radii r=4.3e-7, eta=9.2e-7)")
        print(" j |  side  |     M_j     |  H_j |  |log q*_j| |   rho_0    |   beta*_j   | 16^-j beta*_j | k_*(beta*_j) | k_*-j")
        for j in range(0, 13):
            M,H,lq,rho0,b = composed(j,s)
            k,_ = kstar(b)
            side = 11 if j==0 else 17*2**j-4
            print(f"{j:2d} | {side:6d} | {M:11d} | {H:4d} | {float(-lq):10.1f} | {float(rho0):.3e} | {float(b):.4e} | {float(b/16**j):.4e} | {k:3d} | {k-j:3d}")
        print()
    print("j_max(beta): largest j with beta >= beta*_j (s=5), and k_*(beta)")
    for ex in range(22, 44, 2):
        beta = mpf(10)**ex
        jm = -1
        for j in range(0, 30):
            if composed(j,5)[4] <= beta: jm = j
            else: break
        k,rhoc = kstar(beta)
        aj = 2**jm if jm>=0 else 0
        print(f"  beta=1e{ex}: j_max={jm:2d}  a_jmax/a={aj:8d}  k_*={k:2d}  a_k*/a={2**k:8d}  ratio 2^(k_*-j_max)={2**(k-jm) if jm>=0 else 'n/a'}   (beta/2.9e22)^(1/4)={float((beta/mpf('2.9e22'))**0.25):.3g}")
