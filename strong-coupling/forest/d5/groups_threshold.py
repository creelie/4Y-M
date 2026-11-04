"""Certified constants for the island machinery of Section 13 on a compact
simple group G in a faithful unitary representation of dimension N.

Interval arithmetic (mpmath.iv, 60 digits).  Every printed number is the
midpoint of a certified interval together with the interval's endpoints;
inequalities are checked on the interval endpoints.

Quantities (notation of groups.tex):
  d = dim G, kappa = kappa_G, K = kappa/(d-1), I_d(pi) = int_0^pi sin^{d-1},
  C0(G) = 1/2 + (2/3)[log(d I_d(pi)) + (d/2) log(16/(K N)) - (d-1) log(95/96)],
  f_G(b) = 2 d log b + 6 C0(G) - b rho^2/(8N),
  log q* = -H (log 2 + 2 log D) - log(e M H), M = 11^4, H = 6^4,
  D = (4s+3)^4 - 1 (islands) or (4s+17)^4 - 1 (joining distance 2s+8),
  beta* = largest root of f_G = log q* (- H log Gamma).
"""
from fractions import Fraction as Fr
from mpmath import iv, mp, mpf, log as mlog, findroot
import math

iv.dps = 60
mp.dps = 60

def I_dpi(d):
    """I_d(pi) = int_0^pi sin^{d-1}u du, exactly: pi * rational or rational."""
    if d % 2 == 1:   # d-1 even
        val = Fr(1)  # I_1 = pi
        k = 1
        while k < d:
            k += 2
            val *= Fr(k - 2, k - 1)
        return ('pi', val)
    val = Fr(2)      # I_2 = 2
    k = 2
    while k < d:
        k += 2
        val *= Fr(k - 2, k - 1)
    return ('1', val)

def ivq(q):
    return iv.mpf(q.numerator) / q.denominator

def I_iv(d):
    kind, val = I_dpi(d)
    return ivq(val) * (iv.pi if kind == 'pi' else 1)

def C0(d, kappa, N):
    K = ivq(Fr(kappa)) / (d - 1)
    return (iv.mpf(1) / 2 + iv.mpf(2) / 3 * (iv.log(d * I_iv(d))
            + iv.mpf(d) / 2 * iv.log(16 / (K * N))
            - (d - 1) * iv.log(iv.mpf(95) / 96)))

def C0_paper_su2():
    return iv.mpf(1) / 2 + iv.mpf(2) / 3 * iv.log(24 * iv.pi ** 3)

M, H = 11 ** 4, 6 ** 4

def logqstar(D, Mv=M, Hv=H):
    return -Hv * (iv.log(2) + 2 * iv.log(D)) - iv.log(iv.e * Mv * Hv)

def f(b, d, N, rho, c0):
    return 2 * d * iv.log(b) + 6 * c0 - b * rho ** 2 / (8 * N)

def lo_(x):
    return mpf(x.a)

def hi_(x):
    return mpf(x.b)

def mid(x):
    if not isinstance(x, iv.mpf):
        return mpf(x)
    return (mpf(x.a) + mpf(x.b)) / 2

def betastar(d, N, rho, c0, lq):
    """Certified bracket [lo, hi] of the largest root of f = lq."""
    rho_m, c0_m, lq_m = mid(iv.mpf(rho)), mid(c0), mid(lq)
    fm = lambda b: 2 * d * mlog(b) + 6 * c0_m - b * rho_m ** 2 / (8 * N) - lq_m
    b0 = 8 * N * (-lq_m) / rho_m ** 2
    b = findroot(fm, b0)
    peak = mpf(16) * d * N / rho_m ** 2
    assert b > peak
    lo, hi = b * (1 - mpf(10) ** -25), b * (1 + mpf(10) ** -25)
    flo = f(iv.mpf(lo), d, N, iv.mpf(rho), c0) - lq
    fhi = f(iv.mpf(hi), d, N, iv.mpf(rho), c0) - lq
    assert mpf(flo.a) > 0 and mpf(fhi.b) < 0, "bracket not certified"
    # f is concave with its maximum at 16 d N / rho^2 < lo, so the root in
    # [lo,hi] is the largest root; check the maximum lies above log q*.
    fpk = f(iv.mpf(peak), d, N, iv.mpf(rho), c0) - lq
    assert mpf(fpk.a) > 0
    assert peak < lo
    return lo, hi

def eps_of(lo, hi, N, rho, lq):
    # beta* = 8 N |log q*| rho^-2 (1+eps)
    base = 8 * N * (-lq) / iv.mpf(rho) ** 2
    return iv.mpf([lo, hi]) / base - 1

def show(x, digits=6):
    return mp.nstr(mid(x), digits)

def rng(x, digits=8):
    return f"[{mp.nstr(mpf(x.a), digits)}, {mp.nstr(mpf(x.b), digits)}]"

def sci(lo, hi, digits=5):
    return f"{mp.nstr(mpf(lo), digits)} .. {mp.nstr(mpf(hi), digits)}"

# ---------------------------------------------------------------------------
groups = [
    # name, d, N, kappa, rho_G
    ("SU(2)", 3, 2, Fr(1), "1"),
    ("SU(3)", 8, 3, Fr(3, 2), "1"),
    ("SU(4)", 15, 4, Fr(2), "1"),
    ("SU(5)", 24, 5, Fr(5, 2), "1"),
    ("SU(10)", 99, 10, Fr(5), "1"),
    ("SO(3)", 3, 3, Fr(1, 4), "1"),
    ("SO(5)", 10, 5, Fr(3, 4), "1"),
    ("SO(10)", 45, 10, Fr(2), "1"),
    ("Sp(2)", 10, 4, Fr(3, 2), "1"),
    ("Sp(3)", 21, 6, Fr(2), "1"),
    ("G2", 14, 7, Fr(1), "0.99"),
    ("F4", 52, 26, Fr(3, 4), "0.99"),
    ("E6", 78, 27, Fr(1), "0.99"),
    ("E7", 133, 56, Fr(3, 4), "0.76"),
    ("E8", 248, 248, Fr(1, 4), "0.99"),
]

print("== (A) Haar small-ball constant C0(G) and the hypothesis beta_W >= K N ==")
print("group   d    N   kappa   K N        I_d(pi)               C0(G)")
c0s = {}
for name, d, N, kap, rG in groups:
    c0 = C0(d, kap, N)
    c0s[name] = c0
    KN = Fr(kap) * N / (d - 1)
    assert KN <= 1
    kind, val = I_dpi(d)
    Is = f"{val}*pi" if kind == 'pi' else f"{val}"
    print(f"{name:7s} {d:3d} {N:4d}  {str(kap):5s}  {str(KN):9s}  {Is:20s}  "
          f"{rng(c0)}")
    assert mpf(c0.a) > 0
print("every K N <= 1, so Lemma (Wilson Z) applies at every beta_W >= 1;"
      " every C0(G) > 0")
c0p = C0_paper_su2()
print(f"paper's SU(2) constant 1/2+(2/3)log(24 pi^3) = {rng(c0p)};"
      f" C0(SU(2)) of the general lemma = {rng(c0s['SU(2)'])}")
assert mpf(c0s['SU(2)'].b) < mpf(c0p.a)

print()
print("== (B) log q* (group-free) ==")
lqs = {}
for s in (1, 5):
    D = (4 * s + 3) ** 4 - 1
    Ds = (4 * s + 17) ** 4 - 1
    lqs[(s, 'isl')] = logqstar(D)
    lqs[(s, 'join')] = logqstar(Ds)
    print(f"s={s}: Delta={D}, |log q*|={rng(-lqs[(s,'isl')], 9)};"
          f" N_s={Ds}, |log q*|={rng(-lqs[(s,'join')], 9)}")

print()
print("== (C) consistency: SU(2) with the paper's C0 reproduces"
      " Theorem (island threshold) ==")
for rho, s in ((mpf('1e-8'), 1), (mpf('1e-8'), 5), (mpf('1e-9'), 1),
               (mpf('1e-9'), 5)):
    lo, hi = betastar(3, 2, rho, c0p, lqs[(s, 'isl')])
    print(f"rho={mp.nstr(rho,2)}, s={s}: beta* in {sci(lo, hi, 7)}")
lo, hi = betastar(3, 2, mpf('1e-8'), c0p, lqs[(5, 'isl')] - H * iv.log(iv.mpf(10) ** 30))
print(f"rho=1e-8, s=5, Gamma=1e30: beta* in {sci(lo, hi, 7)}")

print()
print("== (D) SU(3): f(b) = 16 log b + 6 C0 - b rho^2/24 ==")
d, N = 8, 3
c0 = c0s["SU(3)"]
for rho, s in ((mpf('1e-8'), 1), (mpf('1e-8'), 5), (mpf('1e-9'), 1),
               (mpf('1e-9'), 5)):
    lq = lqs[(s, 'isl')]
    lo, hi = betastar(d, N, rho, c0, lq)
    e = eps_of(lo, hi, N, rho, lq)
    print(f"rho={mp.nstr(rho,2)}, s={s}: beta* in {sci(lo, hi, 7)};"
          f"  eps in [{mp.nstr(mpf(e.a),4)}, {mp.nstr(mpf(e.b),4)}]")
    assert 0 < mpf(e.a) and mpf(e.b) < 0.05
lq = lqs[(5, 'isl')] - H * iv.log(iv.mpf(10) ** 30)
lo, hi = betastar(d, N, mpf('1e-8'), c0, lq)
print(f"rho=1e-8, s=5, Gamma=1e30: beta* in {sci(lo, hi, 7)}")
print("ratio beta*_SU(3)/beta*_SU(2,paper C0) at rho=1e-8, s=1:")
lo3, hi3 = betastar(8, 3, mpf('1e-8'), c0, lqs[(1, 'isl')])
lo2, hi2 = betastar(3, 2, mpf('1e-8'), c0p, lqs[(1, 'isl')])
print(f"   {mp.nstr(lo3/hi2, 6)} .. {mp.nstr(hi3/lo2, 6)}")

print()
print("== (E) SU(3): island threshold at scale j (s=5, rho=1e-8) ==")
def Hj(j):
    return 8 ** 4 if j <= 2 else 9 ** 4
worst = mpf(0)
for j in range(1, 13):
    Mj = (17 * 2 ** j - 4) ** 4
    lqj = logqstar((4 * 5 + 3) ** 4 - 1, Mj, Hj(j))
    rho0 = mpf('1e-8') / (mpf('1.01') * mpf('1.023') * 4 ** j)
    # rho0 is rounded: use the interval of the exact quotient
    rho0_iv = iv.mpf('1e-8') / (iv.mpf('1.01') * iv.mpf('1.023') * 4 ** j)
    lo, hi = betastar(8, 3, mpf(rho0_iv.b), c0, lqj)
    # the root is decreasing in rho, so the upper endpoint of rho gives a
    # lower bound; recompute with rho0_iv.a for the upper bound
    lo2_, hi2_ = betastar(8, 3, mpf(rho0_iv.a), c0, lqj)
    r = hi2_ / mpf(16) ** j
    worst = max(worst, r)
    if j <= 3:
        print(f"j={j}: beta*_j in {sci(lo, hi2_, 5)}")
print(f"max_j<=12 16^-j beta*_j <= {mp.nstr(worst, 5)}")

print()
print("== (F) SU(3): two islands, layer bound and the window ==")
r3 = 3
c3 = iv.mpf(15000) * iv.sqrt(r3)        # > dim SU(3) = 8
c6 = iv.mpf(681000) * 8 * c3 ** 2
rho_max = 1 / (iv.mpf(2758) * 195 * 400000 * c3)
print(f"c3 = 1.5e4 sqrt3 in [{mp.nstr(mpf(c3.a),8)}, {mp.nstr(mpf(c3.b),8)}]")
print(f"c6 = 6.81e5*8*c3^2 in [{mp.nstr(mpf(c6.a),6)}, {mp.nstr(mpf(c6.b),6)}]")
print(f"hypothesis (a): rho <= 1/(2758*195*4e5*c3) in"
      f" [{mp.nstr(mpf(rho_max.a),6)}, {mp.nstr(mpf(rho_max.b),6)}]")
bmin = 480 * 3 * c6
print(f"hypothesis (b): beta_W >= 480 N c6 = [{mp.nstr(mpf(bmin.a),6)},"
      f" {mp.nstr(mpf(bmin.b),6)}]")
rho = mpf('1.78e-16')
assert rho < mpf(rho_max.a)
for s in (1, 5):
    lq = lqs[(s, 'join')]
    lo, hi = betastar(8, 3, rho, c0, lq)
    print(f"rho=1.78e-16, s={s}, joining distance 2s+8: beta* in {sci(lo, hi, 7)}")
    assert lo > mpf(bmin.b)
    if s == 1:
        bst = iv.mpf([lo, hi])
print("so hypothesis (b) holds at every beta_W >= beta*")
# c2 of the window theorem at s=1, theta=0, kappa=2 alpha_G
alphaG = iv.mpf('2.5e-4')
r = 43 * iv.mpf(rho)
Lam2 = iv.mpf('66.0001') * bst + iv.mpf(8) / 3 * iv.mpf('4.65e6')
s = 1
c2 = (960 * 3 / bst) * Lam2 ** 2 * r ** 2 * 45 * (2 * s + 12) ** 2 * \
     (2 * s + 7) ** 2 * iv.exp(2 * alphaG * (2 * s + 14))
print(f"Lambda2bar(beta*) <= {mp.nstr(mpf(Lam2.b), 6)}")
print(f"c2 (s=1, theta=0) <= 10^{mp.nstr(mpf(iv.log(c2).b) / math.log(10), 5)}")
# SU(2) value for comparison, paper's numbers
rho2 = iv.mpf('2.17e-16')
lo, hi = betastar(3, 2, mpf('2.17e-16'), c0p, lqs[(1, 'join')])
b2 = iv.mpf([lo, hi])
Lam22 = iv.mpf('99.0001') * b2 + iv.mpf('4.65e6')
c22 = (1920 / b2) * Lam22 ** 2 * (43 * rho2) ** 2 * 45 * 14 ** 2 * 9 ** 2 * \
      iv.exp(2 * alphaG * 16)
print(f"check SU(2): beta*={sci(lo, hi, 6)}, c2 <= 10^"
      f"{mp.nstr(mpf(iv.log(c22).b) / math.log(10), 5)} (paper: 10^22.13)")

print()
print("== (G) the threshold in the other groups (rho = 1e-8 rho_G, s = 1,5) ==")
print("group   rho_G   s=1: beta*            eps         s=5: beta*")
for name, d, N, kap, rG in groups:
    rho = mpf('1e-8') * mpf(rG)
    out = []
    for s in (1, 5):
        lo, hi = betastar(d, N, rho, c0s[name], lqs[(s, 'isl')])
        e = eps_of(lo, hi, N, rho, lqs[(s, 'isl')])
        assert mpf(e.a) > 0
        out.append((lo, hi, e))
    print(f"{name:7s} {rG:5s}  {mp.nstr(out[0][1],5):12s} eps<{mp.nstr(mpf(out[0][2].b),3):9s}"
          f" {mp.nstr(out[1][1],5)}")

print()
print("== (H) the general formula ==")
print("beta*_G = 8 N |log q*| rho^-2 (1+eps), eps = (2 d log beta*_G + 6 C0(G))/|log q*|,")
print("an identity at the root; table (G) gives the values of eps.")

print()
print("== (I) group constants of the layer bound, the mixed derivatives and the site ==")
print("c3 = max(1.5e4 sqrt N, d), c6 = 6.81e5 d c3^2; (b): beta_W >= 480 N c6;")
print("(a) at r=43 rho, eta=92 rho: rho <= 1/(2758*7.8e7*c3);")
print("Lambda2 <= (2/N) 99.0001 beta_W + (d/3) 4.65e6;")
print("site: C <= 197 beta_W + a_G (non-path), 87 beta_W + b_G (inside),")
print("  a_G = 24 sqrt N * 64.94 (d/3) sqrt(kappa) + 1.1986,")
print("  b_G = 24 sqrt N * 32.47 (d/3) sqrt(kappa) + 0.5993;  delta_G = kappa/3")
print("group    c3         (b) beta_W >=   (a) rho <=    cW coeff   cpsi       a_G       b_G")
for name, d, N, kap, rG in groups:
    c3 = max(iv.mpf(15000) * iv.sqrt(N), iv.mpf(d)) if d < 15000 * math.sqrt(N) else iv.mpf(d)
    c6 = iv.mpf(681000) * d * c3 ** 2
    bmin = 480 * N * c6
    rmax = 1 / (iv.mpf(2758) * 195 * 400000 * c3)
    cw = iv.mpf('99.0001') * 2 / N
    cp = iv.mpf('4.65e6') * d / 3
    k = ivq(Fr(kap))
    aG = 24 * iv.sqrt(N) * iv.mpf('64.94') * d / 3 * iv.sqrt(k) + iv.mpf('1.1986')
    bG = 24 * iv.sqrt(N) * iv.mpf('32.47') * d / 3 * iv.sqrt(k) + iv.mpf('0.5993')
    print(f"{name:7s} {mp.nstr(mpf(c3.b),6):10s} {mp.nstr(mpf(bmin.b),4):14s} "
          f"{mp.nstr(mpf(rmax.a),4):12s} {mp.nstr(mpf(cw.b),7):10s} "
          f"{mp.nstr(mpf(cp.b),4):10s} {mp.nstr(mpf(aG.b),6):9s} {mp.nstr(mpf(bG.b),6)}")
print("check SU(2): a_G < 2206 - 0 and b_G < 1103 (the constants of Theorem (site))")
aG2 = 24 * iv.sqrt(2) * iv.mpf('64.94') + iv.mpf('1.1986')
bG2 = 24 * iv.sqrt(2) * iv.mpf('32.47') + iv.mpf('0.5993')
assert mpf(aG2.b) < 2206 and mpf(bG2.b) < 1103
print(f"  a_SU(2) <= {mp.nstr(mpf(aG2.b),7)}, b_SU(2) <= {mp.nstr(mpf(bG2.b),7)}")
print("Haar part of the layer bound: kappa_G/3 * 1e-3 * (1+1e-5) <= 0.84 needs kappa_G <= 2519;")
for name, d, N, kap, rG in groups:
    assert Fr(kap) * Fr(1, 3000) * Fr(100001, 100000) <= Fr(84, 100)
print("  holds for every group of the table, and for SU(N) with N <= 5038.")
