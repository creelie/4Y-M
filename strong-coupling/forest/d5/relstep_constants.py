"""Certified constants for relstep.tex and NOTE.md (interval arithmetic,
mpmath.iv).  Every printed bound is an end point of an interval that
contains the exact value of the stated formula.

  c_*   = (1 + sqrt(21/8))^(-2)                       S_W >= c_* Q_1
  C_n   = 1 + 1/(4 sin^2(pi/(2n)))                    dist(w,G)^2 <= C_n Q_1(w)
  C_n^- = 1/(12 sin^2(pi/(2n+2)))                     no constant below C_n^- works
  eps'  = (eps_1 + delta_G/lambda) C_n / (2 c_*)      relative one-step error
  delta_G = 1/3 for SU(2)  (|Hess sum log j(0)| = |Tr ad_w^2|/12)

eps_1 is the constant of hypothesis (H); the text uses the example value
eps_1 = 1.63e-9 at lambda = 1e100.
"""
from mpmath import iv, mpf
iv.dps = 40
I = iv.mpf
pi = iv.pi

def up(x): return x.b
def lo(x): return x.a

print("== constants of the window lemma")
cstar = 1 / (1 + iv.sqrt(I(21) / 8)) ** 2
print(f"  c_* in [{float(lo(cstar)):.6f}, {float(up(cstar)):.6f}]   (claim c_* > 0.1456)")
assert lo(cstar) > mpf('0.1456')
# ingredients of c_*: entering averages and plaquettes
print("  7/256 * 8 =", 7 * 8 / 256, " (sum |v'|^2 <= (7/32) Q);  12*7/32 =", 12 * 7 / 32)

def Cn(n):  return 1 + 1 / (4 * iv.sin(pi / (2 * n)) ** 2)
def Cmin(n): return 1 / (12 * iv.sin(pi / (2 * n + 2)) ** 2)

print("== C_n (upper, proved) and C_n^- (lower, test field); asymptotics n^2/pi^2 and n^2/(3 pi^2)")
print("   n      C_n<=        1+n^2/4      C_n^->=      (n+1)^2/(3pi^2)<=   ratio C_n/C_n^-")
for n in [1, 2, 3, 4, 5, 6, 8, 16, 32, 64, 128, 1024, 10**4]:
    c, cm = Cn(n), Cmin(n)
    simple = 1 + I(n) ** 2 / 4
    asym = I(n + 1) ** 2 / (3 * pi ** 2)
    assert up(c) <= lo(simple) + mpf('1e-30')
    assert lo(cm) >= up(asym) - mpf('1e-30')
    print(f"  {n:6d} {float(up(c)):12.6g} {float(up(simple)):12.6g} {float(lo(cm)):12.6g} {float(up(asym)):14.6g} {float(up(c/cm)):12.4f}")

# epsilon of hypothesis (H): the example value used in the text
eps1 = I('1.63e-9')
lam = I('1e100'); deltaG = I(1) / 3
epsabs = eps1 + deltaG / lam

print("== relative constant eps' = (eps_1 + delta_G/lambda) C_n / (2 c_*)")
for n in [2, 4, 8, 16, 64, 256, 1024, 4096]:
    e = epsabs * Cn(n) / (2 * cstar)
    print(f"  n={n:5d}: eps' <= {float(up(e)):.4g}")
# largest n (power of two not needed) with eps' <= 0.1, via C_n <= 1 + n^2/4
def ok(n): return up(epsabs * (1 + I(n) ** 2 / 4) / (2 * cstar)) <= mpf('0.1')
lo_n, hi_n = 1, 2
while ok(hi_n): hi_n *= 2
while hi_n - lo_n > 1:
    mid = (lo_n + hi_n) // 2
    lo_n, hi_n = (mid, hi_n) if ok(mid) else (lo_n, mid)
nmax = lo_n          # ok(nmax) and not ok(nmax+1); C_n increases with n
assert ok(nmax) and not ok(nmax + 1)
print(f"  eps' <= 0.1 for every n <= {nmax} (checked through C_n <= 1 + n^2/4)")

print("== across scales: windows of side n_k = n_1 2^(1-k), k = 1..K, n_K >= 1")
print("   sum_k C_{n_k} <= K + n_1^2/3 (from C_n <= 1 + n^2/4 and sum 4^(1-k)/4 <= 1/3)")
for m in [1, 4, 8, 12]:
    n1 = 2 ** m; K = m + 1
    s = sum((Cn(2 ** (m - k)) for k in range(K)), I(0))
    bound = K + I(n1) ** 2 / 3
    assert up(s) <= lo(bound)
    tot = epsabs * s / (2 * cstar)
    print(f"  n_1=2^{m:2d}, K={K:2d}: sum C <= {float(up(s)):.6g} <= K + n_1^2/3 = {float(lo(bound)):.6g};"
          f"  sum eps'_k <= {float(up(tot)):.4g} (constant eps_1)")
K = 13
one = I(4) ** (K - 1) * eps1
allk = sum((I(4) ** (K - k) for k in range(1, K + 1)), I(0)) * eps1
print(f"   the same errors kept absolute (n_1=2^12, K={K}): the error of scale 1 at the last scale,"
      f" 4^(K-1) eps_1 <= {float(up(one)):.4g}; all K errors, sum_k 4^(K-k) eps_1 <= {float(up(allk)):.4g}")
