"""Constants of the block bigwindow.tex (clusters of islands).

Interval arithmetic (mpmath.iv, 60 digits).  A claim "X <= c" is accepted
only when the upper end of the interval enclosing X is <= c.  Large and small
numbers are handled through their logarithms.  The data are those of
Theorem thm:islandwindow: rho = 2.17e-16, r = 43 rho, H = 6^4, M = 11^4,
alpha_G = 2.5e-4, beta_W = beta*.
"""
from mpmath import iv, mp

iv.dps = 60
mp.dps = 60

H, Mbox = 1296, 14641
N5 = 11 ** 4 - 1                      # neighbours at distance <= 5
alphaG = iv.mpf(1) / 4000
rho = iv.mpf(217) / iv.mpf(10) ** 18
r = 43 * rho
C0 = iv.mpf(1) / 2 + iv.mpf(2) / 3 * iv.log(24 * iv.pi ** 3)
LN10 = iv.log(10)
def log10(x): return iv.log(x) / LN10
def up(x): return mp.nstr(mp.mpf(x._mpi_[1]), 12)
def lo(x): return mp.nstr(mp.mpf(x._mpi_[0]), 12)

checks = []
def check(name, ok):
    checks.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name)

def sigma(k):
    x = iv.exp(-k)
    return 1 + 64 * x * (1 + 4 * x + x * x) / (1 - x) ** 4 + 16 * x / (1 - x) ** 2

# M_1 = 10^(10^1700): ln M_1 = 10^1700 ln 10
lnM1 = iv.mpf(10) ** 1700 * LN10
m0 = 10 ** 5

for s in (1, 5):
    print("\n================ s =", s, "================")
    Ns = (4 * s + 17) ** 4 - 1
    logq = -(H * (iv.log(2) + 2 * iv.log(Ns)) + 1 + iv.log(Mbox * H))
    def gpt(b):
        b = mp.mpf(b); rh = mp.mpf(217) / mp.mpf(10) ** 18
        c0 = mp.mpf(1) / 2 + mp.mpf(2) / 3 * mp.log(24 * mp.pi ** 3)
        lq = -(H * (mp.log(2) + 2 * mp.log(Ns)) + 1 + mp.log(Mbox * H))
        return 6 * (mp.log(b) + c0) - b * rh ** 2 / 16 - lq
    a, b = mp.mpf(96) / (mp.mpf(217) / mp.mpf(10) ** 18) ** 2, mp.mpf(10) ** 40
    for _ in range(400):
        m = (a + b) / 2
        if gpt(m) > 0: a = m
        else: b = m
    blo = iv.mpf(a * (1 - mp.mpf(10) ** -20)); bhi = iv.mpf(b * (1 + mp.mpf(10) ** -20))
    g = lambda be: 6 * (iv.log(be) + C0) - be * rho ** 2 / 16 - logq
    check("s=%d: root bracket beta* in [%s, %s]" % (s, lo(blo), up(bhi)),
          g(blo).a > 0 and g(bhi).b < 0)
    # monotonicity in beta: beta* rho^2/16 >= 32481 (s=1) / 38354 (s=5)
    slope = blo * rho ** 2 / 16
    print("beta* rho^2/16 >=", lo(slope))
    check("s=%d: beta* rho^2/16 >= %d" % (s, {1: 32481, 5: 38354}[s]),
          slope.a >= {1: 32481, 5: 38354}[s])
    l10Mq = log10(iv.mpf(Mbox)) + logq / LN10
    l10A = log10(iv.e * Mbox * H) + logq / LN10
    print("log10 Mq* <=", up(l10Mq), "  log10 A_* <=", up(l10A))
    for th, lab, R in ((iv.mpf(0), "0", 10 ** 7), (iv.mpf(1) / 2, "1/2", 2 * 10 ** 7)):
        kap = 2 * (1 - th) * alphaG
        print("--- vartheta =", lab, " kappa =", lo(kap), " R =", R)
        Lam = iv.mpf('99.0001') * bhi + iv.mpf('4.65e6')
        c1 = 1920 * Lam ** 2 * r ** 2 / bhi * 45 * (2 * s + 7) ** 4 * iv.exp(kap * (2 * s + 9))
        l10c1 = log10(c1)
        sk2 = sigma(kap / 2)
        l10tau = l10c1 - kap * R / 2 / LN10 + log10(sk2)
        print("log10 c_1(beta*) <=", up(l10c1), "  sigma(kappa/2) <=", up(sk2))
        print("log10 tau <=", up(l10tau))
        target = {(1, "0"): -1047, (1, "1/2"): -1045, (5, "0"): -1045, (5, "1/2"): -1044}[(s, lab)]
        check("s=%d theta=%s R=%d: tau <= 1e%d" % (s, lab, R, target), l10tau.b <= target)
        # reach l_m = R + (6/kappa) ln m
        def ell(lnm): return R + 6 / kap * lnm
        lM1 = ell(lnM1)
        l10lM1 = log10(lM1)
        print("log10 l_{M_1} in [", lo(l10lM1), ",", up(l10lM1), "]")
        check("s=%d theta=%s: l_{M_1} <= 10^1705" % (s, lab), l10lM1.b <= 1705)
        # nu_l at l = 10^1705, an upper bound for every l <= l_{M_1}
        l10N = 4 * log10(2 * iv.mpf(10) ** 1705 + 1)
        l10NA = l10N + l10A / 2
        print("log10 (N_l A^{1/2}) at l=1e1705 <=", up(l10NA))
        nu = 4 * (N5 + iv.mpf(10) ** l10NA) * iv.mpf(10) ** (l10A / (2 * H))
        l10nu = log10(nu)
        print("nu <=", up(nu), "  log10 nu <=", up(l10nu))
        nub = {1: 0.2130, 5: 0.0222}[s]
        check("s=%d theta=%s: nu_l <= %s for l <= l_{M_1}" % (s, lab, nub), nu.b <= nub)
        # (ii): |Lambda| A^{(H+1)/2H} nu^{M_1}/(1-nu) <= |Lambda| 10^{-c M_1}
        pref = (H + 1) / iv.mpf(2 * H) * l10A - log10(1 - nu)
        print("log10 prefactor A^{(H+1)/2H}/(1-nu) <=", up(pref))
        cM = {1: -0.67, 5: -1.65}[s]
        check("s=%d theta=%s: global bound <= |Lambda| 10^{%s M_1}" % (s, lab, cM),
              pref.b <= 0 and l10nu.b <= cM)
        # (iii): m_0 = 1e5
        lm0 = ell(iv.log(m0))
        near = 4 * log10(2 * (m0 - 1) * lm0 + 1) + 4 * log10(2 * lm0 + 1) + 2 * l10Mq
        tail = pref + m0 * l10nu
        print("l_{m0} <=", up(lm0), " log10 near <=", up(near), " log10 tail <=", up(tail))
        tot = log10(iv.mpf(10) ** near + iv.mpf(10) ** tail)
        print("log10 (near + tail) <=", up(tot))
        tgt = {(1, "0"): -28126, (1, "1/2"): -28124, (5, "0"): -33227, (5, "1/2"): -33224}[(s, lab)]
        check("s=%d theta=%s: P(y bad, U(y) not an island) <= 1e%d + global" % (s, lab, tgt),
              tot.b <= tgt)
        check("s=%d theta=%s: same with the factor 1+1e-6 (beta_W > beta*)" % (s, lab),
              (tot + log10(1 + iv.mpf(10) ** -6)).b <= tgt)
    # beta_W > beta*: (1+1e-7+4e-4 log t)^8 t^12 e^{-64962(t-1)} <= 1+1e-6 for t >= 1
    #   (value at t=1 is (1+1e-7)^8; log-derivative 3.2e-3/t + 12/t - 64962 < 0)
    v1 = (1 + iv.mpf(10) ** -7) ** 8
    check("s=%d: (1+1e-7)^8 <= 1+1e-6" % s, v1.b <= 1 + 1e-6)
    # slack of the stated near terms against a factor 1+1e-6 is checked below via the totals
    # where the counting stops (theta = 0 reach does not matter here)
    # nu_l >= 1 iff N_l >= A^{-1/2}(A^{-1/2H}/4 - N5)
    l10crit = -l10A / 2 + log10(iv.mpf(10) ** (-l10A / (2 * H)) / 4 - N5)
    l10lc = (l10crit) / 4 - log10(iv.mpf(2))
    print("nu_l = 1 at log10 l about", lo(l10lc), up(l10lc))
    l10l1 = (-l10A) / 4 - log10(iv.mpf(2))
    print("N_l A_* = 1 at log10 l about", lo(l10l1), up(l10l1))
    def nu_at(l10l):
        l10Nl = log10((2 * iv.mpf(10) ** l10l + 1) ** 4 - 1)
        return 4 * (N5 + iv.mpf(10) ** (l10Nl + l10A / 2)) * iv.mpf(10) ** (l10A / (2 * H))
    if s == 1:
        n1, n2 = nu_at(iv.mpf('1763.2')), nu_at(iv.mpf('1763.3'))
        print("nu at log10 l = 1763.2 <=", up(n1), "; at 1763.3 >=", lo(n2))
        check("s=1: nu_l < 1 at log10 l = 1763.2 (A = eMHq*)", n1.b < 1)
        check("s=1: nu_l > 1 at log10 l = 1763.3 (A = eMHq*)", n2.a > 1)
        # N_l A_* at log10 l = 3524.5, with A_* the enclosure's upper end replaced by its value
        l10NA = log10((2 * iv.mpf(10) ** iv.mpf('3524.5') + 1) ** 4 - 1) + l10A
        print("log10 N_l A_* at log10 l = 3524.5 in", lo(l10NA), up(l10NA))
        check("s=1: N_l A_* >= 1 at log10 l = 3524.5", l10NA.a >= 0)

print("\nall checks passed:", all(ok for _, ok in checks), "(%d checks)" % len(checks))
