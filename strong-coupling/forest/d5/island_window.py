"""Constants of Theorem thm:islandwindow (the island sums in a window).

Every inequality used in the statement is checked in interval arithmetic
(mpmath.iv): a claim "X <= c" is accepted only when the upper end of the
interval enclosing X is <= c.  Large and small numbers are handled through
their natural or decimal logarithms.
"""
from mpmath import iv, mp

iv.dps = 60
mp.dps = 60

H, M = 1296, 14641                    # overlap 6^4 and witness box 11^4
alphaG = iv.mpf(1) / 4000             # 2.5e-4
rho = iv.mpf(217) / iv.mpf(10) ** 18  # 2.17e-16
r = 43 * rho
C0 = iv.mpf(1) / 2 + iv.mpf(2) / 3 * iv.log(24 * iv.pi ** 3)
LN10 = iv.log(10)

def log10(x):
    return iv.log(x) / LN10

checks = []
stated = []
claims = {}
log10G0_s = {}
def check(name, ok):
    checks.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name)

def up(x):   # upper end of an interval, as a float string
    return mp.nstr(mp.mpf(x._mpi_[1]), 14)
def lo(x):
    return mp.nstr(mp.mpf(x._mpi_[0]), 14)

def sigma(kappa):
    x = iv.exp(-kappa)
    return 1 + 64 * x * (1 + 4 * x + x * x) / (1 - x) ** 4 + 16 * x / (1 - x) ** 2

print("rho in", lo(rho), up(rho), "  r = 43 rho, eta = 92 rho")
print("eta + 62 r = 2758 rho <= 6.0e-13 :", up(2758 * rho))
check("hypothesis (a): 2758 rho <= 6.0e-13", (2758 * rho).b <= 6e-13)
check("Lemma lambdatwo: eta + 30 r = 1382 rho <= 6.0e-13", (1382 * rho).b <= 6e-13)

for s in (1, 5):
    print("\n================ s =", s, "================")
    N = (4 * s + 17) ** 4 - 1
    logqstar = -(H * (iv.log(2) + 2 * iv.log(N)) + 1 + iv.log(M * H))
    print("N_s =", N, "  |log q*| in [", lo(-logqstar), ",", up(-logqstar), "]")

    def g_pt(b):  # point version for bisection
        b = mp.mpf(b)
        rh = mp.mpf(217) / mp.mpf(10) ** 18
        c0 = mp.mpf(1) / 2 + mp.mpf(2) / 3 * mp.log(24 * mp.pi ** 3)
        lq = -(H * (mp.log(2) + 2 * mp.log(N)) + 1 + mp.log(M * H))
        return 6 * (mp.log(b) + c0) - b * rh ** 2 / 16 - lq
    a, b = mp.mpf(96) / (mp.mpf(217) / mp.mpf(10) ** 18) ** 2, mp.mpf(10) ** 40
    for _ in range(400):
        m = (a + b) / 2
        if g_pt(m) > 0: a = m
        else: b = m
    blo, bhi = iv.mpf(a * (1 - mp.mpf(10) ** -20)), iv.mpf(b * (1 + mp.mpf(10) ** -20))
    def g(beta):
        return 6 * (iv.log(beta) + C0) - beta * rho ** 2 / 16 - logqstar
    check("s=%d: g(beta_lo) > 0 > g(beta_hi), beta_lo > 96/rho^2" % s,
          g(blo).a > 0 and g(bhi).b < 0 and (blo - 96 / rho ** 2).a > 0)
    print("beta* in [", lo(blo), ",", up(bhi), "]")
    check("s=%d: beta* > 112/rho^2 (monotonicity of Lambda^2 e^f / beta)" % s,
          (blo - 112 / rho ** 2).a > 0)
    # thm:islandthreshold(i): A <= (2 N^2)^{-H} at beta*, i.e. x = N^2 A^{1/H} <= 1/2
    log10Mq = log10(iv.mpf(M)) + logqstar / LN10
    log10A = log10(iv.e * M * H) + logqstar / LN10
    print("log10(M q*) <=", up(log10Mq), "   log10 A* <=", up(log10A))
    logx = 2 * iv.log(N) + (1 + iv.log(M * H) + logqstar) / H
    print("x = N^2 A*^{1/H}: log x in", lo(logx), up(logx), " (log 1/2 =", mp.nstr(-mp.log(2), 12), ")")
    check("s=%d: x <= 1/2 at beta*" % s, (logx + iv.log(2)).b <= 1e-40)
    # G_0
    log10G0 = log10(iv.mpf(H)) - H * log10(iv.mpf(2)) - log10(iv.mpf(N) ** 2 - 1) \
        - 2 * log10(1 - iv.mpf(2) ** (-H))
    print("log10 G_0 <=", up(log10G0))
    log10G0_s[s] = log10G0
    for theta, label in ((iv.mpf(0), "0"), (iv.mpf(1) / 2, "1/2")):
        kappa = 2 * (1 - theta) * alphaG
        sk, sk2 = sigma(kappa), sigma(kappa / 2)
        print("--- vartheta =", label, " kappa =", lo(kappa))
        print("sigma(kappa) in [", lo(sk), ",", up(sk), "]  sigma(kappa/2) <=", up(sk2))
        # annealed island mass: paper's bound at Gamma=1 vs G_0
        paper = 2 * sk / iv.mpf(N) ** 2
        print("prop:islandmass(ii) at Gamma=1 per |J|:", lo(paper), up(paper))
        print("new mean bound G_0 sigma(kappa): log10 <=", up(log10G0 + log10(sk)))
        # c_2 at beta_hi (c_2 increases with beta)
        Lam = iv.mpf('99.0001') * bhi + iv.mpf('4.65e6')
        c2 = 1920 * Lam ** 2 * r ** 2 / bhi * 45 * (2 * s + 12) ** 2 * (2 * s + 7) ** 2 \
            * iv.exp(kappa * (2 * s + 14))
        print("log10 c_2(beta*) <=", up(log10(c2)))
        check("s=%d theta=%s: log10 c_2 <= %s" % (s, label, {1: "22.13", 5: "23.14"}[s]),
              log10(c2).b <= {1: 22.13, 5: 23.14}[s])
        check("s=%d theta=%s: log10 G_0 sigma(kappa) <= %s" % (s, label, {1: "-380.6", 5: "-382.5"}[s]),
              (log10G0 + log10(sk)).b <= {1: -380.6, 5: -382.5}[s])
        for R in (2 * s + 13, int(mp.ceil(2 * mp.mpf(10) ** 8 / (1 - mp.mpf(theta.a))))):
            near = log10(iv.mpf(11) ** 4) + 4 * log10(iv.mpf(2 * R + 1)) + 2 * log10Mq
            if R == 2 * s + 13:
                ltau = log10(sk)
            else:
                ltau = -kappa * R / 2 / LN10 + log10(sk2)
            farpref = log10Mq + log10(c2) + log10G0 + ltau   # times 1/sigma
            print("R = %d: log10 near/|W| <= %s ; log10 far/|W| <= %s - log10(sigma)"
                  % (R, up(near), up(farpref)))
            for lsig in (-100, -1000):
                tot = up(log10(10 ** near + 10 ** (farpref - lsig)))
                print("    sigma = 1e%d: log10 eps/|W| <= %s" % (lsig, tot))
                claims[(s, label, R, lsig)] = (near, farpref - lsig)
                if R != 2 * s + 13:
                    check("s=%d theta=%s R=%d sigma=1e%d: far <= near"
                          % (s, label, R, lsig), (farpref - lsig - near).b <= 0)
        # reference: decimal statements used in the theorem
    if s == 1:
        check("s=1: log10(M q*) <= -14102.4", log10Mq.b <= -14102.4)
        check("s=1: log10 A* <= -14098.8", log10A.b <= -14098.8)
        check("s=1: log10 G_0 <= -397.6", log10G0.b <= -397.6)
    if s == 5:
        check("s=5: log10(M q*) <= -16652.79", log10Mq.b <= -16652.79)
        check("s=5: log10 G_0 <= -399.56", log10G0.b <= -399.56)
    stated.append((s, blo, bhi))

print("\n================ statements of the theorem ================")
for (s, blo_, bhi_) in stated:
    print("s=%d: beta* in [%s, %s]" % (s, lo(blo_), up(bhi_)))
check("s=1: beta* in [1.1220e37, 1.1221e37]", stated[0][1].a >= 1.1220e37 and stated[0][2].b <= 1.1221e37)
check("s=5: beta* in [1.3216e37, 1.3217e37]", stated[1][1].a >= 1.3216e37 and stated[1][2].b <= 1.3217e37)
def tot(k):
    n, f = claims[k]
    return log10(iv.mpf(10) ** n + iv.mpf(10) ** f)
check("s=1, theta=0, R=15: near <= 1e-28194 |W|", claims[(1, "0", 15, -100)][0].b <= -28194)
check("s=1, theta=0, R=15: far <= 1e-14462 |W| / sigma", (claims[(1, "0", 15, -100)][1] - 100).b <= -14462)
check("s=1, theta=0, R=2e8, sigma>=1e-1000: eps <= 1e-28166 |W|", tot((1, "0", 200000000, -1000)).b <= -28166)
check("s=1, theta=1/2, R=4e8, sigma>=1e-1000: eps <= 1e-28165 |W|", tot((1, "1/2", 400000000, -1000)).b <= -28165)
check("s=5, theta=0, R=2e8, sigma>=1e-1000: eps <= 1e-33267 |W|", tot((5, "0", 200000000, -1000)).b <= -33267)
check("s=1: G_0 sigma(2 alpha_G) <= 1e-381.8", (log10G0_s[1] + log10(sigma(2 * alphaG))).b <= -381.8)
check("s=5: G_0 sigma(2 alpha_G) <= 1e-383.7", (log10G0_s[5] + log10(sigma(2 * alphaG))).b <= -383.7)
check("|W| <= 1e28000 gives eps <= 1e-166 (s=1, theta=0)", (tot((1, "0", 200000000, -1000)) + 28000).b <= -166)
# the paper's value 1.13e37 of cor:twoislands corresponds to rho = 2.16e-16
rho2 = iv.mpf(216) / iv.mpf(10) ** 18
for s, val in ((1, 1.13e37), (5, 1.33e37)):
    N = (4 * s + 17) ** 4 - 1
    lq = -(H * (iv.log(2) + 2 * iv.log(N)) + 1 + iv.log(M * H))
    def g2(beta): return 6 * (iv.log(beta) + C0) - beta * rho2 ** 2 / 16 - lq
    print("rho=2.16e-16, s=%d: g(%.4g) = %s, g(%.4g) = %s" % (s, val * 0.995, lo(g2(iv.mpf(val * 0.995))), val * 1.005, up(g2(iv.mpf(val * 1.005)))))
print("\nsigma(5e-4) (paper: 6.144e15):", lo(sigma(2 * alphaG)), up(sigma(2 * alphaG)))
print("\nall checks passed:", all(ok for _, ok in checks), "(%d checks)" % len(checks))
