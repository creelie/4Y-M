"""Certified constants for the block iterate.tex (mpmath interval arithmetic).

Every printed inequality is checked on intervals; a failed check raises.
"""
from mpmath import iv, mp
iv.prec = 120
mp.prec = 120

def le(a, b, msg):
    a = iv.mpf(a); b = iv.mpf(b)
    assert a.b <= b.a, msg
    print("  OK:", msg)

pi = iv.pi

def asin_encl(y):
    """Enclosure [lo,hi] of arcsin(y) for 0<y<1 by bisection on iv.sin (sin increasing on [0,pi/2])."""
    y = iv.mpf(y)
    lo = iv.mpf(0); hi = iv.mpf('1.5')
    assert iv.sin(hi).a > y.b
    for _ in range(300):
        mid = iv.mpf((lo.a + hi.a) / 2)
        sm = iv.sin(mid)
        if sm.a > y.b: hi = mid
        elif sm.b < y.a: lo = mid
        else: break
    assert iv.sin(lo).b < y.a and iv.sin(hi).a > y.b
    return lo, hi
log = iv.log
print("1. Global comparison constant of Theorem cutoff(ii)")
c2 = (pi / 3) ** 2
le(8 + 2688 * c2, 2956, "8 + 2688 (pi/3)^2 < 2956")
le(iv.mpf('1.6432e6') / 2956, 556, "1.6432e6 / 2956 < 556")
_al0, _ah0 = asin_encl(1 / (4 * iv.sqrt(2)))
kc0 = 2 * iv.sqrt(2) * iv.sin(iv.mpf([_al0.a, _ah0.b]) / 6)
Cstar = 1441 * 8 / kc0 ** 2
print("   C_* = 1441*8/kappa_c^2 =", Cstar)
le(555 * 2956, Cstar, "C_* > 555 * 2956, so 2956 < C_*/555")

print("2. Long-wavelength factor t_L = 3/4 + cos(4 pi/L)/4 at L = 32")
t32 = iv.mpf(3) / 4 + iv.cos(4 * pi / 32) / 4
print("   t_32 =", t32)
le(3.849, 4 * t32 ** 2, "4 t_32^2 > 3.849")
le(1 / t32 ** 2, 1.0393, "t_32^{-2} < 1.0393")
for L in (64, 128, 1024):
    tL = iv.mpf(3) / 4 + iv.cos(4 * pi / L) / 4
    print(f"   L={L}: 4 t_L^2 in", 4 * tL ** 2)

print("3. Norm conversion for SU(2): ||log g|| / delta(g,I) = x/sin x, x = angle/2, sin x = delta/(2 sqrt2)")
d = iv.mpf('1e-3')
_, h = asin_encl(d / (2 * iv.sqrt(2)))     # h >= x at delta = 1e-3; x/sin x increasing in x
r = h / iv.sin(h)
print("   ratio at delta=1e-3 (upper enclosure):", r)
le(r ** 2, 1 + d ** 2 / 20, "(x/sin x)^2 <= 1 + kappa_0^2/20 at kappa_0 = 1e-3")
le(r, 1 + d ** 2 / 40, "x/sin x <= 1 + kappa_0^2/40 at kappa_0 = 1e-3")
# hand bound used in the text: x <= 1.001 sin x, (x/sin x)^2 <= (1-x^2/6)^{-2} <= 1 + x^2/2.98, x^2 <= 1.002 delta^2/8
xx = iv.mpf('1.002') * d ** 2 / 8
le((1 - xx / 6) ** -2, 1 + xx / iv.mpf('2.98'), "(1-u)^{-2} <= 1+2.01u at u = x^2/6, x^2 <= 1.002e-6/8")
le(iv.mpf('1.002') / (8 * iv.mpf('2.98')), iv.mpf(1) / 20, "1.002/(8*2.98) < 1/20")
k0max = iv.mpf('9.99e-4')
le(6030 * k0max * k0max ** 2 / 40, k0max, "6030 kappa_0^3/40 <= kappa_0, so 1+6030 eta <= 1+6031 kappa_0")
le(k0max * (1 + k0max ** 2 / 40), iv.mpf('1e-3'), "eta <= kappa_0(1+kappa_0^2/40) < 1e-3 = rho_G/1000")

print("4. theta(kappa_0) = [(1+kappa_0^2/20)(1+6031 kappa_0)^2]^{-1}, kappa_1 = kappa_0/sqrt(28 theta)")
def theta(k0):
    k0 = iv.mpf(k0)
    return 1 / ((1 + k0 ** 2 / 20) * (1 + 6031 * k0) ** 2)
for k0 in ('1e-4', '1e-5', '1e-6', '1e-8'):
    th = theta(k0); k1 = iv.mpf(k0) / iv.sqrt(28 * th)
    print(f"   kappa_0={k0}: theta in [{float(th.a):.6f},{float(th.b):.6f}],"
          f"  kappa_1/kappa_0 in [{float((k1/iv.mpf(k0)).a):.5f},{float((k1/iv.mpf(k0)).b):.5f}]")
le(theta('1e-4'), iv.mpf('0.3892'), "theta(1e-4) < 0.3892"); le(iv.mpf('0.3890'), theta('1e-4'), "theta(1e-4) > 0.3890")
le(iv.mpf('0.9880'), theta('1e-6'), "theta(1e-6) > 0.9880")
le(iv.mpf('0.99987'), theta('1e-8'), "theta(1e-8) > 0.99987")
for k0 in ('1e-4', '1e-6', '1e-8'):
    le(iv.mpf('0.18'), iv.mpf(1) / iv.sqrt(28 * theta(k0)) , f"kappa_1/kappa_0 > 0.18 at kappa_0={k0}")
le(iv.mpf(1) / iv.sqrt(28 * theta('1e-6')), iv.mpf('0.1902'), "kappa_1/kappa_0 < 0.1902 at kappa_0=1e-6")
print("   loss of theta: log theta <= -2 log(1+6031 k) <= -2*6031 k/(1+6031 k) <= -1716 k for k <= 9.99e-4")
le(1716 * (1 + 6031 * k0max), 2 * 6031, "1716 (1+6031*9.99e-4) <= 12062")

print("5. Large-field thresholds at kappa_c (Theorems lfkp(ii) and sfcalm)")
_al, _ah = asin_encl(1 / (4 * iv.sqrt(2)))
kc = 2 * iv.sqrt(2) * iv.sin(iv.mpf([_al.a, _ah.b]) / 6)
le(iv.mpf('0.0837612'), kc, "kappa_c > 0.0837612"); le(kc, iv.mpf('0.0837614'), "kappa_c < 0.0837614")
print("   kappa_c =", kc)
def c1(lam):   # c_1(lambda,1)
    return 12 * log(iv.mpf('1.957')) + 45 * log(iv.mpf(14)) + iv.mpf('67.5') * (1 + log(11520 * 2 * lam))
def c2f(lam):  # c_2(lambda,1)
    return 24 * log(iv.mpf('1.957')) + 73 * log(iv.mpf(14)) + iv.mpf('109.5') * (1 + log(15979 * 2 * lam))
def a1(lam, k, logM):
    return lam * k ** 2 - 8 * logM - c1(lam)
def a2(lam, k, logM):
    return iv.mpf('0.52') * lam * k ** 2 - 8 * logM - c2f(lam)
logMhi = log(iv.mpf('2.3621e4'))
logMlo = iv.mpf(0)
def thresh(f, k, logM, lo, hi):
    # f increasing in lam on [lo,hi]; return [l,u] with f(l)<9 certified and f(u)>=9 certified
    lo = iv.mpf(lo); hi = iv.mpf(hi)
    assert f(lo, k, logM).b < 9 and f(hi, k, logM).a >= 9
    for _ in range(200):
        mid = iv.mpf((lo.a + hi.a) / 2)
        v = f(mid, k, logM)
        if v.a >= 9: hi = mid
        elif v.b < 9: lo = mid
        else: break
        if (hi.a - lo.a) / hi.a < 1e-9: break
    return lo, hi
for name, f, kap in (("beta^*  (a  >=9), kappa=0.08", a1, iv.mpf('0.08')),
                     ("beta^** (a_2>=9), kappa=0.08", a2, iv.mpf('0.08')),
                     ("beta^*  (a  >=9), kappa=kappa_c", a1, kc),
                     ("beta^** (a_2>=9), kappa=kappa_c", a2, kc)):
    for lm, tag in ((logMhi, "M=2.3621e4"), (logMlo, "M=1")):
        l, u = thresh(f, kap, lm, 1e3, 1e9)
        print(f"   {name}, {tag}: 4*lambda in [{float(4*l.a):.6e}, {float(4*u.a):.6e}]")
# monotonicity of a, a_2 in lambda beyond 67.5/kappa^2 and 109.5/(0.52 kappa^2)
le(iv.mpf('109.5') / (iv.mpf('0.52') * iv.mpf('0.0837') ** 2), 1e5, "a_2 increasing in lambda for lambda >= 1e5 at kappa >= 0.0837")

print("6. One-loop running for SU(2): lambda = beta_W/4 = 1/g^2, lambda_{k+1} = lambda_k - t")
b0 = iv.mpf(22) / (48 * pi ** 2)
t = 2 * b0 * log(iv.mpf(2))
print("   b_0 =", b0, "  t = 2 b_0 log 2 =", t, "  4t =", 4 * t, "  1/t =", 1 / t)
le(t, iv.mpf('0.064378'), "t < 0.064378"); le(iv.mpf('0.064377'), t, "t > 0.064377")
le(1 / t, iv.mpf('15.534'), "1/t < 15.534")
print("7. Scales not covered: (lambda_stop - 1)/t for lambda_stop = beta^**(kappa_c)/4 with M=1 (lower bound)")
l, u = thresh(a2, kc, logMlo, 1e3, 1e9)
print("   lambda^**(kappa_c), M=1 >", float(l.a), "  scales below it:", float(((l - 1) / t).a))
l2, u2 = thresh(a2, kc, logMhi, 1e3, 1e9)
print("   lambda^**(kappa_c), M=2.3621e4 <=", float(u2.a), "  scales below it:", float(((u2 - 1) / t).b))
print("8. Summation along the running")
for cE in (1, 10, 100):
    # sufficient lambda_K: c_E (lambda^-2 + 1/(t lambda)) <= 1/2
    lam = iv.mpf(1)
    while True:
        v = cE * (1 / lam ** 2 + 1 / (t * lam))
        if v.b <= iv.mpf(1) / 2: break
        lam = lam + 1
    print(f"   c_eps={cE}: sum eps_k <= 1/2 once lambda_K >= {int(lam.a)}  (beta_K = 4 lambda_K >= {int(4*lam.a)})")
print("9. Necessary conditions for the rates a >= 9 and a_2 >= 9 (log M >= 0, c_i increasing in lambda)")
c1_100 = c1(iv.mpf(100)); c2_100 = c2f(iv.mpf(100))
print("   c_1(100,1) =", c1_100, "  c_2(100,1) =", c2_100)
le(1190, 9 + c1_100, "9 + c_1(100,1) > 1190")
le(1930, (9 + c2_100) / iv.mpf('0.52'), "(9 + c_2(100,1))/0.52 > 1930")
le(iv.mpf(34.4), iv.sqrt(1190), "sqrt(1190) > 34.4")
C = 1716 * iv.sqrt(14) * iv.mpf('34.4')
print("   1716 sqrt(14) 34.4 =", C)
le(220000, C, "1716 sqrt14 * 34.4 > 2.2e5")
le(iv.mpf('1e11'), (220000 / iv.log(2)) ** 2, "(2.2e5/log 2)^2 > 1e11, so min{log 2, 2.2e5 lambda^-1/2} = log 2 for lambda <= 1e11")
print("10. Below the thresholds at kappa_c the rates are < 9 for every kappa < kappa_c and every lambda >= 100 (log M >= 0)")
le(67.5 / kc ** 2, 1e4, "a(.,kappa_c) increasing for lambda >= 1e4")
le(1e4 * kc ** 2 - c1(iv.mpf(100)), 9, "a(lambda,kappa_c) < 9 for 100 <= lambda <= 1e4")
le(iv.mpf('109.5') / (iv.mpf('0.52') * kc ** 2), 1e5, "a_2(.,kappa_c) increasing for lambda >= 1e5")
le(iv.mpf('0.52') * 1e5 * kc ** 2 - c2f(iv.mpf(100)), 9, "a_2(lambda,kappa_c) < 9 for 100 <= lambda <= 1e5")
lA, _ = thresh(a1, kc, logMlo, 1e3, 1e9)
lB, _ = thresh(a2, kc, logMlo, 1e3, 1e9)
le(245007, lA, "a(lambda,kappa_c) < 9 at lambda = 245007 with M = 1 (certified lower end of bisection)")
assert a1(iv.mpf(245007), kc, logMlo).b < 9
assert a2(iv.mpf(809322), kc, logMlo).b < 9
print("  OK: a(245007,kappa_c) < 9 and a_2(809322,kappa_c) < 9 with log M = 0")
_, uA = thresh(a1, kc, logMhi, 1e3, 1e9)
_, uB = thresh(a2, kc, logMhi, 1e3, 1e9)
assert a1(iv.mpf(256948), kc, logMhi).a >= 9 and a2(iv.mpf(832243), kc, logMhi).a >= 9
print("  OK: a(256948,kappa_c) >= 9 and a_2(832243,kappa_c) >= 9 with M = 2.3621e4")
le(iv.mpf('1.257e7'), (iv.mpf(809322) - 1) / t, "(809322-1)/t > 1.257e7 scales")
