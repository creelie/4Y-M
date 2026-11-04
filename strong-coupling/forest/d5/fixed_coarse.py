"""Certified arithmetic for the fixed coarse field results (end of ssec:labelcov).
Interval arithmetic (mpmath.iv, 50 digits); a check passes only if it holds on the
whole enclosure.  Inputs, as in check_labelcov.py of ssec:labelcov:
alpha_G = 2.5e-4, sigma(alpha_G) in closed form, log10(b_* nu_s) <= 22.084 (s=1),
23.104 (s=5), log10 eps <= -14102.4 (s=1), -16652.7 (s=5).
Part U: Corollary cor:labelstretched (labels, R = ceil(200 sqrt d), t = floor(d/2),
1e8 <= d <= 1e100).
Part L: Theorem thm:fixedobstruction (ii) at eps = eps_s: L = log(e/eps),
H_0 = ceil(L/(2 alpha_G)), the exponent constant kappa_s = (2L+2alpha)/sqrt(H_0),
and the prefactor threshold e^{3L-6alpha}/(2e) at m = 3."""
from mpmath import iv, mp
iv.dps = 50
ln10 = iv.log(10)
def log10(x): return iv.log(x)/ln10
def hi(x): return x.b
def lo(x): return x.a
ok = []
def check(name, cond):
    ok.append((name, bool(cond))); print(("PASS " if cond else "FAIL ") + name)
def show(name, x): print(f"{name} in [{mp.nstr(lo(x),12)}, {mp.nstr(hi(x),12)}]")
def sig(kap):
    x = iv.exp(-kap)
    return 1 + 64*x*(1+4*x+x**2)/(1-x)**4 + 16*x/(1-x)**2

aG = iv.mpf('2.5e-4')
s1 = sig(aG)
check("sigma(alpha_G) <= 9.831e16", hi(s1) <= 9.831e16)
data = {1: dict(bn=iv.mpf('22.084'), L10e=iv.mpf('-14102.4')),
        5: dict(bn=iv.mpf('23.104'), L10e=iv.mpf('-16652.7'))}
d0 = iv.mpf('1e8'); sd0 = iv.sqrt(d0)
print("\n=== Part U: d-independent conditions at d = 1e8 ===")
check("t >= 2R: d/2 - 1 >= 2(200 sqrt d + 1) at d = 1e8 (and d/2 - 400 sqrt d increases)",
      lo(d0/2 - 1 - 2*(200*sd0 + 1)) >= 0)
check("t/R >= sqrt d/401: d/2 >= sqrt d + 401 at d = 1e8", lo(d0/2 - sd0 - 401) >= 0)
check("2R+1 <= 401 sqrt d: sqrt d >= 3", lo(sd0) >= 3)
check("sqrt d/401 - 3 > 0 at d = 1e8", lo(sd0/401 - 3) > 0)
check("main ratio decreases: alpha_G > 1/(40 sqrt d) at d = 1e8", lo(aG - 1/(40*sd0)) > 0)
import math
check("d <= 1e100: ceil(200 sqrt d) <= 200*10^50 + 1", math.isqrt(10**100) == 10**50)
Rmax = 200*iv.mpf(10)**50 + 1
outU = {}
for s, D in data.items():
    print(f"\n=== Part U, s = {s} ===")
    bn = iv.mpf(10)**D['bn']; L10e = D['L10e']; lneps = L10e*ln10
    # far islands: 2 b nu sigma eps^{1/2} e^{alpha(4s+17)} e^{-200 alpha sqrt d}, 200 alpha = 1/20
    from fractions import Fraction
    check("200 alpha_G = 1/20 (exact rationals)", 200*Fraction(25, 100000) == Fraction(1, 20))
    far = log10(2*bn*s1) + L10e/2 + aG*(4*s+17)/ln10
    show("log10 far coefficient", far)
    # main: 11^4 b nu e^{2a(2s+14)} e^{-alpha d} / e^{-sqrt d/20} at d = 1e8
    main = log10(14641*bn) + (2*aG*(2*s+14) - aG*d0 + sd0/20)/ln10
    show("log10 main/e^{-sqrt d/20} at d = 1e8", main)
    lam = -(lneps/4 + iv.log(iv.sqrt(2)) + 4*iv.log(2*Rmax+1))
    show("-ln g_R at R = 200*10^50+1", lam)
    check(f"s={s}: g_R <= 1/2 for R <= 200*10^50+1", lo(lam) >= hi(iv.log(2)))
    check(f"s={s}: lambda/401 > 1/20 + 4/sqrt(1e8) (chain ratio decreases)", lo(lam/401) > hi(iv.mpf(1)/20 + 4/sd0))
    lamlo = iv.mpf(lo(lam))
    chain = log10(4*iv.mpf(401)**4*d0**2) + (3*lamlo - (lamlo/401 - iv.mpf(1)/20)*sd0)/ln10
    show("log10 chain/e^{-sqrt d/20} at d = 1e8", chain)
    tot = iv.mpf(10)**hi(far) + iv.mpf(10)**hi(main) + iv.mpf(10)**hi(chain)
    show("log10 C_U", log10(tot))
    outU[s] = hi(log10(tot))
    show("log10 eps C_U (absolute coefficient)", L10e + log10(tot))
    outU[(s,'abs')] = hi(L10e + log10(tot))

print("\n=== Part L: lower bounds at eps = eps_s ===")
outL = {}
for s, D in data.items():
    L10e = D['L10e']
    L = 1 - L10e*ln10                     # L = log(e/eps)
    show(f"s={s}: L = log(e/eps)", L)
    q = L/(2*aG)
    show("L/(2 alpha_G)", q)
    H0 = math.ceil(float(hi(q)))
    check(f"s={s}: ceil(L/(2 alpha_G)) determined", math.ceil(float(lo(q))) == H0)
    print("H_0 =", H0)
    check(f"s={s}: H_0 >= 2s+9", H0 >= 2*s+9)
    kap = (2*L + 2*aG)/iv.sqrt(iv.mpf(H0))
    show("kappa_s = (2L+2alpha)/sqrt(H_0)", kap)
    thr = (3*L - 6*aG - 1 - iv.log(2))/ln10
    show("log10 e^{3L-6alpha}/(2e)", thr)
    d9 = 9*H0
    print("d at m = 3: 9 H_0 =", d9)
    outL[s] = dict(L=L, H0=H0, kap=hi(kap), thr=lo(thr), d9=d9)
    show("log10 of lower bound prefactor e^{-L}/2", (-L - iv.log(2))/ln10)

print("\n=== stated numbers ===")
check("s=1: C_U <= 10^-7011.7", outU[1] <= -7011.7)
check("s=5: C_U <= 10^-8285.8", outU[5] <= -8285.8)
check("s=1: eps C_U <= 10^-21114.2", outU[(1,'abs')] <= -21114.2)
check("s=5: eps C_U <= 10^-24938.6", outU[(5,'abs')] <= -24938.6)
check("s=1: kappa_1 <= 8.06", outL[1]['kap'] <= 8.06)
check("s=5: kappa_5 <= 8.76", outL[5]['kap'] <= 8.76)
check("s=1: threshold log10 >= 42307", outL[1]['thr'] >= 42307)
check("s=5: threshold log10 >= 49958", outL[5]['thr'] >= 49958)
check("s=1: 9 H_0 <= 5.85e8", outL[1]['d9'] <= 5.85e8)
check("s=5: 9 H_0 <= 6.91e8", outL[5]['d9'] <= 6.91e8)
check("s=1: H_0 >= 1e8/... i.e. 4 H_0 >= 1e8 (m >= 2 lies in [1e8, 1e100])", 4*outL[1]['H0'] >= 1e8)
check("s=5: 4 H_0 >= 1e8", 4*outL[5]['H0'] >= 1e8)
print("\nALL PASS" if all(v for _, v in ok) else "\nSOME FAIL")
