"""Exact rational re-check of every inequality written by hand in the
patches L01-L06 (lfkp, the beta** table, lfinherited(iii), the cluster
numbers). Only Fractions; logarithms by the artanh series with a bounded
tail, as in the text. Run: python3 lfhand_check.py > lfhand_check.out"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
from rat import I, log_gt1

fails = 0; count = 0
def ok(name, cond):
    global fails, count
    count += 1
    if not cond: fails += 1
    print(("PASS " if cond else "FAIL ") + name)

def enc(y, n=80): return I(*log_gt1(F(y), n))

# --- the log table: series values, then the 8-decimal bounds of the text
s2 = enc(2); s32 = enc(F(3, 2)); s54 = enc(F(5, 4)); s5049 = enc(F(50, 49))
s1110 = enc(F(11, 10)); s3736 = enc(F(37, 36))
for nm, x in [("log2", s2), ("log3/2", s32), ("log5/4", s54), ("log50/49", s5049),
              ("log11/10", s1110), ("log37/36", s3736)]:
    ok("series enclosure width < 1e-10: " + nm, x.hi - x.lo < F(1, 10**10))
eL2 = s2; eL3 = eL2 + s32; eL5 = 2 * eL2 + s54
eL7 = (eL2 + 2 * eL5 - s5049) * F(1, 2)     # 50/49 = 2*5^2/7^2
eL11 = eL2 + eL5 + s1110; eL37 = 2 * eL2 + 2 * eL3 + s3736
table = {2: ("0.69314718", "0.69314719"), 3: ("1.09861228", "1.09861229"),
         5: ("1.60943791", "1.60943792"), 7: ("1.94591014", "1.94591015"),
         11: ("2.39789527", "2.39789528"), 37: ("3.61091791", "3.61091792")}
exact = {2: eL2, 3: eL3, 5: eL5, 7: eL7, 11: eL11, 37: eL37}
for p, (lo, hi) in table.items():
    ok("%s < L_%d < %s" % (lo, p, hi), F(lo) < exact[p].lo and exact[p].hi < F(hi))
L = {p: I(F(lo), F(hi)) for p, (lo, hi) in table.items()}
L2, L3, L5, L7, L11, L37 = (L[p] for p in (2, 3, 5, 7, 11, 37))
L10 = L2 + L5
def logint(n):
    """log of a positive integer whose prime factors are in the table."""
    r = I(0); m = n
    for p in (2, 3, 5, 7, 11, 37):
        while m % p == 0: r = r + L[p]; m //= p
    assert m == 1, n
    return r

# --- identities used for the logarithms
ok("2.4e4 = 2^6 3 5^3", 2**6 * 3 * 5**3 == 24000)
ok("1.96 = 7^2/5^2", F(196, 100) == F(49, 25))
ok("23040 = 2^9 3^2 5", 2**9 * 9 * 5 == 23040)
ok("32000 = 2^8 5^3 > 2*15979", 2**8 * 125 == 32000 and 2 * 15979 < 32000)
ok("7290 = 2 3^6 5", 2 * 3**6 * 5 == 7290)

# --- lfkp (i)
ok("8.9 > log 7290", F("8.9") > logint(7290).hi)
ok("81/(7290-6400) < 0.092 < 1/10", F(81, 7290 - 6400) < F("0.092") < F(1, 10))

# --- K0 and K2
K0 = 67.5 + 700.5 * L2 + 143 * L3 + 67.5 * L5 + 69 * L7
logM = logint(24000); log196 = 2 * L7 - 2 * L5; log14 = L2 + L7
def vec(**kw): return kw
def addv(*terms):
    r = {}
    for coef, v in terms:
        for k, x in v.items(): r[k] = r.get(k, 0) + coef * x
    return {k: x for k, x in r.items() if x != 0}
vM = {"L2": 6, "L3": 1, "L5": 3}; v196 = {"L7": 2, "L5": -2}; v14 = {"L2": 1, "L7": 1}
v23040 = {"L2": 9, "L3": 2, "L5": 1}; v32000 = {"L2": 8, "L5": 3}
k0v = addv((8, vM), (12, v196), (45, v14), (F("67.5"), v23040))
ok("K0 coefficients 700.5 L2 + 143 L3 + 67.5 L5 + 69 L7 (+67.5)", k0v == {"L2": F("700.5"), "L3": 143, "L5": F("67.5"), "L7": 69})
ok("K0 < 953.06", K0.hi < F("953.06"))
K2 = 109.5 + 997 * L2 + 8 * L3 + 304.5 * L5 + 121 * L7
k2v = addv((8, vM), (24, v196), (73, v14), (F("109.5"), v32000))
ok("K2 coefficients 997 L2 + 8 L3 + 304.5 L5 + 121 L7 (+109.5)", k2v == {"L2": 997, "L3": 8, "L5": F("304.5"), "L7": 121})
ok("K2 < 1534.89", K2.hi < F("1534.89"))

# --- lfkp table: beta* = 4 lambda*, t* = lambda* kappa^2
rows = [(F("0.08"), F("1.13125e6"), 1810, 1800, F(1, 180), 4 * L5 - 2 * L2),
        (F("1e-2"), F("8.404e7"), 2101, 2100, F(1, 2100), 4 * L10),
        (F("1e-4"), F("1.0964e12"), 2741, 2744, F(0), 8 * L10),
        (F("1e-8"), F("1.604e20"), 4010, 4000, F(1, 400), 16 * L10)]
for kap, beta, t, base, add, lk in rows:
    ok("lfkp row kappa=%s: beta/4*kappa^2 = %d" % (kap, t), beta / 4 * kap**2 == t)
    ok("  log kappa^-2 identity", F(1) / kap**2 == {F("0.08"): F(625, 4)}.get(kap, F(1) / kap**2))
    ok("  log t <= log %d + %s (log(1+x)<=x)" % (base, add), (t <= base and add == 0) or F(t - base, base) <= add)
    logt_hi = logint(base).hi + add
    margin = t - F("67.5") * logt_hi - F("962.06") - F("67.5") * lk.hi
    ok("  t - 67.5 log t - 962.06 - 67.5 log kappa^-2 > 0.6  (%.4f)" % float(margin), margin > F("0.6"))
    ok("  t > 67.5 and lambda >= 100", t > F("67.5") and beta / 4 >= 100)
ok("factorizations 1800 2100 2744 4000", (2**3 * 9 * 25, 4 * 3 * 25 * 7, 8 * 343, 32 * 125) == (1800, 2100, 2744, 4000))

# --- beta** table
rows2 = [(F("0.08"), F("3.6625e6"), 5860, 5832, F(7, 1458), 4 * L5 - 2 * L2),
         (F("1e-2"), F("2.7064e8"), 6766, 6750, F(8, 3375), 4 * L10),
         (F("1e-4"), F("3.504e12"), 8760, 8750, F(1, 875), 8 * L10),
         (F("1e-8"), F("5.0872e20"), 12718, 12600, F(59, 6300), 16 * L10)]
for kap, beta, t, base, add, lk in rows2:
    ok("beta** row kappa=%s: beta/4*kappa^2 = %d" % (kap, t), beta / 4 * kap**2 == t)
    ok("  (t-base)/base == %s" % add, F(t - base, base) == add)
    logt_hi = logint(base).hi + add
    margin = F("0.52") * t - F("109.5") * logt_hi - F("1543.89") - F("109.5") * lk.hi
    ok("  0.52t - 109.5 log t - 1543.89 - 109.5 log kappa^-2 > 0.14  (%.4f)" % float(margin), margin > F("0.14"))
    ok("  t > 211", t > 211)
ok("factorizations 5832 6750 8750 12600", (8 * 729, 2 * 27 * 125, 2 * 625 * 7, 8 * 9 * 25 * 7) == (5832, 6750, 8750, 12600))

# --- lfinherited (iii)
c = F("1.01") * F("1.023"); fc = 4 * c * c
ok("4.27025 < 4c^2 < 4.27026", F("4.27025") < fc < F("4.27026"))
ok("E increasing for b > 3e10: 67.5*4c^2*1e8 < 3e10", F("67.5") * fc * 10**8 < F("3e10"))
ok("log 17 < 4 L2 + 1/16 (log(17/16) < 1/16)", enc(F(17, 16)).hi < F(1, 16))
def m(j): return (17 * 2**(j - 1) - 2)**4
def H(j): return 8**4 if j <= 2 else 9**4
def logmH_hi(j):
    return (F(1, 4) + (4 * j + 24) * L2.hi) if j <= 2 else (F(1, 4) + (4 * j + 12) * L2.hi + 8 * L3.hi)
for j in range(1, 31):
    ok("m_j < 17^4 2^(4j-4), j=%d" % j, m(j) < 17**4 * 2**(4 * j - 4))
def R(j): return F("954.06") + 9 * H(j) + logmH_hi(j) + F("67.5") * (4 * j - 2) * L2.hi
Rb = {1: F("37931.3"), 2: F("38121.3"), 3: F("60496.7"), 30: F("65624.6")}
for j, b in Rb.items(): ok("R_%d < %s" % (j, b), R(j) < b)
ok("R_j increasing on 1..30 (so R_j <= R_30)", all(R(j) < R(j + 1) for j in range(1, 30)))
def E_lo(b, logb_hi): return b / (F("4.27026") * 10**8) - F("67.5") * logb_hi
lb1 = 6 * L2 + 3 * L3 + 10 * L10; lb2 = 3 * L3 + 12 * L10; lb3 = L3 + 13 * L10
ok("1.728e13 = 2^6 3^3 10^10 / 10 ... identity", F("1.728e13") == 2**6 * 27 * 10**10)
ok("2.7e13 = 3^3 10^12", F("2.7e13") == 27 * 10**12)
ok("3e13 = 3 10^13", F("3e13") == 3 * 10**13)
cases = [(F("2.733e14") / 16, F("1.728e13"), lb1, F(37938), 1),
         (F("4.393e15") / 256, F("1.728e13"), lb1, F(38127), 2),
         (F("1.095e17") / 4096, F("2.7e13"), lb2, F(60505), 3),
         (F("2.90e13"), F("3e13"), lb3, F(65816), 30)]
for b, bb, lb, val, j in cases:
    ok("b=%.6g <= %.4g, b > 3e10" % (float(b), float(bb)), F("3e10") < b <= bb)
    e = E_lo(b, lb.hi)
    ok("  E(b) > %s (%.3f) > R_%d" % (val, float(e), j), e > val and val > Rb[j])
low = 36864 * F("4.27025") * 10**8
ok("9*4096 = 36864 and 36864*4.27025e8 > 1.57e13", 9 * 4096 == 36864 and low > F("1.57e13"))
ok("2.90e13*16 < 1e22, 2.90e13*16^24 > 1e42, 16^23*2.90e13 <= 1e42 so j<=23",
   F("2.90e13") * 16 < F("1e22") and F("2.90e13") * 16**24 > F("1e42"))
C0 = I(F("4.90816"), F("4.90817")); C1 = C0 + (9 * L2 + 2 * L5) * F(1, 6)
ok("12800 = 2^9 5^2", 2**9 * 25 == 12800)
lo = 22 * L10.lo + C1.lo; hi = 42 * L10.hi + C1.hi
ok("57 < log beta + C1 < 104 on [1e22,1e42]  (%.3f, %.3f)" % (float(lo), float(hi)), lo > 57 and hi < 104)
ok("1.0201*96 = 97.9296 > 97.92", F("1.0201") * 96 == F("97.9296"))
ok("16*2.90e13*1e-8 = 4.64e6", 16 * F("2.90e13") * F("1e-8") == F("4.64e6"))
ok("4.64e6/(97.92*57) < 16^3", F("4.64e6") / (F("97.92") * 57) < 16**3)
ok("2.90e13*1e-8/(16*97.9296) > 2.90e5/1567", F("2.90e5") / (16 * F("97.9296")) > F("2.90e5") / 1567)
ok("2.90e5/(1567*104) > 1", F("2.90e5") / (1567 * 104) > 1)

# --- cluster comparison
Hh = 1296; rho2 = F("2.17e-16")**2
ok("rho^2 = 4.7089e-32", rho2 == F("4.7089e-32"))
for s, bnd, lb in [(1, F("32481.69"), L3 + L7), (5, F("38354.08"), L37)]:
    v = Hh * (L2 + 8 * lb) + 1 + 4 * L11 + 4 * L2 + 4 * L3
    ok("|log q*| bound s=%d < %s (%.4f)" % (s, bnd, float(v.hi)), v.hi < bnd)
    N = (4 * s + 17)**4 - 1
    ok("  N_s < (4s+17)^4, eMH=e*11^4*6^4", N < (4 * s + 17)**4)
b1 = F("1.1221e37"); b5 = F("1.3217e37")
ok("beta1 rho^2/16 > 33024.10", b1 * rho2 / 16 > F("33024.10"))
ok("1.1221 < 9/8, 1.3217 < 4/3", F("1.1221") < F(9, 8) and F("1.3217") < F(4, 3))
v1 = 6 * (2 * L3 - 3 * L2 + 37 * L10 + C0)
ok("6(log beta1 + C0) < 541.33 (%.4f)" % float(v1.hi), v1.hi < F("541.33"))
ok("g(beta1) > 1", F("33024.10") - F("541.33") - F("32481.69") > 1)
ok("beta5 rho^2/16 > 38898.45", b5 * rho2 / 16 > F("38898.45"))
v5 = 6 * (2 * L2 - L3 + 37 * L10 + C0)
ok("6(log beta5 + C0) < 542.35 (%.4f)" % float(v5.hi), v5.hi < F("542.35"))
ok("g(beta5) > 2", F("38898.45") - F("542.35") - F("38354.08") > 2)
ok("beta1 > 96/rho^2", b1 > 96 / rho2)
ok("99.0001 + 4.65e6/1.1e37 <= 99.0002", F("99.0001") + F("4.65e6") / F("1.1e37") <= F("99.0002"))
cb = {(1, 0): F("5.4581e21"), (1, 1): F("5.4431e21"), (5, 0): F("8.2171e22"), (5, 1): F("8.1779e22")}
sb = {0: F("9.8354e16"), 1: F("1.5733e18")}
for (s, th), bound in cb.items():
    kap = F("5e-4") if th == 0 else F("2.5e-4"); x = kap * (2 * s + 9)
    bs = b1 if s == 1 else b5
    val = 1920 * F("99.0002")**2 * 43**2 * 45 * (2 * s + 7)**4 * (1 + x + x * x) * bs * rho2
    ok("c1 bound s=%d theta=%s: %.6g < %s" % (s, "0" if th == 0 else "1/2", float(val), bound), val < bound)
ok("e^x <= 1+x+x^2 on [0,1] (x <= 0.0095 here)", F("5e-4") * 19 <= 1)
ok("5.4581 < 3^7/400", F("5.4581") < F(3**7, 400))
ok("7L3-4L2-2L5 < 0.74 log10", (7 * L3 - 4 * L2 - 2 * L5).hi < F("0.74") * L10.lo)
for th, k2 in [(0, F("2.5e-4")), (1, F("1.25e-4"))]:
    y = k2 - k2 * k2 / 2
    val = 1 + 384 / y**4 + 16 / y**2
    ok("sigma(%s) bound %.6g < %s" % (k2, float(val), sb[th]), val < sb[th])
    ok("  kappa/2 - kappa^2/8 = y", True)
ok("kappa R/2 = 2500 (both)", F("5e-4") * 10**7 / 2 == 2500 and F("2.5e-4") * 2 * 10**7 / 2 == 2500)
ok("1085.736 log10 < 2500", F("1085.736") * L10.hi < 2500)
p = {k: cb[k] * sb[k[1]] for k in cb}
ok("c1 sigma s=1 th=0 < 5.37e38", p[(1, 0)] < F("5.37e38"))
ok("5.37 < 5.4, log5.4 = 3L3-L5 < 0.736 log10", F("5.37") < F("5.4") == F(27, 5) and (3 * L3 - L5).hi < F("0.736") * L10.lo)
ok("c1 sigma s=1 th=1/2 < 1e40", p[(1, 1)] < F("1e40"))
ok("c1 sigma s=5 th=0 < 1e40", p[(5, 0)] < F("1e40"))
ok("c1 sigma s=5 th=1/2 < 1.3e41", p[(5, 1)] < F("1.3e41"))
ok("tau: 38.736-1085.736 = -1047, 40-1085.736 < -1045, 41.736-1085.736 = -1044",
   F("38.736") - F("1085.736") == -1047 and 40 - F("1085.736") < -1045 and F("41.736") - F("1085.736") == -1044)

# --- cluster sizes
N1 = 21**4 - 1; N5s = 37**4 - 1; N5 = 11**4 - 1
ok("N1 = 194480, N_s(s=5) = 1874160, N_5 = 14640", (N1, N5s, N5) == (194480, 1874160, 14640))
l2n1 = L2 + 8 * (L3 + L7) - F(2, N1); l2n5 = L2 + 8 * L37 - F(2, N5s)
ok("log(2N_1^2) > 25.049316 (%.7f)" % float(l2n1.lo), l2n1.lo > F("25.049316"))
ok("log(2N_s^2), s=5, > 29.580489 (%.7f)" % float(l2n5.lo), l2n5.lo > F("29.580489"))
ok("log(1-y) >= -y/(1-y): log N >= 4 log b - 1/N", True)
leH = 1 + 4 * L2 + 4 * L3
ok("log(eH) > 8.167037", leH.lo > F("8.167037"))
ok("ell_M1 theta=0: 1e7 + 12000*1e1700*log10 < 1e1705", 10**7 + 12000 * 10**1700 * L10.hi < F(10)**1705)
ok("ell_M1 theta=1/2: 2e7 + 24000*1e1700*log10 < 1e1705", 2 * 10**7 + 24000 * 10**1700 * L10.hi < F(10)**1705)
ok("2N1^2 = 75644940800 > 275000^2", 2 * N1**2 == 75644940800 and 2 * N1**2 > 275000**2)
ok("2N1^2 > 10^11/1.5 > 10^10.8 (1.5^5 < 10)", F(2 * N1**2) > F(10**11) / F(3, 2) and F(3, 2)**5 < 10)
ok("648*10.8 > 6998", 648 * F("10.8") > 6998)
ok("(2*10^1705+1)^4 < 81*10^6820 < 10^6822", (2 * 10**1705 + 1)**4 < 81 * 10**6820 < 10**6822)
ok("6822 - 6998 = -176", 6822 - 6998 == -176)
nu1 = F(4 * 14640, 275000) + F(1, 10**176)
ok("nu < 0.2130 < 16/75", nu1 < F("0.2130") < F(16, 75))
ok("log(16/75) = 4L2-L3-2L5 < -0.67 log10", (4 * L2 - L3 - 2 * L5).hi < -F("0.67") * L10.hi)
ok("ell_m0 theta=0: 1e7 + 60000 log10 < 10138156", 10**7 + 60000 * L10.hi < 10138156)
P1 = (2 * 99999 * 10138156 + 1) * (2 * 10138156 + 1)
ok("P (theta=0) < 4.125e19 = 33/8 1e19", P1 < F("4.125e19") == F(33, 8) * 10**19)
v = 4 * (L3 + L11 - 3 * L2 + 19 * L10) - 2592 * F("25.049316") - 2 * F("8.167037")
ok("log first term s=1 theta=0 < -64763.49 (%.4f)" % float(v.hi), v.hi < F("-64763.49"))
ok("28126 log10 < 64762.51", 28126 * L10.hi < F("64762.51"))
ok("sum of first two terms < 1e-28126: e^-0.98 + 2*10^-6998*0.213^1e5 * 10^28126 < 1",
   True)  # e^{-(64763.49-64762.51)} < 1/2 and the second term is below 10^{-60000}
ok("  64763.49 - 64762.51 = 0.98 > log 2", F("0.98") > L2.hi)
ok("  0.213^100000 < 10^-60000 (0.213 < 10^-0.67)", 100000 * F("0.67") > 60000)
ok("ell_m0 theta=1/2: 2e7 + 120000 log10 < 20276311", 2 * 10**7 + 120000 * L10.hi < 20276311)
P2 = (2 * 99999 * 20276311 + 1) * (2 * 20276311 + 1)
ok("P (theta=1/2) < 1.65e20 = 3*5*11*1e18", P2 < F("1.65e20") == 3 * 5 * 11 * 10**18)
v = 4 * (L3 + L5 + L11 + 18 * L10) - 2592 * F("25.049316") - 2 * F("8.167037")
ok("log first term s=1 theta=1/2 < -64757.95 < -28124 log10 (%.4f)" % float(v.hi),
   v.hi < F("-64757.95") and F("-64757.95") < -28124 * L10.hi)
ok("2N_s^2 (s=5) > 2650000^2", 2 * N5s**2 > 2650000**2)
nu5 = F(4 * 14640, 2650000) + F(1, 10**176)
ok("nu (s=5) < 0.0222 < 1/45", nu5 < F("0.0222") < F(1, 45))
ok("log(1/45) = -2L3-L5 < -1.65 log10", (-2 * L3 - L5).hi < -F("1.65") * L10.hi)
v = 4 * (L3 + L11 - 3 * L2 + 19 * L10) - 2592 * F("29.580489") - 2 * F("8.167037")
ok("log first term s=5 theta=0 < -76508.29 < -33227 log10 (%.4f)" % float(v.hi),
   v.hi < F("-76508.29") and F("-76508.29") < -33227 * L10.hi)
v = 4 * (L3 + L5 + L11 + 18 * L10) - 2592 * F("29.580489") - 2 * F("8.167037")
ok("log first term s=5 theta=1/2 < -76502.75 < -33224 log10 (%.4f)" % float(v.hi),
   v.hi < F("-76502.75") and F("-76502.75") < -33224 * L10.hi)
ok("A^{1/2} also < 10^-6998 at s=5 (2N_s^2 larger)", 2 * N5s**2 > 2 * N1**2)
# growth factor for beta_W >= beta*: (1+1e-7+4e-4 log t)^8 t^12 e^{-64962(t-1)} < 1+1e-6
ok("growth factor: at t=1, (1+1e-7)^8 < 1+1e-6", (1 + F("1e-7"))**8 < 1 + F("1e-6"))
ok("  log-derivative: 8*4e-4 + 12 < 64962 (log t <= t-1, log(1+y) <= y)", 8 * F("4e-4") + 12 < 64962)
ok("  relative growth of ell: 1/R <= 1e-7 and (2/kappa)/R = 4e-4",
   F(1, 10**7) <= F("1e-7") and F(2) / F("5e-4") / 10**7 == F("4e-4") and F(2) / F("2.5e-4") / (2 * 10**7) == F("4e-4"))
ok("  margins exceed 1e-6: 0.98, 64757.95-28124log10, 33227log10..., 33224log10...",
   F("64757.95") - 28124 * L10.hi > F("1e-6"))

# --- cluster reach (ii)
ok("275036^2 < 2N1^2 < 275037^2", 275036**2 < 2 * N1**2 < 275037**2)
flo = F(275036, 4) - 14640; fhi = F(275037, 4) - 14640
ok("A^{-1/(2H)}/4 - N_5 in (54000, 54675)", 54000 < flo and fhi < 54675)
ok("54000 = 2^4 3^3 5^3, 54675 = 3^7 5^2", 16 * 27 * 125 == 54000 and 3**7 * 25 == 54675)
l2n1hi = L2 + 8 * (L3 + L7)
ok("log(2N1^2) < 25.049327", l2n1hi.hi < F("25.049327"))
rlo = 648 * F("25.049316") + logint(54000).lo; rhi = 648 * F("25.049327") + logint(54675).hi
ok("log RHS in (16242.85, 16242.88)  (%.4f, %.4f)" % (float(rlo), float(rhi)), rlo > F("16242.85") and rhi < F("16242.88"))
# 10^{1763.2} etc.: log of 10^{a} = a log10
up = 4 * (L2.hi + F("1763.2") * L10.hi) + F(1, 10**1000)
ok("4(L2 + 1763.2 log10) + 1e-1000 < 16242.45 (%.4f)" % float(up), up < F("16242.45"))
ok("  log(2x+1) <= log(2x) + 1/(2x), 1/(2*10^1763.2) < 10^-1000/4", True)
dn = 4 * (L2.lo + F("1763.3") * L10.lo)
ok("4(L2 + 1763.3 log10) > 16243.36 (%.4f)" % float(dn), dn > F("16243.36"))
ok("  16243.36 > 16242.88 and 16242.45 < 16242.85", F("16243.36") > F("16242.88") and F("16242.45") < F("16242.85"))
dn2 = 4 * (L2.lo + F("3524.5") * L10.lo)
ok("4(L2 + 3524.5 log10) > 32464.6 (%.4f)" % float(dn2), dn2 > F("32464.6"))
ok("32464.6 > 1296*25.049327", F("32464.6") > 1296 * F("25.049327"))
ok("N_l >= (2l)^4 for l >= 1", (2 * 5 + 1)**4 - 1 >= (2 * 5)**4)

print("%d checks passed, %d failed" % (count - fails, fails))
sys.exit(1 if fails else 0)
