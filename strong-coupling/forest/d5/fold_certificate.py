"""Rigorous (interval-arithmetic) verification for the fold obstruction,
Proposition fold.  Everything printed as 'CHECK ... OK' is a verified inequality.
Run: python3 fold_certificate.py   (about 1-3 minutes)
"""
import time
from fold_intervals import *
t0 = time.time()
SQ2 = iv.sqrt(iv.mpf(2))
ALPHA = 2 * (iv.pi / 2 - iatan(iv.sqrt(iv.mpf(7))))      # alpha = acos(3/4) = 2 asin(1/(2 sqrt2)); pi/2 - atan(sqrt 7) = asin(1/(2sqrt2))? no: use acos(3/4) = atan(sqrt7/3)
ALPHA = iatan(iv.sqrt(iv.mpf(7)) / 3)                      # acos(3/4) = atan(sqrt(7)/3)
print("alpha = acos(3/4) =", ALPHA)
Q_LO = I('0.49') - ALPHA; Q_HI = I('-0.5') + ALPHA
print("q_lo = 0.49 - alpha =", Q_LO, "  q_hi = -1/2 + alpha =", Q_HI)
assert lo(Q_LO) > -0.233 and hi(Q_HI) < 0.223
print("CHECK [q_lo, q_hi] subset of [-0.233, 0.223]: OK")

# ---------- T1: the torus function f_0 --------------------------------------
SUM = I('2.44')             # sum of the seven angles: -1/2 + 6*0.49
def chiR(q):  # maximiser theta_1 = -1/2 (q >= -0.005)
    return 2 - 4 * SQ2 * iv.sin((q + I('0.5')) / 2)
def chiL(q):  # maximiser the six links at 0.49 (q <= -0.005)
    return 2 - 4 * SQ2 * iv.sin((I('0.49') - q) / 2)
def fR(q): return q + chiR(q) * (SUM - 7 * q) / 8
def fL(q): return q + chiL(q) * (SUM - 7 * q) / 8
def fRp(q): return 1 - 7 * chiR(q) / 8 - 2 * SQ2 * iv.cos((q + I('0.5')) / 2) * (SUM - 7 * q) / 8
def fLp(q): return 1 - 7 * chiL(q) / 8 + 2 * SQ2 * iv.cos((I('0.49') - q) / 2) * (SUM - 7 * q) / 8
def fRpp(q):
    chip = -2 * SQ2 * iv.cos((q + I('0.5')) / 2); chipp = SQ2 * iv.sin((q + I('0.5')) / 2)
    return chipp * (SUM - 7 * q) / 8 - 7 * chip / 4

# consistency of the explicit formula with the quaternion code on the torus
for qq in ['-0.2', '-0.05', '0.0', '0.1', '0.2']:
    V, chi, m = V_of(I(qq), D(0))
    ang = iatan(V[1].v / V[0].v)           # V in T, V_0 > 0 here
    f = fR(I(qq)) if mpf(qq) >= -0.005 else fL(I(qq))
    assert lo(ang) <= hi(f) and lo(f) <= hi(ang), (qq, ang, f)
print("CHECK explicit f_0 agrees with the quaternion V on the torus: OK")

# (a) maximiser switch at q_s = -0.005: q+1/2 = 0.49-q
# (b) convexity on [-0.005, q_hi]: by hand (signs); interval check of f'' > 0 anyway
qs = [mpf(-0.005) + (hi(Q_HI) + mpf(0.005)) * k / 40 for k in range(41)]
assert all(lo(fRpp(iv.mpf([qs[k], qs[k + 1]]))) > 0 for k in range(40))
print("CHECK f_0'' > 0 on [-0.005, q_hi]: OK")
# (c) bracket of the minimiser q*
a, b = I('0.0741'), I('0.0742')
assert hi(fRp(a)) < 0 and lo(fRp(b)) > 0
print("CHECK f_0'(0.0741) < 0 < f_0'(0.0742): OK,  f_0'(0.0741) in", fRp(a), " f_0'(0.0742) in", fRp(b))
gam_lo = lo(fR(a) + fRp(a) * (b - a)); gam_hi = min(hi(fR(a)), hi(fR(b)))
GAM = iv.mpf([gam_lo, gam_hi])
print("gamma* = f_0(q*) in", GAM, "  (by convexity: f_0(q) >= f_0(a) + f_0'(a)(q-a))")
# (d) f_0'(-0.005+) < 0 and f_0'(q_hi) > 0 (so q* is interior); f_0(-0.005) > gamma*
assert hi(fRp(I('-0.005'))) < 0 and lo(fRp(Q_HI)) > 0 and lo(fR(I('-0.005'))) > gam_hi
print("CHECK f_0'(-0.005) in", fRp(I('-0.005')), "< 0;  f_0'(q_hi) in", fRp(Q_HI), "> 0;  f_0(-0.005) in", fR(I('-0.005')), "> gamma*: OK")
# (e) f_0 increasing on [q_lo, -0.005]: f' >= 1 - 7 chi/8 with chi <= chi(-0.005)
chimax = chiL(I('-0.005'))
assert hi(chimax) < 0.6146 and lo(1 - 7 * chimax / 8) > 0.46
print("CHECK chi(-0.005) in", chimax, "<= 0.6146, so f_0' > 0.46 on [q_lo,-0.005]: OK")
# (f) q_L(lambda) for lambda in Lambda = [gamma* - 1e-3, gamma* + 1e-3]
ETA = mpf('0.001')
qL1, qL2 = I('-0.0172'), I('-0.0156')
assert hi(fL(qL1)) < gam_lo - ETA and lo(fL(qL2)) > gam_hi + ETA
print("CHECK f_0(-0.0172) in", fL(qL1), "< gamma*-1e-3 and f_0(-0.0156) in", fL(qL2), "> gamma*+1e-3: OK")
print("  so q_L(lambda) in [-0.0172, -0.0156] for every lambda in Lambda, and q* - q_L >= 0.0741+0.0156 = 0.0897")
GAP = 2 * SQ2 * iv.sin(I('0.0897') / 2)
print("  chordal distance delta(g*, g_L(lambda)) >= 2 sqrt2 sin(0.0897/2) >=", lo(GAP))
assert lo(GAP) > 0.1268
# (g) f_0'(q_L) >= 1.3 (actual) -- use the explicit derivative on the bracket
fpL = fLp(iv.mpf([lo(qL1), hi(qL2)]))
print("  f_0'(q) on [-0.0172,-0.0156] in", fpL)
assert lo(fpL) > 1.3
# (h) the near-side roots q_+- lie in [q*-0.021, q*+0.021] subset (-0.005, q_hi): f'' >= 4.6
assert lo(fRpp(iv.mpf([mpf('0.05'), mpf('0.1')]))) > 4.6
print("CHECK f_0'' > 4.6 on [0.05, 0.1]: OK  (so |q_+- - q*| <= sqrt(2e-3/4.6) < 0.021 for lambda in Lambda)")

# ---------- T2: the strip s in [0, s1]: dV_k/ds > 0 -----------------------------
S1 = mpf('0.02')
def strip(q0, q1, depth=0):
    V, chi, m = V_of(iv.mpf([q0, q1]), D(iv.mpf([0, S1]), 1))
    dk = V[3].d; dj = V[2].d
    if lo(dk) > 0:
        return [(lo(dk), hi(dk), hi(abs(dj)))], 1
    if depth > 30: raise RuntimeError("strip: no positive enclosure at %s" % ((q0, q1),))
    mid = (q0 + q1) / 2
    r1, n1 = strip(q0, mid, depth + 1); r2, n2 = strip(mid, q1, depth + 1)
    return r1 + r2, n1 + n2
res, nbox = strip(mpf('-0.233'), mpf('0.223'))
print("CHECK strip q in [-0.233,0.223], s in [0,%s]: dV_k/ds > 0 on every box (%d boxes): OK" % (S1, nbox))
# finer enclosures on the two q-ranges that enter the constants
def fine(q0, q1, n):
    out = []
    for k in range(n):
        a = q0 + (q1 - q0) * k / n; b = q0 + (q1 - q0) * (k + 1) / n
        V, chi, m = V_of(iv.mpf([a, b]), D(iv.mpf([0, S1]), 1))
        out.append((lo(V[3].d), hi(V[3].d), hi(abs(V[2].d))))
    return min(r[0] for r in out), max(r[1] for r in out), max(r[2] for r in out)
KAPPA, _, _ = fine(mpf('-0.0172'), mpf('-0.0156'), 16)
_, LMAX, DJMAX = fine(mpf('0.05'), mpf('0.1'), 50)
print("  on q in [-0.0172,-0.0156] (the atom g_L): dV_k/ds >= kappa =", KAPPA)
print("  on q in [0.05,0.1] (the atoms g_+-): dV_k/ds <=", LMAX, ", |dV_j/ds| <=", DJMAX)
assert KAPPA > 0.49

# ---------- T3: the region s in [s1, alpha]: V_k > 0 ----------------------------
count = [0]
def region(q0, q1, s0, s1, depth=0):
    count[0] += 1
    V, chi, m = V_of(iv.mpf([q0, q1]), D(iv.mpf([s0, s1])))
    if lo(V[3].v) > 0: return
    if depth > 40: raise RuntimeError("region: cannot exclude box %s" % ((q0, q1, s0, s1),))
    if (q1 - q0) > (s1 - s0):
        mid = (q0 + q1) / 2; region(q0, mid, s0, s1, depth + 1); region(mid, q1, s0, s1, depth + 1)
    else:
        mid = (s0 + s1) / 2; region(q0, q1, s0, mid, depth + 1); region(q0, q1, mid, s1, depth + 1)
region(mpf('-0.233'), mpf('0.223'), S1, hi(ALPHA))
print("CHECK region q in [-0.233,0.223], s in [%s, alpha]: V_k > 0 on every box (%d boxes examined): OK" % (S1, count[0]))
print("Hence V^{-1}(T) = T: every path link g with V(g) in the torus lies in the torus.")

# ---------- constants for the W_1 jump ----------------------------------------
wL_max = 1 / (lo(fpL) * KAPPA ** 2)
print("weight of the atom g_L:  w_L <= 1/(f_0'(q_L) kappa^2) <=", wL_max, "(times e^{-beta S} <= 1)")
print("transverse factor at any torus point of the strip: |lambda_perp| = |dV/ds(q,0)| <= sqrt(LMAX^2+DJMAX^2) <=", mp.sqrt(LMAX**2 + DJMAX**2))
print("done in %.1f s" % (time.time() - t0))
