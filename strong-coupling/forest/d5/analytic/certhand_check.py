"""Re-check of the hand arguments in patches P1 (Lemma foldtransverse,
Proposition fold Steps 1-3) and P2 (Theorem sfprop with the constant 1/52).
Exact rationals for every inequality written in the text; mpmath.iv for the
Step 1 evaluations; a high-precision numerical test of the closed form of
Lemma foldtransverse (the lemma itself is proved by algebra in the text).
The text does not depend on this script."""
from fractions import Fraction as F
from mpmath import mp, iv, mpf
import mpmath, random
mp.dps = 50; iv.dps = 50
ok = 0; bad = 0
def chk(name, cond):
    global ok, bad
    print(('PASS ' if cond else 'FAIL ') + name)
    if cond: ok += 1
    else: bad += 1
PI_UP = F(31416, 10000)
SQ2_LO, SQ2_HI = F(141421356, 10**8), F(141421357, 10**8)
chk('pi < 3.1416 (3.14159265 < 3.1416)', F(314159266, 10**8) < PI_UP)
chk('1.41421356^2 < 2 < 1.41421357^2', SQ2_LO**2 < 2 < SQ2_HI**2)
chk('5.6568^2 < 32 < 5.65686^2', F(56568, 10**4)**2 < 32 < F(565686, 10**5)**2)
# ---------- Step 2 of Proposition fold ----------
x = F(99, 400)
s_lo = x - x**3/6
chi_up = 2 - F(56568, 10**4)*s_lo
chk('Step2: chi <= 2-4sqrt2 sin(99/400) < 0.6143', chi_up < F(6143, 10**4))
Y_up = F(6143, 10**4)*7*PI_UP/32
chk('Step2: 0.6143*7pi/32 < 0.4223', Y_up < F(4223, 10**4))
chk('Step2: 1-0.4223^2/2 > 0.9108', 1 - F(4223, 10**4)**2/2 > F(9108, 10**4))
chk('Step2: 7/8*0.6143 < 0.5376', F(7, 8)*F(6143, 10**4) < F(5376, 10**4))
chk('Step2: 0.9108-0.5376 > 1/3', F(9108, 10**4) - F(5376, 10**4) > F(1, 3))
chk('Step2: cos(pi/4)=1/sqrt2 < 3/4', F(3, 4)**2 > F(1, 2))
# ---------- Step 3 of Proposition fold ----------
dlo, dhi = F(5056, 10**4)/2, F(5072, 10**4)/2      # d/2
chi_hi = 2 - F(56568, 10**4)*(dlo - dlo**3/6)
chi_lo = 2 - F(565686, 10**5)*dhi
chk('Step3: chi/8 < 0.07315', chi_hi/8 < F(7315, 10**5))
chk('Step3: chi/8 > 0.07067', chi_lo/8 > F(7067, 10**5))
far, near = F(4828, 10**4), F(5056, 10**4)
S2 = far**2 + 6*near**2; S4 = far**4 + 6*near**4
chk('Step3: sum phi^2 > 1.7668', S2 > F(17668, 10**4))
chk('Step3: sum phi^4 > 0.4464', S4 > F(4464, 10**4))
chk('Step3: 1.7668/3 > 0.5889 and 0.4464/45 > 0.0099',
    F(17668, 10**4)/3 > F(5889, 10**4) and F(4464, 10**4)/45 > F(99, 10**4))
chk('Step3: 0.07315(7-0.5889-0.0099) < 0.4683',
    F(7315, 10**5)*(7 - F(5889, 10**4) - F(99, 10**4)) < F(4683, 10**4))
D_lo = F(7067, 10**5)*(F(244, 100) + 7*F(156, 10**4))
D_hi = F(7315, 10**5)*(F(244, 100) + 7*F(172, 10**4))
chk('Step3: Delta in [0.1801,0.1873]', D_lo >= F(1801, 10**4) and D_hi <= F(1873, 10**4))
chk('Step3: 1-0.1873^2/2 > 0.9824', 1 - F(1873, 10**4)**2/2 > F(9824, 10**4))
chk('Step3: 1-0.1801^2/6+0.1873^4/120 < 0.99461',
    1 - F(1801, 10**4)**2/6 + F(1873, 10**4)**4/120 < F(99461, 10**5))
chk('Step3: 0.9824-0.99461*0.4683 > 0.5166',
    F(9824, 10**4) - F(99461, 10**5)*F(4683, 10**4) > F(5166, 10**4))
chk('Step3: 1/(1.36*0.5166^2) < 2.78', 1/(F(136, 100)*F(5166, 10**4)**2) < F(278, 100))
# partial-fraction bound x cot x <= 1 - x^2/3 - x^4/45, numerical sanity
chk('x cot x <= 1-x^2/3-x^4/45 on a grid of (0,3)',
    all(mpf(t)*mpmath.cot(mpf(t)) <= 1 - mpf(t)**2/3 - mpf(t)**4/45
        for t in [k/100 for k in range(1, 300)]))
# ---------- Lemma foldtransverse: closed form, numerical test ----------
def qmul(a, b):
    a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
    return [a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
            a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0]
def qlog(h):
    n = mpmath.sqrt(h[1]**2+h[2]**2+h[3]**2)
    if n == 0: return [mpf(0)]*3
    p = mpmath.atan2(n, h[0]); return [p*h[k]/n for k in (1, 2, 3)]
def qexp(Y):
    n = mpmath.sqrt(sum(y*y for y in Y))
    if n == 0: return [mpf(1), 0, 0, 0]
    return [mpmath.cos(n)] + [mpmath.sin(n)*y/n for y in Y]
def phi(t): return min(mpf(1), max(mpf(0), 2-2*t))
def V_and_Xi(thetas, ws, q, s):
    g = [mpmath.cos(s)*mpmath.cos(q), mpmath.cos(s)*mpmath.sin(q), mpf(0), mpmath.sin(s)]
    gi = [g[0], -g[1], -g[2], -g[3]]
    hs = [qmul(gi, [mpmath.cos(t), mpmath.sin(t), 0, 0]) for t in thetas]
    m = max(2*mpmath.sqrt(1-h[0]) for h in hs)
    chi = phi(m)
    if chi == 0: return qmul(g, [1, 0, 0, 0]), mpf(1)
    Y = [chi*sum(w*qlog(h)[k] for w, h in zip(ws, hs)) for k in range(3)]
    n = mpmath.sqrt(sum(y*y for y in Y))
    sig = mpmath.sin(n)/n if n > 0 else mpf(1)
    psis = [mpmath.acos(mpmath.cos(s)*mpmath.cos(t-q)) for t in thetas]
    A = chi*sum(w*(p*mpmath.cot(p) if p > 0 else 1) for w, p in zip(ws, psis))
    return qmul(g, qexp(Y)), mpmath.cos(n) - sig*A
random.seed(1); err = mpf(0); xi_min = mpf(9); xi_max = mpf(-9)
th7 = [mpf(-1)/2] + [mpf(49)/100]*6
for it in range(3000):
    if it < 2000: th, ws = th7, [mpf(1)/8]*7
    else:
        th = [mpf(random.uniform(-0.6, 0.6)) for _ in range(14)]; ws = [mpf(1)/16]*14
    q = mpf(random.uniform(-0.4, 0.4)); s = mpf(random.uniform(-1.6, 1.6))
    v, xi = V_and_Xi(th, ws, q, s)
    err = max(err, abs(v[2]), abs(v[3]-mpmath.sin(s)*xi))
    if it < 2000: xi_min = min(xi_min, xi); xi_max = max(xi_max, xi)
chk('Lemma: |V_j| and |V_k - sin s Xi| < 1e-40 at 3000 random points (err %.1e)' % float(err), err < mpf('1e-40'))
chk('Lemma/Step2: 1/3 < Xi <= 1 at 2000 random points of the fold configuration (min %.4f)' % float(xi_min),
    xi_min > mpf(1)/3 and xi_max <= 1)
lam = [V_and_Xi(th7, [mpf(1)/8]*7, mpf(-0.0172) + k*mpf(0.0016)/100, mpf(0))[1] for k in range(101)]
chk('Step3: lambda_perp on [-0.0172,-0.0156] grid min %.5f > 0.5166' % float(min(lam)), min(lam) > mpf('0.5166'))
# ---------- Step 1 evaluations (interval arithmetic re-check) ----------
def chi_iv(q):
    d = iv.mpf(49)/100 - q if q.b <= -0.005 else q + iv.mpf(1)/2
    return 2 - 4*iv.sqrt(2)*iv.sin(d/2)
def f0(q): return q + chi_iv(q)/8*(iv.mpf(244)/100 - 7*q)
def f0p_far(q):   # q >= -0.005: d = q+1/2, chi' = -2sqrt2 cos(d/2)
    c = 2 - 4*iv.sqrt(2)*iv.sin((q+iv.mpf(1)/2)/2); cp = -2*iv.sqrt(2)*iv.cos((q+iv.mpf(1)/2)/2)
    return 1 + cp/8*(iv.mpf(244)/100-7*q) - 7*c/8
def f0p_near(q):  # q <= -0.005: d = 49/100-q, chi' = +2sqrt2 cos(d/2)
    d = iv.mpf(49)/100 - q
    c = 2 - 4*iv.sqrt(2)*iv.sin(d/2); cp = 2*iv.sqrt(2)*iv.cos(d/2)
    return 1 + cp/8*(iv.mpf(244)/100-7*q) - 7*c/8
P = lambda s: iv.mpf(s)
chk('Step1: alpha=arccos(3/4) in (0.7227,0.7228)', iv.cos(P('0.7227')).a > 0.75 and iv.cos(P('0.7228')).b < 0.75)
chk("Step1: f0'(0.0741)<0<f0'(0.0742)", f0p_far(P('0.0741')).b < 0 and f0p_far(P('0.0742')).a > 0)
lo = f0(P('0.0741')) + P('1e-4')*f0p_far(P('0.0741')); hi1 = f0(P('0.0741')); hi2 = f0(P('0.0742'))
chk('Step1: gamma* bracket inside [0.16978261,0.16978263]', lo.a >= P('0.16978261').b and min(hi1.b, hi2.b) <= P('0.16978263').a)
g_hi = min(hi1.b, hi2.b); g_lo = lo.a
chk('Step1: f0(-0.0172) < gamma*-1e-3', f0(P('-0.0172')).b < g_lo - 1e-3)
chk('Step1: f0(-0.0156) > gamma*+1e-3', f0(P('-0.0156')).a > g_hi + 1e-3)
chk('Step1: f0(0.05), f0(0.1) > gamma*+1e-3', f0(P('0.05')).a > g_hi + 1e-3 and f0(P('0.1')).a > g_hi + 1e-3)
chk("Step1: f0'(-0.005+) < 0 and f0'(q_hi) > 0", f0p_far(P('-0.005')).b < 0 and f0p_far(P('0.2227')).a > 0)
s2 = iv.sqrt(2)
f2_lo = (s2*iv.sin((P('0.05')+P('0.5'))/2))/8*(P('2.44')-P('0.7')) + P(7)/4*2*s2*iv.cos((P('0.1')+P('0.5'))/2)
chk("Step1: endpoint-factor lower bound of f0'' on [0.05,0.1] > 4.8 (%.4f)" % float(f2_lo.a), f2_lo.a > 4.8)
Q = iv.mpf(['-0.0172', '-0.0156'])
fp = f0p_near(Q)
chk("Step1: f0' in [1.36,1.37] on [-0.0172,-0.0156] (naive interval %s)" % fp, fp.a >= 1.36 and fp.b <= 1.37)
chk('Taylor remainders: 0.73^11/11! < 1e-9 and 0.73^12/12! < 1e-9',
    F(73, 100)**11/39916800 < F(1, 10**9) and F(73, 100)**12/479001600 < F(1, 10**9))
# ---------- Theorem sfprop with c = 1/52 ----------
Sbar_hi = F(22740, 100)
chk('sfprop: 168+42sqrt2 < 227.40', 168 + 42*SQ2_HI < Sbar_hi)
chk('sfprop: cosh(0.012) <= 1/(1-0.012^2/2) < 1.0001', 1/(1 - F(12, 1000)**2/2) < F(10001, 10**4))
chk('sfprop: 0.012^2/2*1.0001*227.40 < 0.01638 < 1/52',
    F(12, 1000)**2/2*F(10001, 10**4)*Sbar_hi < F(1638, 10**5) < F(1, 52))
chk('sfprop: cosh(0.006) < 1.0001', 1/(1 - F(6, 1000)**2/2) < F(10001, 10**4))
chk('sfprop: 0.006^2/2*1.0001*227.40 < 0.004094', F(6, 1000)**2/2*F(10001, 10**4)*Sbar_hi < F(4094, 10**6))
c = F(1, 52) - F(4094, 10**6)
chk('sfprop: 1/52-0.004094 > 0.0151', c > F(151, 10**4))
chk('sfprop: e^0.004 <= 1/(1-0.004) < 1.00402, e^0.008 <= 1/(1-0.008) < 1.00807',
    1/(1 - F(4, 1000)) < F(100402, 10**5) and 1/(1 - F(8, 1000)) < F(100807, 10**5))
p1sq_hi = F(49, 8) + F(49, 8)/SQ2_LO + F(49, 64)
chk('sfprop: p1^2 < 11.222', p1sq_hi < F(11222, 1000))
bbar_hi = 168 + 336*SQ2_HI
chk('sfprop: bbar^2 < 413681', bbar_hi**2 < 413681)
cc = F(151, 10**4)
chk('sfprop: 1/0.0151 < 66.3', 1/cc < F(663, 10))
chk('sfprop: 1.00402*11.222/0.0151 < 746.2', F(100402, 10**5)*F(11222, 1000)/cc < F(7462, 10))
chk('sfprop: 1.00807*413681/0.0151 < 2.762e7', F(100807, 10**5)*413681/cc < F(2762, 1000)*10**7)
print('%d checks passed, %d failed' % (ok, bad))
