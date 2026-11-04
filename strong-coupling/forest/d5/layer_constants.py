#!/usr/bin/env python3
"""Explicit constants for the Hessian of the layer law in the link-logarithm
coordinates zeta (thm:layerbound of main.tex), from the constants of
prop:explicitdecay of main.tex.  mpmath at 40 digits; every printed value is
an upper bound rounded up in the last printed digit (the inputs are the
upper bounds of the paper, all formulas are monotone in them).

Part 1 reproduces the chain b, g0, h0, g1, y, sbar, pbar, g2, k1 of
prop:explicitdecay at theta = 0.197 (check against the paper: 0.862, 17.37,
14.96, 5.80, 7.41, 1966).
Part 2 computes the new constants at theta = 0 (lower bound) and 0.197
(weighted norm), and the least C_lambda for which m' >= lambda/60.
Part 3: the island-threshold arithmetic under the corrected calmness
hypothesis (cor:twoislands of main.tex)."""
from mpmath import mp, mpf, sqrt, exp, cosh, sinh, log, pi, e, findroot, ceil
mp.dps = 40

KAPPA   = mpf('1e-8')        # kappa_0 <= 1e-8 (hypothesis of lem:smallframes)
C3S     = mpf('1.29e-8')     # c_3 s   <= 1.29e-8   (from 4e5 c_3 (kappa+s) <= 1/195)
C3KP    = mpf('7.7e-8')      # c_3 kappa' <= 7.7e-8, kappa' = kappa_0 + 6 s
DIMG    = 3                  # SU(2)
H1      = mpf('17.001')      # ||(AB*)^{-1}|| in the norm sum_k ||z_k||
KPAPER  = mpf(2068)          # ||Hess_xi Hbar||_(0.197) <= 2068 lambda (proof of prop:explicitdecay)
C6_OVER_C3SQ = mpf('6.81e5') * DIMG   # c_6 = 6.81e5 dim G c_3^2

def up(x, digits=4):
    """round a positive mpf up to `digits` significant digits"""
    from mpmath import floor, log10
    ex = int(floor(log10(x)))
    sc = mpf(10)**(digits-1-ex)
    return ceil(x*sc)/sc

def chain(th, thpp=mpf('0.6'), need_g2=True):
    """the constants of prop:explicitdecay at weight theta = th"""
    eth = exp(th)
    b   = eth*sqrt(mpf('0.5')+mpf('1e-4'))
    epsA = mpf('7.35')*C3S
    AB  = sqrt(48)*eth*epsA
    sig = (sqrt(11)/8 + 260*KAPPA)/(mpf(1)/8 - 260*KAPPA)
    sbar = sqrt(1 + 2*sig**2*eth**2)
    g0  = 1/(mpf(1)/16 - 700*KAPPA - (cosh(2*th)-1)*(mpf(1)/16 + 454*KAPPA))
    E   = AB*b
    h0  = g0/(1 - g0*E)
    g1  = g0/(1 - g0*(2*E + AB**2))
    y   = b*h0
    pbar = sbar*(1 + y*AB)
    if th > 0 and need_g2:
        d = chain(thpp, need_g2=False)                       # ||M||_(theta'') with M = Y(A-B)S
        Mpp = d['y']*d['AB']*d['sbar']
        g2 = 1/(1 - 2*sig**2*sinh(th)**2 - 2*(th/thpp)**2*Mpp**2)
    else:
        g2 = mpf(1)
    x1 = 2*(2+mpf('1e-4'))*H1*eth**2*pbar**2
    k1 = 48*eth*pbar**2 + 12*C3KP*x1
    return dict(th=th, eth=eth, b=b, epsA=epsA, AB=AB, sig=sig, sbar=sbar, g0=g0,
                h0=h0, g1=g1, y=y, pbar=pbar, g2=g2, x1=x1, k1=k1)

print("=== Part 1: the chain of prop:explicitdecay at theta = 0.197 (paper: 0.862, 17.37, 14.96, 5.80, 7.41, 1966)")
d = chain(mpf('0.197'))
for k in ('b','g0','h0','g1','y','sbar','pbar','g2','k1','x1'):
    print(f"  {k:5s} <= {up(d[k],5)}")

print("\n=== Part 2: constants of the zeta-coordinate Hessian (thm:layerbound, Step 4)")
# l^1 norms (theta-independent)
S11  = 1 + 2*mpf('8.01')*(mpf(1)/8 + 260*KAPPA)      # ||S||_{1->1}
Y11  = (2+mpf('1e-4'))*H1                           # ||Y||_{1->1}
epsA = mpf('7.35')*C3S
AB11 = 2*epsA                                       # ||A-B||_{1->1}
q1   = S11*(1 + Y11*AB11)                           # ||Pbar||_{1->1}
dg11 = Y11*AB11*S11                                 # ||Dg||_{1->1}
j1   = 1/(1-dg11)                                   # ||J||_{1->1}
print(f"  ||S||_1->1 <= {up(S11,6)}, ||Pbar||_1->1 = q1 <= {up(q1,6)}, ||Dg||_1->1 <= {up(dg11,3)}, j1 <= {up(j1,8)}")

def zeta_terms(th):
    d  = chain(th)
    eth, y, pbar, x1 = d['eth'], d['y'], d['pbar'], d['x1']
    dg = y*d['AB']*d['sbar']           # ||Dg||_(theta)
    j  = 1/(1-dg)                      # ||J||_(theta)
    n1 = 12*q1                         # ||grad_xi Hbar||_inf <= lambda*(c3 kappa')*n1/c3 + n0 c3
    n0 = (mpf('0.84') + 30 + 6*y)*q1   # lambda-independent part of ||grad Hbar||_inf, times c_3
    # middle term  <grad Hbar, D^2 xi[p,q]>, per |u||u'|:
    mid_lam = 12*C3KP*q1*j1*x1*j**2    # times lambda
    mid_c3  = n0*j1*x1*j**2            # times c_3^2
    # third term (Hessian of log det(1+Dg) o xi), per |u||u'|, times c_3^2:
    W3a = 2*DIMG*eth**3*pbar**4*y**2*j**4
    W3b = 4*DIMG*eth**3*pbar**3*y**2*j*j**2
    W3c = 2*DIMG*eth**2*pbar**3*y*j*j**2          # (times c_3 only; bounded by c_3^2)
    W3d = 2*DIMG*pbar*y*j*x1*j**2
    W3e = 2*DIMG*j*y*pbar*q1*j1*x1*j**2
    third_c3 = W3a+W3b+W3c+W3d+W3e
    return dict(d=d, dg=dg, j=j, n1=n1, n0=n0, mid_lam=mid_lam, mid_c3=mid_c3,
                W3=(W3a,W3b,W3c,W3d,W3e), third_c3=third_c3)

z0 = zeta_terms(mpf(0))
print("  theta = 0:")
print(f"    ||Dg|| <= {up(z0['dg'],3)},  j0 = ||J|| <= {up(z0['j'],8)},  y0 <= {up(z0['d']['y'],5)},  pbar0 <= {up(z0['d']['pbar'],5)},  x1(0) <= {up(z0['d']['x1'],5)}")
print(f"    n0 <= {up(z0['n0'],4)} (x c_3),  12 q1 (c3 kappa') <= {up(12*q1*C3KP,4)} (x lambda / c_3)")
print(f"    middle term: {up(z0['mid_lam'],4)} lambda + {up(z0['mid_c3'],4)} c_3^2")
print(f"    third term : W3a..W3e = " + ", ".join(f"{up(w,4)}" for w in z0['W3']) + f"  sum = {up(z0['third_c3'],4)} c_3^2")
lam_indep0 = z0['mid_c3'] + z0['third_c3']
print(f"    lambda-independent total (besides c_6): {up(lam_indep0,4)} c_3^2 = {up(lam_indep0/C6_OVER_C3SQ,4)} c_6")

def mprime_over_lambda(Cl):
    """lower bound of Hess_zeta U / lambda under lambda >= Cl c_6"""
    z = z0
    lb = (mpf(1)/30 - 1/mpf(Cl))/(1+z['dg'])**2
    return lb - z['mid_lam'] - lam_indep0/(mpf(Cl)*C6_OVER_C3SQ)

for Cl in (60, 120, 180, 200, 220, 240, 300, 400):
    v = mprime_over_lambda(Cl)
    print(f"    C_lambda = {Cl:4d}:  m'/lambda >= {float(v):+.5f}   (1/60 = {1/60:.5f})  {'OK' if v >= mpf(1)/60 else '--'}")
# least integer C_lambda with m' >= lambda/60
Cl = 60
while mprime_over_lambda(Cl) < mpf(1)/60: Cl += 1
print(f"    least integer C_lambda with m' >= lambda/60: {Cl}")

CL = 240
zw = zeta_terms(mpf('0.197'))
print("  theta = 0.197:")
print(f"    ||Dg||_(th) <= {up(zw['dg'],3)},  j_th <= {up(zw['j'],8)}")
print(f"    middle term: {up(zw['mid_lam'],4)} lambda + {up(zw['mid_c3'],4)} c_3^2")
print(f"    third term : W3a..W3e = " + ", ".join(f"{up(w,4)}" for w in zw['W3']) + f"  sum = {up(zw['third_c3'],4)} c_3^2")
lam_indep_w = zw['mid_c3'] + zw['third_c3']
Kp = zw['j']**2*KPAPER + zw['mid_lam'] + lam_indep_w/(CL*C6_OVER_C3SQ)
print(f"    with lambda >= {CL} c_6:  K'/lambda <= {up(Kp,5)}")
alpha_max = mpf('0.197')/sqrt(240*Kp)
print(f"    admissible alpha <= theta'(lambda/(240 K'))^(1/2) = {up(alpha_max,4)}  (need >= 2.5e-4: {'OK' if alpha_max >= mpf('2.5e-4') else 'FAIL'})")
print(f"    with lambda >= 60 c_6 instead: K'/lambda <= {up(zw['j']**2*KPAPER + zw['mid_lam'] + lam_indep_w/(60*C6_OVER_C3SQ),5)}")

print("\n=== Part 2b: chart inclusion (thm:layerbound, Step 2) and the size of h")
c3r = C3S/32
t   = mpf('31.5')            # ||X||_inf <= 31.5 r for zeta in D_{1.25 r}
zk  = 17*12*c3r*t**2         # max_k ||z_k|| / r : 204 (c_3 r) t^2 r
print(f"  c_3 r <= {up(c3r,3)};  max_k||z_k|| <= {up(zk,3)} r;  ||eta||_inf <= {up(t*(1+ mpf('0.88')*17*12*c3r*t),6)} r < 32 r;  ||h||_inf <= {up(mpf('0.88')*zk,3)} r")

print("\n=== Part 3: island threshold arithmetic (cor:twoislands)")
C0 = mpf(1)/2 + mpf(2)/3*log(24*pi**3)
M, H = mpf(11)**4, mpf(6)**4
def beta_star(rho, s, join):
    Dn = (2*join+1)**4 - 1
    logq = -H*(log(2)+2*log(Dn)) - log(e*M*H)
    f = lambda b: 6*(log(b)+C0) - b*rho**2/16 - logq
    return findroot(f, 16*abs(logq)/rho**2*mpf('1.02')), logq
for rho, s in ((mpf('1e-8'),1),(mpf('1e-8'),5),(mpf('1e-9'),1)):
    b, lq = beta_star(rho, s, 2*s+1)
    print(f"  check paper: rho={float(rho):.0e} s={s}: |log q*|={float(-lq):.1f} beta*={float(b):.4e}")
# calmness: kappa_0 <= eta + 17 r, s = 32 r, r = 43 rho, eta = 92 rho:
coef = 92 + 17*43 + 32*43
print(f"  kappa_0 + s <= {coef} rho;  SU(2) chart hypothesis kappa_0 + s <= 6.0e-13  =>  rho <= {up(mpf('6.0e-13')/coef,3)}")
rho_c = mpf('6.0e-13')/coef
for s in (1,5):
    b, lq = beta_star(rho_c, s, 2*s+8)
    print(f"  rho={float(rho_c):.3e} s={s} join={2*s+8}: |log q*|={float(-lq):.1f} beta*={float(b):.4e}")
print(f"  beta_W >= 4 C_lambda c_6 = {float(4*CL*mpf('9.3e14')):.3e}  (C_lambda = {CL}, c_6 <= 9.3e14 for SU(2))")
