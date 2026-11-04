"""Certified SU(2) constants of prop:smallframessharp (mpmath.iv).  Usage: python3 witness_radius.py
Conventions: angle(g) in [0,pi] with eigenvalues e^{+-i angle}; delta(g,I)=2*sqrt2*sin(angle/2);
||log g||=sqrt2*angle; p = psi_0(kappa) = 2*arcsin(kappa/(2*sqrt2)), the largest plaquette angle."""
from mpmath import iv, mp
iv.dps = 30; mp.dps = 30
I = iv.mpf; s2 = iv.sqrt(I(2)); pi = iv.pi
def asin(x):
    """certified enclosure of arcsin on [0,1): checked through iv.sin (increasing)."""
    x = I(x); eps = mp.mpf(10)**-25
    lo = mp.asin(x.a) - eps; hi = mp.asin(x.b) + eps
    assert iv.sin(I(lo)).b <= x.a and iv.sin(I(hi)).a >= x.b and hi < mp.pi/2
    return I([lo, hi])
iv.asin = asin
def ang(kappa): return 2*iv.asin(I(kappa)/(2*s2))
def show(name, x): print(f"{name}: [{mp.nstr(mp.mpf(x.a),10)}, {mp.nstr(mp.mpf(x.b),10)}]")
# 1. sharp threshold of smallframes (iv): 2*sqrt2*sin(3p) = 1/2
pc = iv.asin(1/(4*s2))/3; kc = 2*s2*iv.sin(pc/2)
show("psi_c", pc); show("kappa_c", kc)
assert kc.a > I(1)/12
# 2. K_iv(kappa)/kappa <= (p/kappa)*G(p), both factors increasing; evaluate at kappa_c
def G(p):
    a6, a12 = 6*p, 12*p
    psih = a6/2*iv.exp(a6)                      # >= (e^a-1-a)/a >= sup|Psi(i t)-1|, |t|<=a
    phih = a12/2 + I('0.09')*a12**2             # >= sup|Phi(i t)-1|, |t|<=a<=1 (paper's Bernoulli bound)
    assert a12.b <= 1
    return 2*3 + 2*6 + (1+iv.sqrt(I(3))/4)*(psih*(1+phih)+phih) / p * 1 + iv.sqrt(8*(18*p)**2+8*(6*p)**2)/16/p
pkc = ang(kc.b)
ratio = (pkc/I(kc.a)) * G(pkc)
show("sup_{kappa<=kappa_c} K_iv(kappa)/kappa <=", ratio)
show("limit coefficient as kappa->0 (G(0+)/sqrt2)", G(ang(I('1e-20')))/s2)
assert ratio.b < 260
show("p/kappa at kappa_c", pkc/I(kc.a))
# 3. coarse witness and link witness ratios at kappa_c (delta <= sqrt2*n*p <= c*kappa)
r = pkc/I(kc.a)
show("link witness 12*sqrt2*p/kappa", 12*s2*r); show("coarse plaquette 16*sqrt2*p/kappa", 16*s2*r)
assert (8*pkc).b < pi.a/2
# 4. explicit example: L=32, n=10, phi=2*pi*n/L^2
phi = 2*pi*10/I(1024)
kap = 2*s2*iv.sin(phi/2); d6 = 2*s2*iv.sin(3*phi); chi = 2-2*d6
show("example kappa", kap); show("example max delta(P_x,v)", d6); show("example chi", chi)
assert d6.a > 0.5 and chi.b < 1 and (48*phi).b < pi.a
show("example 1/16 sum of angles of w_x = 3*phi (angle of M)", 3*phi)
# 5. branchbound for SU(2): ||DF|| <= (7/8)*x/sin x, x = 8*sqrt2*r
rr = I(1)/32; x = 8*s2*rr; c = I(7)/8*x/iv.sin(x); show("branchbound contraction at r=1/32", c); assert c.b < 0.9
# 6. Cauchy chain of prop:explicitconstants
def chain(a, R, x0=I('1e-3')):
    e0 = a + 3*x0; ph = iv.exp(R)-1; ps = (1+ph)**3-1; d = ps/(1-e0-ps); et = (1+ph)*iv.exp(d)-1
    g0 = x0 + 3*x0/(1-e0); assert (e0+ps).b < 1 and (g0+et).b < 1
    return et/(1-g0-et)
for a in (I('1e-6'), I('6e-4')):
    print("a =", mp.nstr(a.a,3))
    for m, R in ((1, I('1e-3')), (2, I('0.05')), (3, I('0.07'))):
        S = chain(a, R); show(f"  m={m}: m! S(R)/R^m", iv.factorial(m)*S/R**m)
    e0 = a + 3*I('1e-3'); ph = iv.exp(I('0.1'))-1; ps = (1+ph)**3-1; d = ps/(1-e0-ps); et = (1+ph)*iv.exp(d)-1
    g0 = I('1e-3') + 3*I('1e-3')/(1-e0)
    show("  |z|<0.1: e0+psi", e0+ps); show("  gamma0+eta", g0+et); show("  gamma0", g0)
    assert (e0+ps).b < 0.36 and (g0+et).b < 0.91 and g0.b < 0.0041
c2 = s2*I(307); show("c_2 = sqrt2*307", c2)
# 7. radius of fibrechart (iv) with c_2 in (i),(ii): s/4+(5.4e4+9000 c_2)kappa' <= 1/195, s<=kappa'/6
Riv = (I(1)/195)/(I('5.4e4')+9000*I(435)+I(1)/24); show("R_iv (kappa+6s <=)", Riv)
Rpaper = (I(1)/195)/(4*10**5*I('2.13e4')); show("paper: kappa+s <=", Rpaper)
# 8. smallfibre (i) literal constant: (1+k)^2/(1-1.31e6 k (1+k)^2) <= 1+1.4e6 k
def ok(k):
    k = I(k); den = 1-I('1.31e6')*k*(1+k)**2
    return den.a > 0 and ((1+k)**2/den).b <= (1+I('1.4e6')*k).a
lo, hi = mp.mpf('1e-9'), mp.mpf('1e-6')
assert ok(lo) and not ok(hi)
for _ in range(80):
    m = (lo+hi)/2
    if ok(m): lo = m
    else: hi = m
print("smallfibre(i) literal 52(1+1.4e6 k) holds (paper's proof) up to k =", mp.nstr(lo,6))
# 9. flat background: sup ||D^2 c(0)[Y,Y]|| over N(Y)=1, mu=0: explicit Y attains 1/16
print("flat second derivative lower bound: 1/16 = 0.0625 (exact, prop:flatsecond); route bound kappa' <=",
      mp.nstr(((I(1)/195)/(12*I('10.6')*34*2*I('1.03')*I(1)/16)).b, 6))
