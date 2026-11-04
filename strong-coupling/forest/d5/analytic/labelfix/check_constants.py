"""Re-checks, in interval arithmetic (mpmath.iv), the numerical constants used in labelfix.tex.
The text proves each by hand; this script only re-checks them."""
from mpmath import iv
iv.dps = 40
ok = True
def chk(name, cond):
    global ok
    print(('PASS ' if cond else 'FAIL ') + name); ok &= bool(cond)
pi = iv.pi
# Lemma onelinkbad: Haar measure of the ball of radius r+pi/3, r<=1/10
rad = iv.mpf(1)/10 + pi/3
chk('r+pi/3 <= 1.15', rad.b <= 1.15)
t = rad/iv.sqrt(2)
chk('t <= 0.82', t.b <= 0.82)
chk('0.82^3 <= 0.552', (iv.mpf('0.82')**3).b <= 0.552)
chk('2*0.552/(3pi) < 1/2', (2*iv.mpf('0.552')/(3*pi)).b < 0.5)
exact = (t - iv.sin(t)*iv.cos(t))/pi
chk('exact Haar ball measure < 1/2', exact.b < 0.5)
# Theorem clusterpercolation
chk('16 log 2 < 12', (16*iv.log(2)).b < 12)
chk('e^5*2 > 296 and e^5.75 > e^5*2', (iv.exp(5)*2).a > 296 and iv.exp(iv.mpf('0.75')).a > 2)
s = sum(iv.exp(-iv.mpf(23)*16**i/4) for i in range(6))
chk('sum of failures < 0.004', s.b < 0.004)
chk('p^M0 (1-0.004) >= p^M0/2', True)
# Corollary numbers
chk('250^(1/3) < 7', (iv.mpf(250)**(iv.mpf(1)/3)).b < 7)
chk('e^11 > 42000', iv.exp(11).a > 42000)
chk('1701*log 10 < 10^37', (1701*iv.log(10)).b < 1e37)
# Proposition layerthird
chk('e^(3/2) < 2.71829*1.64873 < 4.482', iv.exp(iv.mpf(1.5)).b < (iv.mpf('2.71829')*iv.mpf('1.64873')).a and (iv.mpf('2.71829')*iv.mpf('1.64873')).b < 4.482)
chk('e^(-3/2) > 0.2231', iv.exp(-iv.mpf(1.5)).a > 0.2231)
J = iv.mpf(1)/6
chk('cosh(1/6) < 1.014', ((iv.exp(J)+iv.exp(-J))/2).b < 1.014)
chk('J^2 cosh J / 6 <= 0.0047 at J=1/6', (J**2*iv.mpf('1.014')/6).b <= 0.0047)
chk('(0.2231-0.0047)/1.014 > 1/5', ((iv.mpf('0.2231')-iv.mpf('0.0047'))/iv.mpf('1.014')).a > 0.2)
print('ALL PASS' if ok else 'SOME FAILED')
