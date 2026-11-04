"""Re-check of the hand arithmetic of nuprove.tex (exact algebra and mpmath interval arithmetic).
The text does not depend on this script."""
import sympy as sp
from mpmath import iv
iv.dps = 30
ok = True
def rep(name, cond):
    global ok
    ok &= bool(cond)
    print(('PASS ' if cond else 'FAIL ') + name)

# Lemma loops(a): trace identity for SU(2), with generic SU(2) matrices in quaternion form.
a0,a1,a2,a3,b0,b1,b2,b3 = sp.symbols('a0:4 b0:4', real=True)
I = sp.I
def su2(c0,c1,c2,c3):
    return sp.Matrix([[c0+I*c3, c2+I*c1],[-c2+I*c1, c0-I*c3]])
V, W = su2(a0,a1,a2,a3), su2(b0,b1,b2,b3)
nV = a0**2+a1**2+a2**2+a3**2; nW = b0**2+b1**2+b2**2+b3**2
Vi = sp.Matrix([[a0-I*a3, -a2-I*a1],[a2-I*a1, a0+I*a3]])  # adjugate = inverse when det=1
Wi = sp.Matrix([[b0-I*b3, -b2-I*b1],[b2-I*b1, b0+I*b3]])
rep('V*adj(V) = det V', sp.simplify(V*Vi - nV*sp.eye(2)) == sp.zeros(2))
x, y, z = V.trace(), W.trace(), (V*W).trace()
c = (V*W*Vi*Wi).trace()
# on det V = det W = 1 the identity reads x^2+y^2+z^2-xyz-2 = Tr(VWV^-1W^-1)
expr = sp.expand(c - (x**2+y**2+z**2-x*y*z-2))
red = sp.expand(expr.subs(a0**2, 1-a1**2-a2**2-a3**2).subs(b0**2, 1-b1**2-b2**2-b3**2))
red = sp.expand(red.subs(a0**2, 1-a1**2-a2**2-a3**2).subs(b0**2, 1-b1**2-b2**2-b3**2))
rep('Fricke identity on SU(2) x SU(2)', sp.simplify(red) == 0)

# Lemma loops(a) inequality: for S<3, |xyz| <= (S/3)^{3/2} <= S/3, so 4-q <= 4S/3.
S = sp.symbols('S', positive=True)
rep('(S/3)^(3/2) <= S/3 on (0,3)', all((sp.Rational(k,100)/3)**sp.Rational(3,2) <= sp.Rational(k,100)/3 for k in range(1,300)))

# Lemma boxenergy constants.
K = 1 + iv.log(iv.mpf(20)/3) + iv.mpf(2)/3*iv.log(3*iv.pi**3/8)
print('K =', K)
rep('K < 5', K.b < 5)
rep('log(20/3) < log 7 < 2', iv.log(iv.mpf(7)).b < 2)
rep('pi^3 < 31.3 < 32', (iv.pi**3).b < 31.3)
rep('e^2 > 7.38', (iv.e**2).a > 7.38 and iv.mpf(2.718)**2 > 7.38)
rep('2.718^5 > 148', (iv.mpf('2.718')**5).a > 148)
rep('(2/3) log 12 < 5/3', (iv.mpf(2)/3*iv.log(iv.mpf(12))).b < iv.mpf(5)/3)
rep('30 r^2 b = 9/2 at r^2=3/(20b)', sp.Rational(30)*sp.Rational(3,20) == sp.Rational(9,2))
bb = sp.symbols('b', positive=True)
r3 = (sp.Rational(3,20)/bb)**sp.Rational(3,2)
d = sp.expand_log(3*sp.log(3*sp.pi**3/(8*r3)) - 3*sp.log(3*sp.pi**3/8) - sp.Rational(9,2)*sp.log(20*bb/3), force=True)
rep('3 log(3pi^3/(8 r^3)) = 3log(3pi^3/8)+(9/2)log(20b/3)', sp.simplify(d) == 0)
rep('8 r^2 + 2 r^2 per link: 3L^3*8+3L^3*2 = 30L^3', 3*8+3*2 == 30)
rep('r^2 = 3/20 <= (pi/2)^2', (iv.pi/2)**2 > 0.15)

# Theorem centreodd: norm step.
rep('1.22^2 = 1.4884', sp.Rational(122,100)**2 == sp.Rational(14884,10000))
rep('3/2 - 3/6400 >= 1.4884', sp.Rational(3,2)-sp.Rational(3,6400) >= sp.Rational(14884,10000))
rep('(3/8)*2 = 3/4', sp.Rational(3,8)*2 == sp.Rational(3,4))
rep('6 sqrt2 / 40 < 0.22', (6*iv.sqrt(2)/40).b < 0.22)
rep('Lipschitz: sum L_e^2 <= 2sqrt2 * 3sqrt2 L = 12 L', sp.sqrt(2)*2*3*sp.sqrt(2) == 12)
rep('flow: 12 * 2 * 3 = 72', 12*2*3 == 72)
rep('energy: 2*(3 sqrt2)^2 = 36', 2*(3*sp.sqrt(2))**2 == 36)
rep('sqrt2 * 3 * 2 = 6 sqrt2 (3L^2 plaquettes, |U_p-I| norm 2 eps^1/2)', sp.sqrt(2)*3*2 == 6*sp.sqrt(2))
print('ALL PASS' if ok else 'SOME FAIL')
