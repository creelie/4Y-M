"""Repeats, in 50-digit interval arithmetic, the hand arithmetic of the text:
the count |N_1(S)| = 243, the eta table of thm:sfgaugeobst, the constants of
thm:sfprop, and the margin of thm:sfperp.  The text derives each by hand."""
import itertools
from mpmath import iv
iv.dps = 50
# |N_1(S)| on a coarse torus of side 9, by levels in e_0
N = 9
y = (4, 4, 4, 4)
def sh(b, mu, k):
    b = list(b); b[mu] = (b[mu] + k) % N; return tuple(b)
z = sh(y, 0, 1)
S = {y, z, sh(y, 0, -1)} | {sh(y, m, -1) for m in (1, 2, 3)} | {sh(z, m, -1) for m in (1, 2, 3)}
def dist(a, b):
    return max(min((a[i] - b[i]) % N, (b[i] - a[i]) % N) for i in range(4))
N1 = [b for b in itertools.product(range(N), repeat=4) if min(dist(b, s) for s in S) <= 1]
lev = {}
for b in N1:
    l = (b[0] - y[0] + N // 2) % N - N // 2
    lev[l] = lev.get(l, 0) + 1
print("|S| =", len(S), " |N_1(S)| =", len(N1), " by e_0-level:", dict(sorted(lev.items())))
# eta table
def eta2(phi):
    s, c = iv.sin(phi), iv.cos(phi)
    cp = c * c; sp = iv.sqrt(1 - cp * cp); psi = iv.atan2(sp, cp); k = psi / sp
    return 4 * (k * s * c - phi) ** 2 + 2 * (k * s * s) ** 2
print("phi=pi/2: |eta|^2 =", eta2(iv.pi / 2), " 3pi^2/2 =", 3 * iv.pi ** 2 / 2, " 5.4/pi^2 =", iv.mpf('5.4') / iv.pi ** 2)
for p in ['0.5', '0.1', '0.05']:
    ph = iv.mpf(p); e = eta2(ph)
    print("phi", p, "|eta|^2", e, " 2 sin phi", 2 * iv.sin(ph), " 8.1/|eta|^2", iv.mpf('8.1') / e,
          " 4.05/sin^4", iv.mpf('4.05') / iv.sin(ph) ** 4)
print("4.05/0.04997^4 =", iv.mpf('4.05') / iv.mpf('0.04997') ** 4)
# thm:sfprop
r2 = iv.sqrt(2)
p1 = 7 / (2 * r2) + iv.mpf(7) / 8
Sb = 168 + 42 * r2; bb = 168 + 336 * r2
print("p1", p1, " 48 sqrt2 p1", 48 * r2 * p1, " Sbar", Sb, " bbar = 192 p1", 192 * p1, bb)
c0 = iv.mpf(537) / 5000
def ca(a):
    return c0 - ((iv.exp(3 * a) + iv.exp(-3 * a)) / 2 - 1) * Sb
for a in ['0.005', '0.01']:
    print("alpha", a, " c_alpha", ca(iv.mpf(a)))
a = iv.mpf('0.005'); c = ca(a)
print("hand: (0.015^2/2)*1.0002*227.40 =", iv.mpf('0.015') ** 2 / 2 * iv.mpf('1.0002') * iv.mpf('227.40'),
      " (0.03^2/2)*1.0005*227.40 =", iv.mpf('0.03') ** 2 / 2 * iv.mpf('1.0005') * iv.mpf('227.40'))
print("1/0.0818 =", 1 / iv.mpf('0.0818'), " e^{2a} p1^2/0.0818 =", iv.exp(2 * a) * p1 ** 2 / iv.mpf('0.0818'),
      " e^{4a} bbar^2/0.0818 =", iv.exp(4 * a) * bb ** 2 / iv.mpf('0.0818'))
# thm:sfperp margin
k = iv.mpf('1e-8')
D = (21 + 8 * iv.pi)
print("21+8pi =", D, " 46.14^2 =", iv.mpf('46.14') ** 2, " 3*46.14^2 =", 3 * iv.mpf('46.14') ** 2,
      " 3.54*46.14^3 =", iv.mpf('3.54') * iv.mpf('46.14') ** 3)
lhs = 6 * iv.sin(11 * k / (15 * r2)) / k
print("6 sin(11k/(15 sqrt2))/k at k=1e-8:", lhs, " 11 sqrt2/5 =", 11 * r2 / 5)
eps4 = 3 * (D * k) ** 2 + iv.mpf('3.54') * (D * k) ** 3
marg = lhs - 2 * r2 - 12 * eps4 / k - 540 * iv.mpf('1e-4') / r2
print("12 eps4/k", 12 * eps4 / k, " 540 rho/(sqrt2 k)", 540 * iv.mpf('1e-4') / r2, " margin/kappa", marg)
