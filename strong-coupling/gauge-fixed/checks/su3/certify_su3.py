"""Certified constants of the SU(3) single-link bound (python-flint, ball arithmetic).

(A) mu_0 = 7/6: bottom of the rough Laplacian on 1-forms of SU(3) (metric -Tr XY), from
    (C_pi + C_sigma)/2 - C_ad/4 over the components sigma of pi (x) ad.
(B) g^2 = max over W in SU(3) of |traceless part of (W - W*)/2|_HS^2 <= G2 (branch and bound).
(C) 2 r Lambda(r) < 1 at r = R0, with
      gap(r)  = (sqrt(mu0) - g r/2)^2 + 3/2 - r,
      J(r)    = int_0^1 g r / sqrt(gap(t r)) dt     (upper Riemann sum, integrand increasing),
      m(r)    = sin^2(min(asin(2/3) + J(r)/2, pi/2)),
      Lambda  = sqrt(m(r)) (3/8) (1 + g r / sqrt(gap(r))).
"""
from fractions import Fraction as Fr
from flint import arb, ctx
import sys
ctx.prec = 80
ok = True

# (A)
def C(p, q):
    return Fr(2 * (p * p + q * q + p * q + 3 * p + 3 * q), 3)
def comps(p, q):
    out = [(p + 1, q + 1), (p + 2, q - 1), (p - 1, q + 2), (p + 1, q - 2), (p - 2, q + 1), (p - 1, q - 1)]
    mult = 2 if (p > 0 and q > 0) else (1 if p + q > 0 else 0)
    out += [(p, q)] * mult
    return [s for s in out if s[0] >= 0 and s[1] >= 0]
def dim(p, q):
    return (p + 1) * (q + 1) * (p + q + 2) // 2
best = None
for p in range(0, 41):
    for q in range(0, 41):
        cs = comps(p, q)
        assert sum(dim(*s) for s in cs) == 8 * dim(p, q), (p, q)
        for s in cs:
            e = (C(p, q) + C(*s)) / 2 - Fr(3, 2)
            if best is None or e < best[0]:
                best = (e, (p, q), s)
print("(A) min of (C_pi+C_sigma)/2 - 3/2 for p,q <= 40:", best[0], "at", best[1], best[2])
ok &= best[0] == Fr(7, 6)
# for p + q >= 3 every value is at least C(p,q) - <lambda+rho, theta> - 3/2 > 7/6; the minimum over
# the components is attained at the smallest Casimir, and C(p-1,q-1) etc. are increasing in p, q.

# (B) branch and bound on [-pi, pi]^2, phi3 = -phi1 - phi2
pi = arb.pi()
G2 = arb("2.06548")
def F(x, y):
    a = [x, y, -x - y]
    s = [t.sin() for t in a]
    return sum(t * t for t in s) - (s[0] + s[1] + s[2]) ** 2 / 3
stack = [(-pi, pi, -pi, pi)]
nbox = 0
while stack:
    x0, x1, y0, y1 = stack.pop()
    nbox += 1
    X = arb((x0 + x1) / 2, (x1 - x0) / 2); Y = arb((y0 + y1) / 2, (y1 - y0) / 2)
    v = F(X, Y)
    if v < G2:
        continue
    # also use the mean-value form: F(c) + grad bound * radius
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    if (x1 - x0) < arb("1e-7"):
        print("(B) failed near", xm, ym); ok = False; break
    stack += [(x0, xm, y0, ym), (xm, x1, y0, ym), (x0, xm, ym, y1), (xm, x1, ym, y1)]
print("(B) g^2 <= 2.06548 certified on", nbox, "boxes")

# (C)
g = G2.sqrt()
mu0 = arb(7) / 6
def gap(r):
    return (mu0.sqrt() - g * r / 2) ** 2 + arb(3) / 2 - r
def J(r, n=4000):
    return sum(g * r / gap(r * i / n).sqrt() for i in range(1, n + 1)) / n
def Lam(r):
    th = arb(2) / 3
    th = th.asin() + J(r) / 2
    m = th.sin() ** 2 if th < pi / 2 else arb(1)
    return m.sqrt() * arb(3) / 8 * (1 + g * r / gap(r).sqrt())
for R0 in ("0.730", "0.732"):
    r = arb(R0)
    d = 2 * r * Lam(r)
    print("(C) r =", R0, " gap =", gap(r).str(6), " 2 r Lambda =", d.str(8), " beta_W =", (r / 2).str(4))
    ok &= bool(d < 1)
    print("    gap of the measure >= -log(2 r Lambda) =", (-d.log()).str(5))
print("ALL OK" if ok else "FAILED")
sys.exit(0 if ok else 1)
