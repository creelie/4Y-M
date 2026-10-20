"""Constants of the von Mises-Fisher partition function on S^3, in ball arithmetic.

Z0(c) = int_{S^3} e^{<c,x>} dx = sum_j (|c|^2/4)^j / (j! (j+1)!),  log Z0(c) = phi(|c|^2),
phi = log g.  For c != 0 the Hessian of log Z0 has the eigenvalues 2 phi'(s) (tangential,
three times) and 2 phi'(s) + 4 s phi''(s) (radial), s = |c|^2.  The third derivative in a
unit direction v is D^3(v,v,v) = a t + g t^3 with t = <c/|c|, v>, a = 12 phi'' r,
g = 8 phi''' r^3, r = |c|; its supremum over t in [-1,1] is max(|a+g|, 2|a|/3), and by
Banach's theorem on symmetric forms this is the norm of D^3 as a trilinear form.

scan(R, n) covers 0 <= r <= R by n closed subintervals and returns rigorous bounds
lam (lower bound of every Hessian eigenvalue), top (upper bound of every eigenvalue) and
M3 (upper bound of the norm of D^3) on the ball |c| <= R.  Requires R^2 <= 5.
Usage: python3 constants.py R n   (R a decimal string)."""
from flint import arb, ctx
import sys
ctx.prec = 120
J = 40

def tail_bound():
    # For s <= 5, k <= 3 and j >= J the j-th term of the k-th derivative of g is at most
    # t_j = (5/4)^j j^3 / (j! (j+1)!), and t_{j+1}/t_j = (5/4)(1+1/j)^3/((j+1)(j+2)) < 1/2,
    # so the tail is at most 2 t_J.
    tJ = (arb(5)/4)**J * arb(J)**3 / (arb(J).fac() * arb(J+1).fac())
    ratio = (arb(5)/4) * (1 + arb(1)/J)**3 / ((J+1)*(J+2))
    assert ratio.upper() < 0.5
    return 2*tJ.upper()

TAIL = tail_bound()

def gser(s):
    """g, g', g'', g''' at the ball s (s <= 5), with the tail added as a radius."""
    g = [arb(0)]*4
    for j in range(J):
        c = arb(1)/(arb(4)**j*arb(j).fac()*arb(j+1).fac())
        for k in range(4):
            if j >= k:
                g[k] += c*arb(j).fac()/arb(j-k).fac()*s**(j-k)
    return [x + arb(0, TAIL) for x in g]

def derivs(s):
    g0, g1, g2, g3 = gser(s)
    p1 = g1/g0
    p2 = g2/g0 - p1**2
    p3 = g3/g0 - 3*g2*g1/g0**2 + 2*p1**3
    return p1, p2, p3

def scan(R, n):
    R = arb(R)
    assert (R*R).upper() <= 5
    lam = None; top = None; M3 = None
    for i in range(n):
        r = (R*i/n).union(R*(i+1)/n)
        s = r*r
        p1, p2, p3 = derivs(s)
        tang = 2*p1
        rad = 2*p1 + 4*s*p2
        a = 12*p2*r
        g = 8*p3*r**3
        m3 = max(abs(a+g).upper(), (2*abs(a)/3).upper())
        lo = min(tang.lower(), rad.lower())
        hi = max(tang.upper(), rad.upper())
        lam = lo if lam is None else min(lam, lo)
        top = hi if top is None else max(top, hi)
        M3 = m3 if M3 is None else max(M3, m3)
    return lam, top, M3

if __name__ == "__main__":
    lam, top, M3 = scan(sys.argv[1], int(sys.argv[2]))
    print("|c| <=", sys.argv[1], ": every Hessian eigenvalue in [", lam.str(10, radius=False),
          ",", top.str(10, radius=False), "],  |D^3| <=", M3.str(10, radius=False))
