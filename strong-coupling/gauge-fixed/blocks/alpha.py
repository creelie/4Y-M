"""Column sum of the plaquette-block Dobrushin matrix for SU(2), in ball arithmetic.

beta = beta_W/4.  With lam, top, M3 the constants of constants.py on |c| <= 24 beta:
  w   = 8 beta^2 (top - lam)                 oscillation of the linearisation remainder
  s1  = 1/4 + (e^w - 1) e^w sqrt(top)        sensitivity of the two-site conditional law to b1
  C   = 20 beta (4 beta top s1 + 16 beta^2 M3)   coupling of the two sites
  K   = (1 + beta) s1 / (1 - C)              sensitivity of the block law to one field
  alpha = 40 beta K                          column sum (ten outside partners, 4 beta each)
Usage: python3 alpha.py 1/3 0.34 0.345 ...  (decimals or fractions, read exactly)."""
from flint import arb, ctx
from constants import scan
import sys
ctx.prec = 120

def parse(bW):
    if "/" in bW:
        p, q = bW.split("/")
        return arb(p)/arb(q)
    return arb(bW)

def alpha(bW, n=20000):
    beta = parse(bW)/4
    lam, top, M3 = scan((24*beta).upper(), n)
    w = 8*beta**2*(top - lam)
    s1 = arb(1)/4 + (w.exp() - 1)*w.exp()*top.sqrt()
    C = 20*beta*(4*beta*top*s1 + 16*beta**2*M3)
    assert C.upper() < 1
    K = (1 + beta)*s1/(1 - C)
    return dict(lam=lam, top=top, M3=M3, w=w, s1=s1, C=C, K=K, alpha=40*beta*K)

if __name__ == "__main__":
    for bW in sys.argv[1:]:
        d = alpha(bW)
        print("beta_W =", bW, " ".join("%s=%s" % (k, v.str(8, radius=False)) for k, v in d.items()),
              " alpha <=", d["alpha"].upper().str(8, radius=False))
