"""Column sum of the plaquette-block Dobrushin matrix for SU(3), in ball arithmetic.

N beta = beta_W/3, r = 2 beta_W.  Every middle-link field c of a block lies in r conv SU(3).  On that set the
certified covers give
  2 r Lambda <= B      (cover2.py of ../../checks/su3/deg2 with --thr: one-link sensitivity, cover_log.txt),
  tau, lam, M          (consts.py: Hessian eigenvalues of log Z in [lam, tau], third cumulant <= M, consts_log.txt),
and the spectral gap of nu_c is at least gamma = 283/128 - r/2 (Lemma 19).  Then
  w     = (3/2) (N beta)^2 (tau - lam)                oscillation of the Taylor remainder (|x|^2 = 3)
  s     = Lambda + (e^w - 1) e^w sqrt(tau/gamma)      sensitivity of the two-site conditional law to b1
  C     = 5 (N beta)^2 (tau s + sqrt(3) N beta M e^w / gamma)   coupling of the two outer links
  K     = (1 + N beta Lambda) s / (1 - C)             sensitivity of the block law to one outside field
  alpha = 10 N beta K                                 column sum (ten outside partners, coupling N beta each)
All quantities increase with beta_W, so constants certified at r cover every beta_W <= r/2.
Usage: python3 alpha3.py COVERDIR CONSTSDIR beta_W ...   (cover_log.txt and consts_log.txt of the same r; beta_W <= r/2),
  e.g. python3 alpha3.py ../../checks/su3/deg2/r3_2 r3_2 3/4."""
import os, re, sys
from fractions import Fraction
from flint import arb, ctx
ctx.prec = 120

def read(dc, dk):
    cv = open(os.path.join(dc, 'cover_log.txt')).read()
    m = re.search(r"r = (\d+)/(\d+)", cv); rq = Fraction(int(m.group(1)), int(m.group(2))); r = parse(str(rq))
    B = re.findall(r'COVER COMPLETE: \d+ boxes, largest bound ([0-9.]+)', cv)[-1]
    cs = open(os.path.join(dk, 'consts_log.txt')).read()
    t, l, M = re.findall(r'tau <= ([0-9.]+), lambda >= ([0-9.]+), M <= ([0-9.]+)', cs)[-1]
    u = arb('1e-6')                            # the logs print six decimals rounded to nearest
    return rq, r, arb(B) + u, arb(t) + u, arb(l) - u, arb(M) + u

def parse(x):
    q = Fraction(x)                           # decimals and fractions are read exactly
    return arb(q.numerator) / arb(q.denominator)

def alpha(dc, dk, bW):
    r0q, r0, B, tau, lam, M = read(dc, dk)
    assert r0q == Fraction(re.search(r'r = (\d+/\d+)', open(os.path.join(dk, 'consts_log.txt')).read()).group(1)), 'covers at different r'
    assert 2 * Fraction(bW) <= r0q, "beta_W beyond the certified range"
    bW = parse(bW)
    Lam = B / (2 * r0)                        # sup over r0 conv SU(3), which contains 2 beta_W conv SU(3)
    Nb = bW / 3
    gam = arb(283) / 128 - r0 / 2             # the gap at r0 is the smallest on the range
    w = arb(3) / 2 * Nb ** 2 * (tau - lam)
    e = w.exp()
    s = Lam + (e - 1) * e * (tau / gam).sqrt()
    C = 5 * Nb ** 2 * (tau * s + arb(3).sqrt() * Nb * M * e / gam)
    assert C.upper() < 1
    K = (1 + Nb * Lam) * s / (1 - C)
    return dict(Lambda=Lam, tau=tau, lam=lam, M=M, gamma=gam, w=w, s=s, C=C, K=K, alpha=10 * Nb * K)

if __name__ == '__main__':
    for bW in sys.argv[3:]:
        a = alpha(sys.argv[1], sys.argv[2], bW)
        print('beta_W =', bW, ' '.join('%s=%s' % (k, v.str(8, radius=False)) for k, v in a.items()),
              ' alpha <=', a['alpha'].upper().str(8, radius=False))
