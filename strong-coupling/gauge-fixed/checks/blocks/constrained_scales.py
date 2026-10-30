"""Arithmetic of Theorem "constrainedscales": the constrained cube condition at
every scale for the effective forms S_k, with constants independent of k.
The inputs are C = 1e26 and eta = 1e-36 of Theorem "effectivelocal", M_0 = 105
and the weights of Theorem "constrainedcube". Every inequality of the proof is
evaluated at 80 significant digits, and the counting bounds it uses are
checked by enumeration on small cases. Run: python3 constrained_scales.py
(a few seconds). Needs mpmath."""
from mpmath import mp, mpf, exp, coth, factorial, log10
from itertools import product
mp.dps = 80
C, eta, M0 = mpf(10)**26, mpf(10)**-36, 105
h = 256*24; s = 4*h; lam = mpf(10)**-215
ok = [0, 0]
def check(name, cond):
    ok[0 if cond else 1] += 1
    print(('PASS ' if cond else 'FAIL ') + name)
def tail(p, x):  # bound for sum_{m>=1} m^p e^{-x m}
    return exp(x)*factorial(p)/x**(p+1)
def main():
    # counting facts, by enumeration
    for r in range(1, 9):
        nr = sum(1 for z in product(range(-r, r+1), repeat=4) if sum(map(abs, z)) == r)
        check('points with |z|_1 = %d: %d <= 8 r^3' % (r, nr), nr <= 8*r**3)
    for p in range(5):
        for x in (mpf('0.3'), mpf(1), mpf(2)):
            tot = sum(mpf(m)**p*exp(-x*m) for m in range(1, 4000))
            check('sum m^%d e^{-%s m} <= e^x p!/x^(p+1)' % (p, x), tot <= tail(p, x))
    for n in (2, 4, 10, 30):
        for m in range(1, 8):
            check('shell n=%d m=%d' % (n, m),
                  4*((n+2*m)**4-(n+2*m-2)**4) <= 128*(n**3+8*m**3))
    for m in range(1, 8):
        sites = (2*m+1)**4-(2*m-1)**4
        check('sites at l-infinity distance %d: %d = 64m^3+16m' % (m, sites), sites == 64*m**3+16*m)
    # Lemma saddle: the bounds for (h,s) do not exceed those for (24,96)
    check('h/s unchanged', mpf(h)/s == mpf(24)/96)
    # row sums of Lemma combeslong
    rowH = 24*C*lam*exp(lam)*(1+16*24*exp(eta/2)*(2/eta)**5)
    rows = max(rowH + s*(exp(lam)-1)/8, 2*s*(exp(lam)-1))
    check('lambda_1 <= eta/2', lam <= eta/2)
    check('row sums %.3g < 3e-4 < 1/(2 M_0)' % float(rows), rows < mpf(3)/10**4 < mpf(1)/(2*M0))
    fw = 24*C*exp(lam)*coth((eta-lam)/2)**4
    check('|e^{lam d} f| < 4e172 (10^%.3f)' % float(log10(fw)), fw < 4*mpf(10)**172)
    R2 = (2*M0*(4*mpf(10)**172+exp(lam)*s/16))**2
    check('weighted norm squared < 1e350 (10^%.3f)' % float(log10(R2)), R2 < mpf(10)**350)
    Z = lambda x: 4*exp(x)*(1+384*exp(x)/x**4+16*exp(x)/x**2)
    check('Z(lambda_1) < 1e864 (10^%.3f)' % float(log10(Z(lam))), Z(lam) < mpf(10)**864)
    check('l1 bound: 1e350 * Z < 1e1214', mpf(10)**350*Z(lam) < mpf(10)**1214)
    Lw = 2*M0*s*Z(2*lam)**0.5
    check('L_w < 1e438 (10^%.3f)' % float(log10(Lw)), Lw < mpf(10)**438)
    q = lam/2
    for alpha in (mpf('0.99'), mpf('0.5'), mpf('1e-3'), mpf('1e-12')):
        n = mpf(10)**827/alpha
        K = mpf(10)**607*128*exp(2*q)*(n**3/q+48/q**4)
        check('K_c <= alpha N_n at n = 1e827/alpha, alpha=%s (ratio %.3g)' % (alpha, float(K/(alpha*(n/2)**4))),
              K <= alpha*(n/2)**4)
    print('%d checks passed, %d failed' % tuple(ok))
if __name__ == '__main__':
    main()
