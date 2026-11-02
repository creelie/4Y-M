#!/usr/bin/env python3
"""Sharpened plaquette-block constant alpha(beta_W) = 40 |beta| K for SU(2), in ball arithmetic.

g(c) = log Z0(c) = phi(|c|^2) with Z0(c) = sum_j (|c|^2/4)^j / (j! (j+1)!), the normalisation of
the von Mises-Fisher law on S^3 = SU(2).  On the ball |c| <= R the program encloses

    lambda_*  <= every eigenvalue of Hess g,      tau >= every eigenvalue of Hess g,
    M         >= the operator norm of D^3 g,      A   >= |grad g(c)| = |c| 2 phi'(|c|^2) on |c| <= R20,

by python-flint at 128 bits on a subdivision of the radius, with the series of phi and its first
three derivatives summed to J terms and the tail bounded by twice the first omitted term (the
terms decrease by a factor below one half on the whole ball).  The eigenvalues of Hess g at c are
2 phi' (tangential, multiplicity three) and 2 phi' + 4 |c|^2 phi'' (radial); D^3 g has norm at most
max(|A + G|, 2|A|/3) with A = 12 r phi'' and G = 8 r^3 phi'''.

The chain of the paper (Lemma onesite (ii) with the exponents interpolated, and the Poincare
constant 1/3 of the von Mises-Fisher law):

    w = 8 beta^2 (tau - lambda_*),      s = 1/4 + (e^w - 1) e^w sqrt(tau/3),
    C = 4 |beta| s A + (320/3) |beta|^3 M e^w,     K = (1 + |beta|) s / (1 - C),     alpha = 40 |beta| K,

with beta = beta_W / 4.  The old chain (s with sqrt(tau), C = 20|beta|(4|beta| tau s + 16 beta^2 M))
is available as --mode paper.  With --lean the four constants are rounded outward to 10^-12 and
written as Lean definitions for lean/BlockCheck.lean.

    python3 sharp.py --n 200000 --lean lean/BlockData.lean 1/3 0.345 0.35 0.351
"""
import argparse
import time
from fractions import Fraction

from flint import arb, ctx

ctx.prec = 128
J = 48


def coeffs(m):
    """S^(m)(t) = sum_{j >= m} t^(j-m) / (4^j (j-m)! (j+1)!) for S(t) = sum_j t^j / (4^j j! (j+1)!)."""
    return [arb(1) / (arb(4) ** j * arb.fac_ui(j - m) * arb.fac_ui(j + 1)) for j in range(m, J)]


COEF = [coeffs(m) for m in range(4)]


def tail(m, T):
    """bound for the omitted terms j >= J of the m-th derivative series on 0 <= t <= T."""
    term = arb(T) ** (J - m) / (arb(4) ** J * arb.fac_ui(J - m) * arb.fac_ui(J + 1))
    ratio = arb(T) / (4 * (J - m + 1) * (J + 2))
    assert ratio.upper() < 0.5
    return 2 * term


def horner(cs, t):
    s = arb(0)
    for c in reversed(cs):
        s = s * t + c
    return s


def phis(t, T):
    S = [horner(COEF[m], t) + arb(0, tail(m, T).upper()) for m in range(4)]
    p1 = S[1] / S[0]
    p2 = S[2] / S[0] - p1 * p1
    p3 = S[3] / S[0] - 3 * S[2] * S[1] / (S[0] * S[0]) + 2 * p1 ** 3
    return p1, p2, p3


def ball(lo, hi):
    lo, hi = arb(lo), arb(hi)
    return arb((lo + hi) / 2, ((hi - lo) / 2).upper())


def constants(R, R20, n):
    R = arb(R)
    T = (R * R).upper() * 1.0001
    lam = tau = M = A = None
    for k in range(n):
        rlo = R * k / n if k else arb(0)
        rhi = R * (k + 1) / n
        r = ball(rlo.lower(), rhi.upper())
        t = ball(arb(rlo.lower()) ** 2, arb(rhi.upper()) ** 2)
        p1, p2, p3 = phis(t, T)
        assert p2.upper() <= 0, "phi'' must be nonpositive"
        tan = 2 * p1
        rad = 2 * p1 + 4 * t * p2
        lo = min(tan.lower(), rad.lower())
        hi = max(tan.upper(), rad.upper())
        a3 = 12 * r * p2
        g3 = 8 * r ** 3 * p3
        d3 = max(abs(a3 + g3).upper(), (2 * abs(a3) / 3).upper())
        lam = lo if lam is None else min(lam, lo)
        tau = hi if tau is None else max(tau, hi)
        M = d3 if M is None else max(M, d3)
        if r.lower() <= arb(R20).upper():
            a = (r * tan).upper()
            A = a if A is None else max(A, a)
    return lam, tau, M, A


def chain(bW, consts, mode):
    lam, tau, M, A = consts
    b = abs(arb(bW)) / 4
    w = 8 * b * b * (tau - lam)
    ew = w.exp()
    if mode == "paper":
        s = arb(1) / 4 + (ew - 1) * ew * tau.sqrt()
        C = 20 * b * (4 * b * tau * s + 16 * b * b * M)
    else:
        s = arb(1) / 4 + (ew - 1) * ew * (tau / 3).sqrt()
        C = 4 * b * s * A + arb(320) / 3 * b ** 3 * M * ew
    assert C.upper() < 1
    K = (1 + b) * s / (1 - C)
    return w.upper(), s.upper(), C.upper(), K.upper(), (40 * b * K).upper()


def show(x):
    return x.str(10, radius=False)


def parse(bW):
    f = Fraction(bW)
    return arb(f.numerator) / arb(f.denominator)


def rnd(x, up, D=10 ** 12):
    y = x * D
    return int((y.ceil() if up else y.floor()).unique_fmpz())


def write_lean(path, lam, tau, M, A):
    vals = [("lamLo", rnd(lam, False), "lower bound for every eigenvalue of Hess g on |c| <= 24 beta"),
            ("tauHi", rnd(tau, True), "upper bound for every eigenvalue of Hess g on |c| <= 24 beta"),
            ("mHi", rnd(M, True), "upper bound for the norm of D^3 g on |c| <= 24 beta"),
            ("aHi", rnd(A, True), "upper bound for |grad g(c)| on |c| <= 20 beta")]
    lines = ["/- Written by sharp.py: Arb enclosures rounded outward, denominator 10^12. -/", "",
             "def den : Nat := 1000000000000"]
    for name, v, doc in vals:
        lines += [f"/-- {doc} -/", f"def {name} : Nat := {v}"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return [v for _, v, _ in vals]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["paper", "sharp"], default="sharp")
    ap.add_argument("--n", type=int, default=200000)
    ap.add_argument("--ball", default=None, help="radius of the ball (default 24 beta of the last coupling)")
    ap.add_argument("--lean", default=None)
    ap.add_argument("betaW", nargs="+")
    a = ap.parse_args()
    bmax = parse(a.betaW[-1]) / 4
    R = parse(a.ball) if a.ball else 24 * bmax
    R20 = 20 * bmax
    t0 = time.time()
    lam, tau, M, A = constants(R, R20, a.n)
    print(f"# mode={a.mode} ball r <= {show(R.upper())} subintervals={a.n} terms={J} prec={ctx.prec} [{time.time()-t0:.1f}s]")
    print(f"# lambda_* >= {show(lam)}  tau <= {show(tau)}  M <= {show(M)}  A(|c| <= {show(R20.upper())}) <= {show(A)}")
    consts = (lam, tau, M, A)
    if a.lean:
        rounded = write_lean(a.lean, lam, tau, M, A)
        D = 10 ** 12
        consts = tuple(arb(v) / D for v in rounded)
        print(f"# wrote {a.lean}: {rounded}")
    for bW in a.betaW:
        w, s, C, K, al = chain(parse(bW), consts, a.mode)
        print(f"beta_W={bW}: w<={show(w)} s<={show(s)} C<={show(C)} K<={show(K)} alpha<={show(al)} "
              + ("alpha < 1" if al < 1 else "FAILS"))


if __name__ == "__main__":
    main()
