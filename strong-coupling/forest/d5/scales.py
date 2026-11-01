"""Item 5 scoping (D4_D8.md): the effective coupling at the coarsest scale of (D).

Heuristic, one-loop running only. SU(2): beta_W(mu) = 2N/g^2 = 8 b0 log(mu/Lambda), b0 = 11N/(48 pi^2).
(D) at physical size ell forces a gap >= log(1/alpha)/(2 ell) (paper, 'continuum step').
If the true gap is m, then ell >= log(1/alpha)/(2m).  The cube has n0 blocks per side,
so the coarsest blocks have size ell/n0 and effective coupling beta_top = 8 b0 log(n0/(ell Lambda)).
"""
from math import log, pi
N = 2
b0 = 11*N/(48*pi**2)
c = 8*b0
print(f"beta_W = {c:.4f} * log(mu/Lambda)")
alpha = 0.9
for m_over_L in (1.0, 10.0, 100.0):           # gap in units of Lambda (lattice vs MSbar differ by ~20)
    ell_L = log(1/alpha)/(2*m_over_L)          # minimal ell*Lambda
    for n0, name in ((5.67e14, "chart C_G"), (2e15/alpha, "Gaussian constrained"), (None, "Gaussian all scales")):
        lg = log(n0) if n0 is not None else 827*log(10)-log(alpha)
        bt = c*(lg - log(ell_L))
        print(f"m/Lambda={m_over_L:6.1f}  n0={name:22s} beta_top <= {bt:9.1f}")
print()
for thr, name in ((3.8e38, "chartdecay/fibrechart"), (3.8e97, "off-chart window"), (8.3e163, "islands")):
    print(f"threshold {name:22s} beta >= {thr:.1e}: needs log(n0/(ell Lambda)) >= {thr/c:.2e}")
