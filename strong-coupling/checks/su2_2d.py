"""2D SU(2) Wilson theory: explicit finite-box constants.

f = exp(beta cos theta)/c_0 has character coefficients a_j = I_{2j+1}(beta)/I_1(beta);
f^{*m} = sum_j d_j a_j^m chi_j, so eps_m := sup|f^{*m}-1| = sum_{j>=1/2} d_j^2 a_j^m
(attained at the identity).  eta_k := 2 eps_k/(1-eps_k^2) bounds the total-variation
influence of the boundary holonomy of a disk of k plaquettes on everything outside it.

DS box: V = interior links of an L x L square (2L^2-2L links), dV = 4L boundary links.
Coupling of Prop. 2D: k(y) <= 2 s(s-1) + eta_{s^2} * |V|.  DS condition: 4L k < |V|.
Run: python3 -I su2_2d.py
"""
import math
from scipy.special import ive

def eps(beta, m):
    s = 0.0
    for j2 in range(1, 400):           # j2 = 2j, d_j = j2+1
        a = ive(j2 + 1, beta) / ive(1, beta)
        t = (j2 + 1) ** 2 * a ** m
        s += t
        if t < 1e-300 or (j2 > 5 and t < 1e-18 * s):
            break
    return s

def eta(beta, k):
    e = eps(beta, k)
    return float('inf') if e >= 1 else 2 * e / (1 - e * e)

def ds_box(beta, Lmax=100000):
    """smallest L for which some s gives 4L*(2s(s-1)+eta_{s^2}|V|) < |V|."""
    L = 2
    while L <= Lmax:
        V = 2 * L * L - 2 * L
        for s in range(1, L + 1):
            k = 2 * s * (s - 1) + eta(beta, s * s) * V
            if 4 * L * k < V:
                return L, s, k
            if 2 * s * (s - 1) * 4 * L >= V:
                break
        L = L + 1 if L < 200 else int(L * 1.05) + 1
    return None

if __name__ == "__main__":
    print("beta   a_1/2=I2/I1   sigma=-log a   L_DS  s   k(y)bound   Conj-G delta(L)<=eta_{L^2-1} at L_DS")
    for beta in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0]:
        a = ive(2, beta) / ive(1, beta)
        res = ds_box(beta)
        L, s, k = res
        print(f"{beta:5.1f}  {a:.6f}     {-math.log(a):.5f}      {L:5d} {s:3d}  {k:9.3f}   {max(eta(beta, L*L-1),1e-300):.1e}")
