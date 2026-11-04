"""Certified constants of Proposition prop:energybeta (ssec:slabenergy); A, B, D here are K_1, K_2, K_3 of the paper.

Interval arithmetic with mpmath.iv at 60 decimal digits.  With lambda = 1/5 and
c = 5/8 the per-plaquette energy bound of the proposition reads
    eps <= ( A + B/L^4 + (15/16)(1/L + 1/(3L^4)) log beta + D/beta ) / beta
for beta >= 1 and L >= 2, where
    A = (5/4) K,
    K = (1/4)log(2/pi) - (3/4)log(lam) - (1/2)log(2/(3pi)) - (3/4)log(c) + (6/5)c,
    B = (5/24)( -log(2/(3pi)) - (3/2)log(c) ),
    D = (5/4)(1 + 1/48) * (c/6) / (1 - c/6).
The script also checks the sign conditions used in the proof and scans the
choice of lambda and c on a grid (floating point, only to report that the
chosen pair is near the optimum; no statement rests on the scan).
"""
from mpmath import iv, mp, mpf, log, pi
iv.dps = 60
mp.dps = 30

lam = iv.mpf(1) / 5
c = iv.mpf(5) / 8
PI = iv.pi
L = iv.log

def show(name, x):
    print(f"{name}: {iv.nstr(x, 15)}")

half_log_2pi = L(2 / PI) / 2
ball = L(2 / (3 * PI))
K = half_log_2pi / 2 - 3 * L(lam) / 4 - ball / 2 - 3 * L(c) / 4 + 6 * c / 5
A = 5 * K / 4
B = 5 * (-ball - 3 * L(c) / 2) / 24
y = c / 6                      # delta^2/6 at beta = 1, the largest value
D = 5 * (1 + iv.mpf(1) / 48) / 4 * y / (1 - y)
dropped = half_log_2pi - 3 * L(lam) / 2   # coefficient of the dropped 1/L term

show("K", K)
show("A", A)
show("B", B)
show("D", D)
show("(1/2)log(2/pi) - (3/2)log(lambda)  [must be >= 0]", dropped)
show("delta^2 at beta=1 = c  [must be <= pi^2/4]", c)
show("pi^2/4", PI**2 / 4)
assert dropped.a > 0
assert c.b < (PI**2 / 4).a
assert A.b < mpf("3.7148") and B.b < mpf("0.4699") and D.b < mpf("0.149")
print("certified: A < 3.7148, B < 0.4699, D < 0.149")

# Gaussian heuristic comparison constant (not used in any statement)
def Kf(l, cc):
    return (log(2 / pi) / 4 - 0.75 * log(l) - 0.5 * log(2 / (3 * pi))
            - 0.75 * log(cc) + 1.2 * cc) / (1 - l)
best = min((Kf(l / 1000, cc / 1000), l / 1000, cc / 1000)
           for l in range(20, 900, 5) for cc in range(100, 2000, 5))
print("grid minimum of A over (lambda, c), floating point:", best)
