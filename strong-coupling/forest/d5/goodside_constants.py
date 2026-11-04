"""Constants of Lemma 2 of NOTE.md (extended good side r <= 1/32)."""
from mpmath import mp, mpf, sin, pi, sqrt
mp.dps = 30
r = mpf(1)/32
print("16 r =", 16*r, " (<= 1/2: cut-off equals one at every branch)")
theta = 16*r/sqrt(2)                      # largest eigenangle of Y = log(e^{-X}A_x), ||Y|| <= 16 r
L = theta/sin(theta)
print("theta = 16r/sqrt2 =", theta, "  theta/sin theta =", L)
k = mpf(7)/8*L
print("contraction constant (7/8) theta/sin theta =", k, " < 0.8936:", k < mpf('0.8936'))
print("singular values of the path-link derivative in [1-k, 1+k] = [", 1-k, ",", 1+k, "]")
print("|det| in [", (1-k)**3, ",", (1+k)**3, "]")
print("response to one non-path private link: (2/16) L /(1-k) =", (2*L/16)/(1-k))
print("response to one inside link u_a or w_a: (1/16) L /(1-k) =", (L/16)/(1-k))
print("pi/3 =", pi/3, " (conversion d <= (pi/3) delta for delta <= 1)")
