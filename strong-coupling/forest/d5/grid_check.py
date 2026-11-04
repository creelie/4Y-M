"""k_*(beta) - j_max(beta) on a fine logarithmic grid of couplings (s=5, rho=1e-8)."""
from composed_threshold import composed, kstar
from mpmath import mpf, log
thr = [composed(j, 5)[4] for j in range(0, 20)]   # beta*_j, increasing in j
worst = 0; best = 99
import math
for t in range(0, 2001):
    beta = mpf(10)**(mpf(22) + mpf(20)*t/2000)
    jm = max(j for j in range(20) if thr[j] <= beta)
    k, _ = kstar(beta)
    worst = max(worst, k - jm); best = min(best, k - jm)
print("beta in [1e22,1e42], 2001 log-spaced points: k_* - j_max ranges over", best, "..", worst)
print("16^-j beta*_j max over j<=12:", float(max(thr[j]/16**j for j in range(1,13))))
