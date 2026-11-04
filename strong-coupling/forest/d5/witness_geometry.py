"""Brute-force check of the witness-box counts M_j=(17 2^j-4)^4, H_j (per coordinate),
and the composed threshold for the chart-calm choice rho=1e-10 (eta=9.2e-9, r=4.3e-9)."""
from composed_threshold import composed, kstar, witness_counts
from mpmath import mpf
for j in range(1, 7):
    d = 3*(2**j - 1)                      # spreading, Theorem defectscales(ii), fine units
    corners = set()
    for xp in range(-5, 6):              # scale-j witness corner, Theorem islandthreshold proof (box of side 11)
        lo, hi = 2**j*xp - d, 2**j*xp + 2**j + d     # fine coordinates reachable from the closed plaquette P
        for x in range(lo - 5, hi + 5):
            if x + 1 >= lo and x <= hi:  # [x,x+1] meets [lo,hi]
                corners.add(x)
    side = max(corners) - min(corners) + 1
    # H: a fine corner x lies in the box of y iff x - 2^{j+1} y in corners
    step = 2**(j+1)
    H1 = max(sum(1 for y in range(-30, 31) if (x - step*y) in corners) for x in range(step))
    M, H = witness_counts(j)
    print(f"j={j}: corner range [{min(corners)},{max(corners)}], side={side} (formula {17*2**j-4}), per-coordinate overlap {H1} -> H_j={H1**4} (formula {H})")
    assert side**4 == M and H1**4 == H and min(corners) == -8*2**j+2 and max(corners) == 9*2**j-3
print()
print("chart-calm labels rho=1e-10 (eta=9.2e-9<=1e-8), s=5:")
for j in range(0, 9):
    M,H,lq,rho0,b = composed(j, 5, rho=mpf('1e-10'))
    k,_ = kstar(b)
    print(f"  j={j}: beta*_j={float(b):.4e}  16^-j beta*_j={float(b/16**j):.4e}  k_*(beta*_j)={k}  k_*-j={k-j}")
