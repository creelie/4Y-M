"""First-order ceiling of a block shape (estimate, not a bound).

For an active link j, every plaquette p containing j and meeting other blocks moves the fields of its other
active links by 4 beta per unit of j.  A block B met in p responds, at first order, like a one-link vMF law
(sensitivity 1/4) amplified by its internal plaquettes: a link i in m_i plaquettes lying inside B passes
the response to the rest of each such plaquette with weight beta (the z = y2 y3 route of the plaquette block,
which gives the known (1 + beta)/4 for plaquette blocks).  Two links of p in the same block count once.
    alpha_j(beta_W) ~ (beta_W/4) sum_{p ni j} sum_{B != B(j) met in p} (1 + beta * max_{i in p cap B} m_i),
beta = beta_W/4.  Prints the beta_W where max_j alpha_j = 1.  For plaquette blocks this is 10(beta_W/4)(1+beta)
= 1, i.e. 0.366, the known ceiling.  Ignores chains of length >= 2 inside a block, which only raise K.

    python3 ceiling.py [P]
"""
import sys
import numpy as np
from geometry import Lattice


def ceiling(L, lab):
    m = L.internal_counts(lab)
    # coef[j] = (n_j, s_j) with alpha_j = (bW/4) (n_j + (bW/4) s_j)
    n = np.zeros(len(L.active)); s = np.zeros(len(L.active))
    for _, ls in L.plaqs:
        for j in ls:
            blocks = {}
            for e in ls:
                if lab[e] != lab[j]:
                    blocks[lab[e]] = max(blocks.get(lab[e], 0), m[e])
            n[j] += len(blocks); s[j] += sum(blocks.values())
    lo, hi = 0.0, 1.0
    for _ in range(60):
        bW = (lo + hi) / 2
        a = (bW / 4) * (n + (bW / 4) * s)
        lo, hi = (bW, hi) if a.max() < 1 else (lo, bW)
    return lo, n, s


def main(P=8):
    L = Lattice(P)
    base = L.plaquette_partition()
    shapes = [(1, 1, 1, 1), (2, 1, 1, 1), (2, 2, 1, 1), (2, 1, 2, 1), (2, 2, 2, 1), (2, 2, 2, 2), (4, 2, 2, 2),
              (4, 4, 2, 2), (4, 4, 4, 4)]
    for c in shapes:
        cell = {}; lab_b = []
        for (x, mu, nu), ls in L.block_plaq:
            lab_b.append(cell.setdefault(tuple(x[i] // c[i] for i in range(4)), len(cell)))
        lab = np.array([lab_b[b] for b in base])
        b, n, s = ceiling(L, lab)
        print(f"cell {c}: first-order ceiling beta_W ~ {b:.4f} (max n {int(n.max())}, max internal-weighted sum {int(s.max())})",
              flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
