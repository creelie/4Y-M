"""Regular super-blocks: group the plaquette blocks by the coarse cell of their base point.

A block plaquette at x goes to the cell (x_0 // c_0, ..., x_3 // c_3).  For each cell shape we print the
distribution of the leading-order influence count n_j (distinct other blocks met by the plaquettes of j),
the links per super-block, and the number m_e of plaquettes lying inside the block of e (these cost an
internal contraction).  Periodic box of side P (a multiple of every c_i).

    python3 structured.py [P]
"""
import sys
from collections import Counter
import numpy as np
from geometry import Lattice


def counts(L, lab):
    n = np.zeros(len(L.active), int)
    for _, ls in L.plaqs:
        for j in ls:
            n[j] += len({lab[e] for e in ls if lab[e] != lab[j]})
    return n


def main(P=8):
    L = Lattice(P)
    base = L.plaquette_partition()
    shapes = [(1, 1, 1, 1), (2, 1, 1, 1), (2, 2, 1, 1), (1, 1, 2, 2), (2, 1, 2, 1), (2, 2, 2, 1), (2, 2, 2, 2),
              (4, 4, 1, 1), (4, 2, 2, 2), (4, 4, 2, 2), (4, 4, 4, 4), (P, P, 1, 1), (P, P, 2, 2), (P, 1, 1, 1)]
    for c in shapes:
        cell = {}
        lab_b = []
        for (x, mu, nu), ls in L.block_plaq:
            k = tuple(x[i] // c[i] for i in range(4))
            lab_b.append(cell.setdefault(k, len(cell)))
        lab = np.array([lab_b[b] for b in base])
        n = counts(L, lab)
        m = L.internal_counts(lab)
        size = Counter(Counter(lab).values())
        print(f"cell {c}: links/block {dict(size)}, max n {n.max()}, mean n {n.mean():.3f}, "
              f"share at max {np.mean(n == n.max()):.3f}, ceiling 4/max {4 / n.max():.3f}, "
              f"internal plaquettes per link: max {m.max()} mean {m.mean():.2f}", flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 8)
