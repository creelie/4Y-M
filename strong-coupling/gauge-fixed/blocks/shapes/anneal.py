"""Search partitions of the active links into blocks of at most k links that minimise max_j n_j.

n_j (leading-order influence count) = sum over the plaquettes p containing j of the number of blocks, other
than the block of j, that meet p.  (A plaquette whose two other active links lie in one block B enters the
law of B through the single 1-Lipschitz function y_e y_e' and counts once.)  For plaquette blocks n_j = 10.

    python3 anneal.py K [P] [SWEEPS] [SEED]
"""
import sys
import numpy as np
from geometry import Lattice


def counts(L, lab):
    n = np.zeros(len(L.active), int)
    for _, ls in L.plaqs:
        for j in ls:
            n[j] += len({lab[e] for e in ls if lab[e] != lab[j]})
    return n


def anneal(L, k, sweeps=400, seed=0, init=None):
    rng = np.random.default_rng(seed)
    N = len(L.active)
    lab = (L.plaquette_partition() if init is None else init).copy()
    nb = int(lab.max()) + 1 + N // 3            # spare (empty) block labels
    size = np.bincount(lab, minlength=nb)
    inc = [[] for _ in range(N)]                 # plaquettes containing each link
    for q, (_, ls) in enumerate(L.plaqs):
        for l in ls:
            inc[l].append(q)
    P = [ls for _, ls in L.plaqs]

    def local_n(j, lab):
        return sum(len({lab[e] for e in P[q] if lab[e] != lab[j]}) for q in inc[j])

    def affected(j):
        s = {j}
        for q in inc[j]:
            s.update(P[q])
        return s

    n = np.array([local_n(j, lab) for j in range(N)])

    def cost(n):
        return n.max() + 0.02 * np.sort(n)[-N // 10:].mean() + 0.001 * n.mean()

    cur = cost(n)
    best = (n.max(), lab.copy())
    T0 = 0.5
    for sw in range(sweeps):
        T = T0 * (1 - sw / sweeps) + 1e-3
        for _ in range(N):
            j = rng.integers(N)
            # candidate target: block of a neighbour, or a fresh block
            if rng.random() < 0.9:
                q = inc[j][rng.integers(len(inc[j]))]
                e = P[q][rng.integers(3)]
                t = lab[e]
            else:
                t = int(np.argmin(size))
            s = lab[j]
            if t == s or size[t] >= k:
                continue
            aff = affected(j)
            old = {i: n[i] for i in aff}
            lab[j] = t
            for i in aff:
                n[i] = local_n(i, lab)
            new = cost(n)
            if new <= cur or rng.random() < np.exp((cur - new) / T):
                cur = new; size[s] -= 1; size[t] += 1
                if n.max() < best[0]:
                    best = (n.max(), lab.copy())
            else:
                lab[j] = s
                for i, v in old.items():
                    n[i] = v
    return best


if __name__ == "__main__":
    k = int(sys.argv[1]); P = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    sweeps = int(sys.argv[3]) if len(sys.argv) > 3 else 200; seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    L = Lattice(P)
    print("plaquette blocks: max n_j =", counts(L, L.plaquette_partition()).max(), flush=True)
    m, lab = anneal(L, k, sweeps, seed)
    n = counts(L, lab)
    sizes = np.bincount(lab); sizes = sizes[sizes > 0]
    print(f"k={k} P={P}: max n_j = {n.max()}  mean n_j = {n.mean():.3f}  block sizes {np.bincount(sizes)[1:]}"
          f"  leading-order ceiling 4/max n = {4 / n.max():.3f}")
    np.save(f"best_k{k}_P{P}_s{seed}.npy", lab)
