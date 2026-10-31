"""Group plaquette blocks into super-blocks of at most m plaquette blocks (3m links), minimising max_j n_j.

Same count as anneal.py; the moves act on whole plaquette blocks, so that the three links of a block plaquette
always stay together.  Simulated annealing on a periodic box of side P.

    python3 merge.py M [P] [SWEEPS] [SEED]
"""
import sys
import numpy as np
from geometry import Lattice


def main(m, P=4, sweeps=2000, seed=0):
    rng = np.random.default_rng(seed)
    L = Lattice(P)
    base = L.plaquette_partition()               # link -> plaquette block
    nB = base.max() + 1
    plist = [ls for _, ls in L.plaqs]
    # n_j depends on group labels of links; precompute for each link the plaquettes (as link triples)
    inc = [[] for _ in range(len(L.active))]
    for q, ls in enumerate(plist):
        for l in ls:
            inc[l].append(q)
    links_of = [np.where(base == b)[0] for b in range(nB)]
    nbrs = [set() for _ in range(nB)]
    for ls in plist:
        bs = {base[l] for l in ls}
        for b in bs:
            nbrs[b] |= bs - {b}
    nbrs = [sorted(s) for s in nbrs]

    grp = np.arange(nB)                           # plaquette block -> group
    gsize = np.ones(2 * nB, int); gsize[nB:] = 0

    def nj(j):
        gj = grp[base[j]]
        return sum(len({grp[base[e]] for e in plist[q] if grp[base[e]] != gj}) for q in inc[j])

    n = np.array([nj(j) for j in range(len(L.active))])

    def cost(n):
        s = np.sort(n)
        return n.max() + 0.05 * s[-len(n) // 8:].mean() + 0.002 * n.mean()

    cur = cost(n); best = (n.max(), n.mean(), grp.copy())
    for sw in range(sweeps):
        T = 0.6 * (1 - sw / sweeps) + 0.01
        for _ in range(nB):
            b = rng.integers(nB)
            t = grp[nbrs[b][rng.integers(len(nbrs[b]))]] if rng.random() < 0.95 else nB + rng.integers(nB)
            s = grp[b]
            if t == s or gsize[t] >= m:
                continue
            aff = set()
            for l in links_of[b]:
                for q in inc[l]:
                    aff.update(plist[q])
            # links in the old and new groups that share plaquettes with b's links are covered by aff
            old = {i: n[i] for i in aff}
            grp[b] = t
            for i in aff:
                n[i] = nj(i)
            new = cost(n)
            if new <= cur or rng.random() < np.exp((cur - new) / T):
                cur = new; gsize[s] -= 1; gsize[t] += 1
                if (n.max(), n.mean()) < best[:2]:
                    best = (n.max(), n.mean(), grp.copy())
            else:
                grp[b] = s
                for i, v in old.items():
                    n[i] = v
    mx, mean, g = best
    sizes = np.bincount(np.bincount(g)[np.bincount(g) > 0])
    print(f"m={m} P={P} seed={seed}: max n_j = {mx}, mean n_j = {mean:.3f}, groups by size {dict(enumerate(sizes)) }"
          f", leading-order ceiling 4/max = {4 / mx:.3f}", flush=True)
    np.save(f"merge_m{m}_P{P}_s{seed}.npy", g)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(int(a[0]), *(int(x) for x in a[1:]))
