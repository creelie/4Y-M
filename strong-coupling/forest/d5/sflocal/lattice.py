"""Exact lattice model of the block-constrained quadratic form on a periodic
torus of fine side n (even).  Integer sparse matrices (scaled), no floats."""
import itertools
import numpy as np
import scipy.sparse as sp

def build(n):
    N = n // 2
    sites = list(itertools.product(range(n), repeat=4))
    sidx = {s: i for i, s in enumerate(sites)}
    def add(x, mu, k=1):
        y = list(x); y[mu] = (y[mu] + k) % n; return tuple(y)
    def tree(x, mu):
        o = [c % 2 for c in x]
        return o[mu] == 0 and all(o[nu] == 0 for nu in range(mu + 1, 4))
    links = [(x, mu) for x in sites for mu in range(4)]
    active = [l for l in links if not tree(*l)]
    aidx = {l: i for i, l in enumerate(active)}
    def block(x):
        return tuple(c // 2 for c in x)
    # path links
    path = {}
    for y in itertools.product(range(N), repeat=4):
        for mu in range(4):
            x = tuple(2 * c for c in y); x = add(x, mu)
            path[(y, mu)] = (x, mu)
    pathset = set(path.values())
    free = [l for l in active if l not in pathset]
    fidx = {l: i for i, l in enumerate(free)}
    coarse = [(y, mu) for y in itertools.product(range(N), repeat=4) for mu in range(4)]
    cidx = {k: i for i, k in enumerate(coarse)}
    # P2 = 2P : psi -> 2 phi ; L2 = 2L : gamma -> 2 phi
    rows, cols, vals = [], [], []
    for l in free:
        rows.append(aidx[l]); cols.append(fidx[l]); vals.append(2)
    Lr, Lc, Lv = [], [], []
    for (y, mu), d in path.items():
        r = aidx[d]
        Lr.append(r); Lc.append(cidx[(y, mu)]); Lv.append(16)
        # other private links of the pair
        for x in itertools.product(*[(2 * c, 2 * c + 1) for c in y]):
            if x[mu] % 2 == 1:
                l = (x, mu)
                if l != d:
                    rows.append(r); cols.append(fidx[l]); vals.append(-2)
        # S_y and S_{y+mu}: links (x,mu), x in B(y) or B(y+mu), o_mu=0, active
        for yy in (y, tuple((y[i] + (1 if i == mu else 0)) % N for i in range(4))):
            for x in itertools.product(*[(2 * c, 2 * c + 1) for c in yy]):
                if x[mu] % 2 == 0 and not tree(x, mu):
                    rows.append(r); cols.append(fidx[(x, mu)]); vals.append(-1)
    P2 = sp.csr_matrix((vals, (rows, cols)), shape=(len(active), len(free)), dtype=np.int64)
    L2 = sp.csr_matrix((Lv, (Lr, Lc)), shape=(len(active), len(coarse)), dtype=np.int64)
    # d : plaquettes x active links
    pr, pc, pv = [], [], []
    plaq = []
    for x in sites:
        for mu in range(4):
            for nu in range(mu + 1, 4):
                i = len(plaq); plaq.append((x, mu, nu))
                for (l, s) in (((x, mu), 1), ((add(x, mu), nu), 1), ((add(x, nu), mu), -1), ((x, nu), -1)):
                    if l in aidx:
                        pr.append(i); pc.append(aidx[l]); pv.append(s)
    dmat = sp.csr_matrix((pv, (pr, pc)), shape=(len(plaq), len(active)), dtype=np.int64)
    # V : coarse x active, scaled by 16
    vr, vc, vv = [], [], []
    for (y, mu) in coarse:
        i = cidx[(y, mu)]
        for x in itertools.product(*[(2 * c, 2 * c + 1) for c in y]):
            for l in ((x, mu), (add(x, mu), mu)):
                if l in aidx:
                    vr.append(i); vc.append(aidx[l]); vv.append(1)
    V16 = sp.csr_matrix((vv, (vr, vc)), shape=(len(coarse), len(active)), dtype=np.int64)
    V16.sum_duplicates()
    return dict(n=n, N=N, sites=sites, active=active, aidx=aidx, free=free, fidx=fidx,
                coarse=coarse, cidx=cidx, P2=P2, L2=L2, d=dmat, V16=V16, block=block,
                plaq=plaq, path=path)

def torus_dist(a, b, N):
    return max(min((a[i] - b[i]) % N, (b[i] - a[i]) % N) for i in range(4))
