"""Walsh parametrisation of the fibre {V(phi)=gamma}: exact integer matrices.
16*phi = P16 psi' + L16 gamma, psi' = (inside links, unnormalised Walsh
coordinates xi'_w, w != 0, of the eight private links of each pair)."""
import itertools, numpy as np, scipy.sparse as sp
from lattice import build, torus_dist

def tree(x, mu):
    o = [c % 2 for c in x]
    return o[mu] == 0 and all(o[nu] == 0 for nu in range(mu + 1, 4))

WALSH = [w for w in itertools.product((0, 1), repeat=3) if any(w)]

def build_walsh(n):
    M = build(n); N = M['N']
    aidx, active, coarse, cidx = M['aidx'], M['active'], M['coarse'], M['cidx']
    inside = [l for l in active if l[0][l[1]] % 2 == 0]
    iidx = {l: i for i, l in enumerate(inside)}
    nv = len(inside) + 7 * len(coarse)
    r, c, v = [], [], []
    lr, lc, lv = [], [], []
    for l in inside:
        r.append(aidx[l]); c.append(iidx[l]); v.append(16)
    for kk, (y, mu) in enumerate(coarse):
        trans = [nu for nu in range(4) if nu != mu]
        priv = []
        for x in itertools.product(*[(2 * q, 2 * q + 1) for q in y]):
            if x[mu] % 2 == 1:
                priv.append((tuple(x[nu] % 2 for nu in trans), (x, mu)))
        S = []
        for yy in (y, tuple((y[i] + (1 if i == mu else 0)) % N for i in range(4))):
            for x in itertools.product(*[(2 * q, 2 * q + 1) for q in yy]):
                if x[mu] % 2 == 0 and not tree(x, mu):
                    S.append((x, mu))
        base = len(inside) + 7 * kk
        for a, l in priv:
            row = aidx[l]
            for s in S:
                r.append(row); c.append(iidx[s]); v.append(-1)
            for wi, w in enumerate(WALSH):
                sgn = (-1) ** sum(ai * wj for ai, wj in zip(a, w))
                r.append(row); c.append(base + wi); v.append(16 * sgn)
            lr.append(row); lc.append(cidx[(y, mu)]); lv.append(16)
    P16 = sp.csr_matrix((v, (r, c)), shape=(len(active), nv), dtype=np.int64)
    P16.sum_duplicates()
    L16 = sp.csr_matrix((lv, (lr, lc)), shape=(len(active), len(coarse)), dtype=np.int64)
    blk = [M['block'](l[0]) for l in inside] + [coarse[k][0] for k in range(len(coarse)) for _ in range(7)]
    lam = [1] * len(inside) + [8] * (7 * len(coarse))   # |phi|^2 >= sum lam_i psi'_i^2
    kind = ['in'] * len(inside) + ['w'] * (7 * len(coarse))
    M.update(dict(P16=P16, L16=L16, inside=inside, blk=blk, lam=lam, kind=kind, nv=nv))
    return M
