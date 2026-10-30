"""Independent exact check of the certificate for Proposition prop:blockgapsharp of the forest paper:
for every set of blocks W and every phi in Phi_W (Def. def:blockgauge, Thm. thm:blockgap),
    sum_e phi_e^2 <= (1/c) Q(phi),   c = 537/5000.

The certificate (file certificate_c537_5000.json) gives c, multipliers l (4 x 49 rationals) and a symmetric
245 x 245 matrix G = G_upper/den on the window S(y) = B(y) u B(y+e_0) u ... u B(y+e_3). The program

  (a) builds, from the definitions of the paper alone, the translation-invariant symmetric kernels of
      Q(phi) = sum_p (d phi)_p^2,  sum_e phi_e^2,  sum_{y,mu} V_{y,mu}(phi) L_{y,mu}(phi)  and  sum_y G(phi|S(y)),
      with L_{y,mu}(phi) = sum_b l[mu][b] phi_{(y,b)}, and checks in exact rational arithmetic that
      Q - c sum phi^2 - sum V L = sum G  entry by entry;
  (b) checks that G is positive definite: all leading principal minors of the integer matrix den * G are
      positive (fraction-free Gaussian elimination);
  (c) checks the reduction: on Phi_W every term V_{y,mu} L_{y,mu} vanishes (V = 0 for y in W, and L_{y,mu}
      depends only on links with base in B(y), where phi = 0 for y not in W) on random fields in exact arithmetic.

Local links of a block: offsets o in {0,1}^4 (itertools.product order) times mu = 0..3, tree links removed; window
index = slot * 49 + local index, slots (0, e_0, e_1, e_2, e_3).

Usage: python3 blockgap_certificate.py certificate_c537_5000.json
"""

import itertools
import json
import sys
from fractions import Fraction as Fr

checks = []


def check(name, ok):
    checks.append((name, bool(ok)))
    print(("PASS " if ok else "FAIL ") + name, flush=True)


def is_tree(o, mu):
    return o[mu] == 0 and all(o[n] == 0 for n in range(mu + 1, 4))


LOCAL = [(o, mu) for o in itertools.product((0, 1), repeat=4) for mu in range(4) if not is_tree(o, mu)]
INDEX = {l: i for i, l in enumerate(LOCAL)}
NL = len(LOCAL)
E = [tuple(int(i == mu) for i in range(4)) for mu in range(4)]
SLOTS = [(0, 0, 0, 0)] + E


def vadd(a, b):
    return tuple(p + q for p, q in zip(a, b))


def vsub(a, b):
    return tuple(p - q for p, q in zip(a, b))


def active(x, mu):
    """(cell, local index) of the link (x, mu), or None for a tree link."""
    cell = tuple(v // 2 for v in x)
    i = INDEX.get((tuple(v % 2 for v in x), mu))
    return None if i is None else (cell, i)


class Kernel:
    """A translation-invariant symmetric bilinear form sum_{y} B(T_y phi), stored by its matrix entries
    K[(delta, a, b)] = coefficient of phi_{(0,a)} phi_{(delta,b)} in the symmetric matrix of the form."""

    def __init__(self):
        self.K = {}

    def add_local(self, vec_pairs):
        """Add sum_y of the rank-one or general local form sum_{i,j} w_ij phi_{(c_i,a_i)} phi_{(c_j,a_j)};
        vec_pairs: list of ((c_i,a_i), (c_j,a_j), w_ij) over ordered pairs (the symmetric matrix of a local form)."""
        for (c1, a), (c2, b), w in vec_pairs:
            k = (vsub(c2, c1), a, b)
            self.K[k] = self.K.get(k, 0) + w


def rank_one(v):
    return [(l1, l2, s1 * s2) for l1, s1 in v for l2, s2 in v]


def sym_product(u, v):
    out = []
    for l1, s1 in u:
        for l2, s2 in v:
            out.append((l1, l2, Fr(s1 * s2, 2) if isinstance(s1 * s2, int) else s1 * s2 / 2))
            out.append((l2, l1, Fr(s1 * s2, 2) if isinstance(s1 * s2, int) else s1 * s2 / 2))
    return out


def main(path):
    cert = json.load(open(path))
    c = Fr(cert["c"])
    lmul = [[Fr(x) for x in cert["L"][mu * NL:(mu + 1) * NL]] for mu in range(4)]
    den = int(cert["den"])
    NW = 5 * NL
    Gint = [[0] * NW for _ in range(NW)]
    for i, row in enumerate(cert["G_upper"]):
        for k, v in enumerate(row):
            Gint[i][i + k] = Gint[i + k][i] = int(v)
    check(f"certificate: c = {c}, {NL} active links per block, window of {NW} links, "
          f"{len(cert['L'])} multipliers", NL == 49 and len(cert["L"]) == 4 * NL and len(cert["G_upper"]) == NW)

    # Q: every plaquette, counted once per translation class of blocks, as the plaquettes with base in B(0)
    KQ = Kernel()
    for o in itertools.product((0, 1), repeat=4):
        for mu, nu in itertools.combinations(range(4), 2):
            x = o
            sides = [(x, mu, 1), (vadd(x, E[mu]), nu, 1), (vadd(x, E[nu]), mu, -1), (x, nu, -1)]
            v = [(active(s, d), sg) for s, d, sg in sides if active(s, d) is not None]
            KQ.add_local(rank_one(v))
    KS = Kernel()
    KS.add_local([((SLOTS[0], a), (SLOTS[0], a), 1) for a in range(NL)])
    # sum_mu V_{0,mu} L_{0,mu}, translated
    KV = Kernel()
    for mu in range(4):
        vrow = {}
        for o in itertools.product((0, 1), repeat=4):
            for s in (o, vadd(o, E[mu])):
                l = active(s, mu)
                if l is not None:
                    vrow[l] = vrow.get(l, 0) + Fr(1, 16)
        lrow = [((SLOTS[0], b), lmul[mu][b]) for b in range(NL) if lmul[mu][b] != 0]
        KV.add_local(sym_product(list(vrow.items()), lrow))
    KG = Kernel()
    win = [(SLOTS[s], a) for s in range(5) for a in range(NL)]
    KG.add_local([(win[i], win[j], Fr(Gint[i][j], den)) for i in range(NW) for j in range(NW) if Gint[i][j]])

    keys = set(KQ.K) | set(KS.K) | set(KV.K) | set(KG.K)
    bad = [k for k in keys if KQ.K.get(k, 0) - c * KS.K.get(k, 0) - KV.K.get(k, 0) != KG.K.get(k, 0)]
    check(f"(a) Q - c sum phi^2 - sum V L = sum G as translation-invariant kernels, exactly, "
          f"on all {len(keys)} entries (failing: {len(bad)})", not bad)

    # (b) positive definiteness by fraction-free elimination
    A = [row[:] for row in Gint]
    prev, ok = 1, True
    for k in range(NW):
        if A[k][k] <= 0:
            ok = False
            break
        for i in range(k + 1, NW):
            for j in range(k + 1, NW):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    check(f"(b) all {NW} leading principal minors of den * G are positive", ok)

    # (c) the multiplier terms vanish on Phi_W (W = a 2 x 1 x 1 x 1 row of blocks), exactly
    import random
    rnd = random.Random(5)
    W = [(0, 0, 0, 0), (1, 0, 0, 0)]
    var = [(cell, a) for cell in W for a in range(NL)]
    vrows = []
    for y in W:
        for mu in range(4):
            r = {}
            for o in itertools.product((0, 1), repeat=4):
                x = tuple(2 * y[i] + o[i] for i in range(4))
                for s in (x, vadd(x, E[mu])):
                    l = active(s, mu)
                    if l is not None and l[0] in W:
                        r[l] = r.get(l, 0) + Fr(1, 16)
            vrows.append(r)
    # exact projection of a random rational field onto {V_{y,mu} = 0, y in W} (Gaussian elimination)
    phi = {v: Fr(rnd.randint(-9, 9), rnd.randint(1, 9)) for v in var}
    rows = [dict(r) for r in vrows]
    piv = []
    for r in rows:
        for v, p in piv:
            if v in r and r[v] != 0:
                f = r[v] / p[v]
                for u, w in p.items():
                    r[u] = r.get(u, 0) - f * w
        nz = [v for v in r if r[v] != 0]
        if nz:
            piv.append((nz[0], r))
    for v, p in reversed(piv):
        # set phi[v] so that p(phi) = 0
        s = sum(w * phi[u] for u, w in p.items() if u != v)
        phi[v] = -s / p[v]
    for v, p in piv:
        s = sum(w * phi[u] for u, w in p.items())
        assert s == 0
    allzero = True
    for y in itertools.product(range(-1, 3), (-1, 0, 1), (-1, 0, 1), (-1, 0, 1)):
        for mu in range(4):
            V = Fr(0)
            for o in itertools.product((0, 1), repeat=4):
                x = tuple(2 * y[i] + o[i] for i in range(4))
                for s in (x, vadd(x, E[mu])):
                    l = active(s, mu)
                    if l is not None:
                        V += Fr(1, 16) * phi.get(l, 0)
            Lv = sum(lmul[mu][b] * phi.get((y, b), 0) for b in range(NL))
            if V * Lv != 0:
                allzero = False
    check("(c) on a field of Phi_W (two blocks) every term V_{y,mu} L_{y,mu} vanishes, exactly", allzero)
    print(f"{sum(ok for _, ok in checks)} of {len(checks)} checks passed")


if __name__ == "__main__":
    main(sys.argv[1])
