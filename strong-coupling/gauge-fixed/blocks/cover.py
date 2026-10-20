"""Plaquette blocks for the line-field gauge of SU(2).

The line field sigma of the forest-gauge paper (period two, Table I) fixes exactly one link in
every plaquette, so every plaquette has three active links.  A *block cover* is a set of
plaquettes whose active triples partition the active links.  This script finds a block cover
invariant under translation by 2Z^4 (exact cover on the residues modulo two) and checks, on the
tori of side P = 4 and P = 6 (so that no plaquette wraps onto itself), that
  (a) the triples of the block plaquettes partition the active links;
  (b) every plaquette that is not a block plaquette meets three distinct blocks;
  (c) every active link lies in exactly five plaquettes that are not block plaquettes, so it has
      ten (plaquette, partner) incidences outside its block;
  (d) every link j and every block B that j influences (j shares a non-block plaquette with a link
      of B) satisfy |j - i| <= 2 for every link i of B, distances between link midpoints.
Usage: python3 cover.py"""
import itertools
import numpy as np

SIG = {"0000": 0, "0001": 1, "0010": 3, "0011": 3, "0100": 2, "0101": 1, "0110": 2, "0111": 0,
       "1000": 0, "1001": 2, "1010": 1, "1011": 2, "1100": 3, "1101": 3, "1110": 1, "1111": 0}

def sigma(x):
    return SIG["".join(str(c % 2) for c in x)]

def plaquettes(P):
    def sh(x, m):
        y = list(x); y[m] = (y[m] + 1) % P; return tuple(y)
    out = []
    for x in itertools.product(range(P), repeat=4):
        for m, n in itertools.combinations(range(4), 2):
            links = [(x, m), (sh(x, m), n), (sh(x, n), m), (x, n)]
            act = tuple(l for l in links if sigma(l[0]) != l[1])
            assert len(act) == 3
            out.append(((x, m, n), act))
    return out

def find_cover():
    """Exact cover of the active links of the torus of side two by plaquette triples."""
    plaq = plaquettes(2)
    items = sorted({l for _, a in plaq for l in a})
    cols = {l: [i for i, (_, a) in enumerate(plaq) if l in a] for l in items}
    sol, covered = [], set()
    def search():
        free = [l for l in items if l not in covered]
        if not free:
            return True
        l = min(free, key=lambda l: sum(1 for i in cols[l] if not set(plaq[i][1]) & covered))
        for i in cols[l]:
            a = plaq[i][1]
            if set(a) & covered:
                continue
            covered.update(a); sol.append(i)
            if search():
                return True
            sol.pop(); covered.difference_update(a)
        return False
    assert search()
    return sorted(plaq[i][0] for i in sol)

def check(cover, P):
    keys = set(cover)
    blk, plaq = {}, []
    for (x, m, n), act in plaquettes(P):
        inb = (tuple(c % 2 for c in x), m, n) in keys
        plaq.append((act, inb))
        if inb:
            for l in act:
                assert l not in blk, "two blocks share a link"
                blk[l] = (x, m, n)
    active = {l for act, _ in plaq for l in act}
    assert set(blk) == active, "(a) fails"
    outside = {l: 0 for l in active}
    for act, inb in plaq:
        if inb:
            continue
        assert len({blk[l] for l in act}) == 3, "(b) fails"
        for l in act:
            outside[l] += 1
    assert set(outside.values()) == {5}, "(c) fails"
    mid = lambda l: np.array(l[0], float) + 0.5*np.eye(4)[l[1]]
    def dist(a, b):
        d = mid(a) - mid(b); d = (d + P/2) % P - P/2; return np.linalg.norm(d)
    members = {}
    for l, b in blk.items():
        members.setdefault(b, []).append(l)
    dmax = 0.0
    for act, inb in plaq:
        if inb:
            continue
        for j in act:
            for i in act:
                if i != j:
                    dmax = max(dmax, max(dist(e, j) for e in members[blk[i]]))
    return len(active), len(members), dmax

if __name__ == "__main__":
    cover = find_cover()
    print("block plaquettes modulo two (base point x mod 2, plane {m, n}):")
    for c in cover:
        print("  ", "".join(map(str, c[0])), c[1:])
    for P in (4, 6):
        na, nb, dmax = check(cover, P)
        print("P = %d: %d active links in %d blocks; (a), (b), (c) hold; largest distance in (d) = %g"
              % (P, na, nb, dmax))
