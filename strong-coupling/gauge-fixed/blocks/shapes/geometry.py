"""Geometry of the line field and of unions of plaquette blocks (direction A3 of forest/DIRECTIONS.md).

Builds, on a periodic box of side P (P even), the line field F of Table "sigma", the active links, the block
plaquettes of Sec. "plaquette blocks", and for any partition of the active links into blocks the external
count n_j = number of pairs (plaquette p containing j, active link e of p) with e in a block other than the
block of j, and p not inside one block.  At small beta the Dobrushin sum of the block method is
alpha_j = (beta_W/4) * n_j * (1 + O(beta_W)), so 4/n_j is the leading-order ceiling of a block shape.
"""
import itertools
import numpy as np

# Table sigma: x mod 2 (digits x0x1x2x3) -> sigma(x)
SIG = {"0000": 0, "1000": 0, "0100": 2, "1100": 3, "0010": 3, "1010": 1, "0110": 2, "1110": 1,
       "0001": 1, "1001": 2, "0101": 1, "1101": 3, "0011": 3, "1011": 2, "0111": 0, "1111": 0}


def sigma(x):
    return SIG["".join(str(c % 2) for c in x)]


class Lattice:
    def __init__(self, P=4):
        self.P = P
        self.sites = list(itertools.product(range(P), repeat=4))
        self.links = [(x, mu) for x in self.sites for mu in range(4)]
        self.F = {(x, sigma(x)) for x in self.sites}
        self.active = [l for l in self.links if l not in self.F]
        self.aidx = {l: i for i, l in enumerate(self.active)}
        self.plaqs = []
        for x in self.sites:
            for mu, nu in itertools.combinations(range(4), 2):
                ls = [(x, mu), (self.sh(x, mu), nu), (self.sh(x, nu), mu), (x, nu)]
                nF = sum(l in self.F for l in ls)
                assert nF == 1, "every plaquette contains exactly one F link"
                self.plaqs.append(((x, mu, nu), [self.aidx[l] for l in ls if l not in self.F]))
        self.block_plaq = [p for p in self.plaqs if self.in_B(*p[0])]

    def sh(self, x, mu):
        y = list(x); y[mu] = (y[mu] + 1) % self.P; return tuple(y)

    @staticmethod
    def in_B(x, mu, nu):
        if (mu, nu) == (0, 1):
            return (x[0] - x[1]) % 2 == 0
        if (mu, nu) == (2, 3):
            return (x[2] - x[3]) % 2 == 0
        return False

    def plaquette_partition(self):
        lab = -np.ones(len(self.active), int)
        for b, (_, ls) in enumerate(self.block_plaq):
            for l in ls:
                assert lab[l] < 0; lab[l] = b
        assert (lab >= 0).all()
        return lab

    def external_counts(self, lab):
        n = np.zeros(len(self.active), int)
        for _, ls in self.plaqs:
            for j in ls:
                for e in ls:
                    if e != j and lab[e] != lab[j]:
                        n[j] += 1
        return n

    def internal_counts(self, lab):
        """m_e: plaquettes whose active links all lie in the block of e (and contain e)."""
        m = np.zeros(len(self.active), int)
        for _, ls in self.plaqs:
            if len({lab[l] for l in ls}) == 1:
                for l in ls:
                    m[l] += 1
        return m


if __name__ == "__main__":
    L = Lattice(4)
    lab = L.plaquette_partition()
    n = L.external_counts(lab)
    print("active links", len(L.active), "plaquettes", len(L.plaqs), "blocks", lab.max() + 1)
    print("plaquette blocks: n_j in", set(n.tolist()), " (paper: 10)")
