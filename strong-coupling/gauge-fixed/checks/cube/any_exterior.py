"""Checks of the bound for every exterior (Proposition prop:anyexterior and
Corollary cor:anyexterior of the forest paper, Sec. ssec:largefields).

For boxes [0,l]^4 with l = 2, ..., 6 this program

  (a) enumerates the inner links and plaquettes and checks the counts
      n_L = 4l(l-1)^3, n_P = 6l^2(l-1)^2, the 24(l-1)^3 inner plaquettes
      that contain a boundary link, and that every plaquette containing an
      inner link is inner;
  (b) checks that the unit cubes x + [0,1]^4 with x in {1,...,l-2}^4 have
      only inner plaquettes, that l^4 - (l-2)^4 unit cubes meet the
      boundary, and that a plaquette lies in at most four unit cubes;
  (c) checks the constants: 24 pi^3 < e^7, the Haar mass of the cap of
      angle t in SU(2) against 8t^3/(3 pi^3) by quadrature, and the
      inequality 6x + 31 <= e^(x/4)/128 for x >= 44 on a fine grid and at
      x = 44 exactly in rational arithmetic, with the final bracket 1/64.

Usage: python3 any_exterior.py
"""

import itertools
import math
from fractions import Fraction as Fr

D = 4
checks = []


def check(name, ok):
    checks.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name)


def unit(i):
    return tuple(1 if k == i else 0 for k in range(D))


def add(a, b):
    return tuple(p + q for p, q in zip(a, b))


def box(l):
    def link_bd(x, m):
        return any(x[i] in (0, l) for i in range(D) if i != m)

    def plaq_bd(x, m, n):
        return any(x[i] in (0, l) for i in range(D) if i not in (m, n))

    def plinks(x, m, n):
        return [(x, m), (add(x, unit(m)), n), (add(x, unit(n)), m), (x, n)]

    sites = list(itertools.product(range(l + 1), repeat=D))
    links = [(x, m) for x in sites for m in range(D) if x[m] < l]
    inner = {e for e in links if not link_bd(*e)}
    plaqs = [(x, m, n) for x in sites for m in range(D) for n in range(m + 1, D)
             if x[m] < l and x[n] < l]
    ip = [p for p in plaqs if not plaq_bd(*p)]
    bp = sum(1 for p in ip if any(link_bd(*e) for e in plinks(*p)))
    closed = all(not plaq_bd(*p) for p in plaqs if any(e in inner for e in plinks(*p)))
    check(f"(a) l={l}: n_L={len(inner)}, n_P={len(ip)}, boundary-link plaquettes {bp}, "
          f"plaquettes with an inner link are inner",
          len(inner) == 4 * l * (l - 1) ** 3 and len(ip) == 6 * l * l * (l - 1) ** 2
          and bp == 24 * (l - 1) ** 3 and closed)
    cubes = list(itertools.product(range(l), repeat=D))
    meet = sum(1 for x in cubes if any(c in (0, l - 1) for c in x))
    ok = True
    incid = {}
    for x in cubes:
        for m in range(D):
            for n in range(m + 1, D):
                for o in itertools.product((0, 1), repeat=D):
                    if o[m] or o[n]:
                        continue
                    p = (add(x, o), m, n)
                    incid[p] = incid.get(p, 0) + 1
                    if not any(c in (0, l - 1) for c in x):
                        ok &= not plaq_bd(*p)
    check(f"(b) l={l}: {meet} unit cubes meet the boundary, interior cubes have inner "
          f"plaquettes, at most {max(incid.values())} cubes per plaquette",
          meet == l ** 4 - (l - 2) ** 4 and ok and max(incid.values()) <= 4)


def constants():
    check("(c) 24 pi^3 < e^7", 24 * math.pi ** 3 < math.exp(7))
    ok = True
    for t in (0.01, 0.05, 0.1, 0.25):
        n = 20000
        h = t / n
        mass = sum((2 / math.pi) * math.sin((k + 0.5) * h) ** 2 * h for k in range(n))
        ok &= mass >= 8 * t ** 3 / (3 * math.pi ** 3)
    check("(c) Haar mass of the cap of angle t is at least 8t^3/(3 pi^3)", ok)
    # e^11 > (8/3)^11 since e > 8/3
    check("(c) at x = 44: 6x + 31 = 295 <= (8/3)^11/128 < e^11/128",
          Fr(8, 3) ** 11 / 128 >= 295)
    ok = all(6 * x + 31 <= math.exp(x / 4) / 128 for x in (44 + k / 100 for k in range(20000)))
    check("(c) 6x + 31 <= e^(x/4)/128 on [44, 244]", ok)
    check("(c) 1/32 - 1/128 - 1/128 = 1/64", Fr(1, 32) - Fr(2, 128) == Fr(1, 64))
    check("(c) 48 beta/l <= beta^(1/4)/128 when l >= 6144 beta^(3/4)", Fr(48, 6144) == Fr(1, 128))


if __name__ == "__main__":
    for l in range(2, 7):
        box(l)
    constants()
    print(f"{sum(checks)} of {len(checks)} checks passed")
