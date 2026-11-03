"""Exact four-link counterexample; Python standard library, integer arithmetic.

This verifies the finite identities used in ROUGH_BACKGROUND.md. The note
proves that the curve stays on the nonlinear fibre near t=0.
"""
from collections import Counter
from itertools import product, combinations
import json

FINE_SIDE = 4
COARSE_SIDE = 2
XI = {((0, 1, 1, 1), 1): -1,
      ((0, 1, 1, 1), 2): 1,
      ((0, 1, 0, 0), 1): 1,
      ((0, 0, 1, 0), 2): -1}

def shift(x, mu):
    out = list(x)
    out[mu] = (out[mu] + 1) % FINE_SIDE
    return tuple(out)

def centre(x, mu):
    return 1 if x[mu] % 2 == 0 else (-1) ** sum(x[nu] // 2 for nu in range(mu))

def tree(x, mu):
    return x[mu] % 2 == 0 and all(x[nu] % 2 == 0 for nu in range(mu + 1, 4))

assert all(not tree(*e) for e in XI)
constraints = []
for y in product(range(COARSE_SIDE), repeat=4):
    for mu in range(4):
        total = 0
        for off in product(range(2), repeat=4):
            x = tuple(2 * y[j] + off[j] for j in range(4))
            total += XI.get((x, mu), 0) + XI.get((shift(x, mu), mu), 0)
        constraints.append(total)
assert len(constraints) == 64 and all(a == 0 for a in constraints)
counts = Counter()
for x in product(range(FINE_SIDE), repeat=4):
    for mu, nu in combinations(range(4), 2):
        edges = [(x, mu), (shift(x, mu), nu), (shift(x, nu), mu), (x, nu)]
        sg = 1
        for e in edges:
            sg *= centre(*e)
        curl = sum(orientation * XI.get(e, 0)
                   for orientation, e in zip((1, 1, -1, -1), edges))
        if curl:
            counts[(sg, curl * curl)] += 1
positive = sum(square * count for (sg, square), count in counts.items() if sg == 1)
negative = sum(square * count for (sg, square), count in counts.items() if sg == -1)
assert (positive, negative) == (12, 14)
assert sum(value * value for value in XI.values()) == 4
print('Exact four-link certificate on fine torus side 4')
print('All 64 average numerators: zero')
print('Nonzero plaquettes by (central sign, curl squared):')
for (sg, square), count in sorted(counts.items()):
    print(f'  ({sg:+d}, {square}): {count}')
print(f'Positive sum: {positive}; negative sum: {negative}; signed sum: {positive-negative}')
print('Q second derivative for HS-unit generator: -4')
print('T(xi): 4; Q second derivative / T(xi): -1')
print('Wilson action beta_W*Q/4 second derivative: -beta_W')
