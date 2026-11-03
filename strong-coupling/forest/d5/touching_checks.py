"""Computation checks for TOUCHING_CHART.md; not a Hessian certificate.

Run with Python's standard library. The geometric proof is in the note;
these finite enumerations check offsets and periodic separation cases.
"""
from itertools import product
from math import asin, sin, sqrt, pi

alpha = 2 * asin(1 / (2 * sqrt(2)))
a = 7 * alpha / (8 * sin(alpha))
assert a < 0.9561
assert 46 * a < 44
assert 2 * a / (1 - a) < 44
assert 8 * pi * 1e-8 + 46e-4 < 0.01

# One coordinate: a vertex is in a neighbouring block, and a plaquette
# base can be the vertex itself or one step behind it.
bases = {2 * block + offset - backward
         for block in (-1, 0, 1)
         for offset in (0, 1) for backward in (0, 1)}
assert bases == set(range(-3, 4))
assert len(list(product(bases, repeat=4))) == 2401

# On a periodic coarse cycle, shared fine-base boxes force distance <=3.
# Taking products gives the four-dimensional separation statement.
for side in range(1, 65):
    origin = {v % (2 * side) for v in bases}
    for y in range(side):
        shifted = {(2 * y + v) % (2 * side) for v in bases}
        if origin & shifted:
            assert min(y, side - y) <= 3

print(f'a = {a:.12f}')
print(f'46a = {46*a:.12f} < 44')
print(f'2a/(1-a) = {2*a/(1-a):.12f} < 44')
print('M = Dsep = 2401; 6M = 14406')
print('Offset enumeration and periodic separation (sides 1..64): passed')
print('No Hessian, label-domain, or island-threshold claim is checked here.')
