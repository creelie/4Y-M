#!/usr/bin/env python3
"""Exact finite diagnostic for ENVIRONMENT_COMPARISON.md; no gauge input."""
from fractions import Fraction as F
from itertools import product


def inverse(a):
    n = len(a)
    b = [list(row) + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for j in range(n):
        k = next(i for i in range(j, n) if b[i][j])
        b[j], b[k] = b[k], b[j]
        pivot = b[j][j]
        b[j] = [x / pivot for x in b[j]]
        for i in range(n):
            if i != j:
                factor = b[i][j]
                b[i] = [x - factor*y for x, y in zip(b[i], b[j])]
    return [row[n:] for row in b]


def evaluate(weights):
    states = list(weights)
    n = len(states[0])
    total = sum(weights.values())
    mu = {x: w/total for x, w in weights.items()}
    def flip(x, i):
        return x[:i] + (1-x[i],) + x[i+1:]
    def cond(x, i):
        a = x[:i] + (0,) + x[i+1:]
        b = flip(a, i)
        return mu[b]/(mu[a]+mu[b])
    c = [[F(0) if i == j else max(abs(cond(x, i)-cond(flip(x, j), i))
          for x in states) for j in range(n)] for i in range(n)]
    assert max(map(sum, c)) < 1
    inv = inverse([[F(i == j)-c[i][j] for j in range(n)] for i in range(n)])
    for i in range(n):
        for j in range(n):
            cov = sum(mu[x]*x[i]*x[j] for x in states) - (
                sum(mu[x]*x[i] for x in states)*sum(mu[x]*x[j] for x in states))
            assert abs(cov) <= inv[i][j]/4
    return c, inv


def main():
    states = list(product((0, 1), repeat=3))
    # Asymmetric chain: orientation C_ij and C_ji differ because fields differ.
    weights = {x: F(2)**(x[0]+2*x[1]+3*x[2]+x[0]*x[1]+x[1]*x[2])
               for x in states}
    c, inv = evaluate(weights)
    assert c[0][1] != c[1][0]
    assert inv[0][1] != inv[1][0]
    # Hidden-bit mixture, q=1/4: pair covariance q^2/16 at any separation.
    q = F(1, 4)
    mean = (q/2+q)/2
    joint = ((q/2)**2+q**2)/2
    assert joint-mean**2 == q**2/16
    print('exact diagnostic passed; row sums:', [str(sum(row)) for row in c])
    print('asymmetric entries C01,C10:', c[0][1], c[1][0])
    print('hidden-phase covariance:', joint-mean**2)


if __name__ == '__main__':
    main()
