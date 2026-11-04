#!/usr/bin/env python3
"""Floating quadrature check of the SU(2) radial hazard bound."""
from math import cos, exp, expm1, pi, sin


def density(s, k, alpha):
    z = k * sin(alpha) * sin(s)
    base = k * (cos(alpha) * cos(s) - 1)
    angular = exp(base) if abs(z) < 1e-12 else exp(base + z) * (-expm1(-2*z)) / (2*z)
    return sin(s)**2 * angular


def integrate(f, a, b, n=8192):
    h = (b-a)/n
    return h/3 * (f(a)+f(b)+4*sum(f(a+h*i) for i in range(1,n,2))
                  + 2*sum(f(a+h*i) for i in range(2,n,2)))


if __name__ == '__main__':
    largest = (0.0, None)
    for k in (0, 1, 10, 100):
        for alpha in (0, pi/4, pi/2, 3*pi/4, pi):
            for r in (0.01, 0.3, 1.5):
                f = lambda s: density(s,k,alpha)
                ratio = f(r) / ((k+1)*integrate(f,r,pi))
                largest = max(largest, (ratio,(k,alpha,r)))
    print('maximum f / ((kappa+1) tail):', largest)
    print('proved upper bound:', exp(1)/cos(1)**2)
    assert largest[0] < exp(1)/cos(1)**2
