#!/usr/bin/env python3
"""CPU diagnostic for the scalar crossover bounds, not a gauge RG test."""

from math import log


def run(u0, b, c, threshold, sign):
    u = u0
    steps = 0
    accumulated = 0.0
    while u > threshold:
        accumulated += 1.0 / u
        u += -b + sign * c / u
        steps += 1
    delta = b - c / threshold
    bound = 1 / threshold + log(u0 / threshold) / delta
    error = abs(steps - (u0 - u) / b)
    assert accumulated <= bound + 1e-10
    assert error <= c * bound / b + 1e-10
    assert threshold - b - c / threshold < u <= threshold
    return steps, u, accumulated, bound


if __name__ == "__main__":
    for initial in (100.0, 1000.0, 10000.0):
        for sign in (-1, 0, 1):
            n, endpoint, total, upper = run(initial, 1.0, 0.5, 10.0, sign)
            print(f"u0={initial:g} sign={sign:+d} n={n} "
                  f"exit={endpoint:.6f} sum_inverse={total:.6f} bound={upper:.6f}")
