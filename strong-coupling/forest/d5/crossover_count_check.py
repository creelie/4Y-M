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


def propagated_sum(n, threshold, b, amplification, exponent):
    # Traverse backwards from the crossover so that weight is the product
    # of all later multipliers. The final error has an empty product.
    weight = 1.0
    total = 0.0
    for m in range(1, n + 1):
        u = threshold + m * b
        total += weight / u**exponent
        weight *= 1 + amplification / u
    return total


if __name__ == "__main__":
    for initial in (100.0, 1000.0, 10000.0):
        for sign in (-1, 0, 1):
            n, endpoint, total, upper = run(initial, 1.0, 0.5, 10.0, sign)
            print(f"u0={initial:g} sign={sign:+d} n={n} "
                  f"exit={endpoint:.6f} sum_inverse={total:.6f} bound={upper:.6f}")
    for exponent in (1.5, 2.0, 2.5):
        values = [propagated_sum(n, 10.0, 1.0, 1.0, exponent)
                  for n in (100, 1000, 10000)]
        print(f"p={exponent:g} propagated=" + ", ".join(f"{x:.6f}" for x in values))
