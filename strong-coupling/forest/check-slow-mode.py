"""Exact finite-torus controls for the Gaussian comparison, not a Wilson gap proof."""
from itertools import product
import json


def check(block, coarse_side):
    side = block * coarse_side

    def value(base):
        return int(any(x % block for x in base[1:]))

    norm_squared = curl_energy = divergence_energy = 0
    for base in product(range(side), repeat=4):
        v = value(base)
        norm_squared += v
        for direction in range(1, 4):
            neighbor = list(base)
            neighbor[direction] = (neighbor[direction] + 1) % side
            curl_energy += (v - value(neighbor)) ** 2
        if sum(x % block != 0 for x in base) > 1:
            neighbor = list(base)
            neighbor[0] = (neighbor[0] - 1) % side
            divergence_energy += (v - value(neighbor)) ** 2
    assert norm_squared == (block**4 - block) * coarse_side**4
    assert curl_energy == 6 * block * coarse_side**4
    assert divergence_energy == 0
    return dict(block=block, coarse_side=coarse_side, norm_squared=norm_squared,
                curl_energy=curl_energy, divergence_energy=divergence_energy)


if __name__ == '__main__':
    print(json.dumps([check(block, 3) for block in (3, 4, 5)], indent=2))
