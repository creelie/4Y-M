"""Finite checks of the overlap constants; proof is in COARSE_LABELS.md."""
from math import ceil
from itertools import combinations, product

for radius, H1 in [(3,4),(5,6)]:
    worst=0
    for side in range(1,65):
        fine=2*side
        boxes=[{(2*y+j)%fine for j in range(-radius,radius+1)} for y in range(side)]
        for x in range(fine):
            count=sum(x in box for box in boxes)
            worst=max(worst,count)
            assert count<=H1
    print(f'radius {radius}: one-coordinate overlap <= {H1}; attained {worst}')
    print(f'  four-coordinate M={(2*radius+1)**4}, H={H1**4}')

# Exhaustive tiny incidence examples: the deterministic coverage inequality
# is what permits the rare-event union bound. This does not test chessboard.
for labels in range(1,4):
    subsets=[set(i for i in range(3) if mask & (1<<i)) for mask in range(1,8)]
    for family in product(subsets,repeat=labels):
        H=max(sum(x in q for q in family) for x in range(3))
        k=ceil(labels/H)
        for mask in range(8):
            bad={i for i in range(3) if mask & (1<<i)}
            if all(q & bad for q in family):
                assert len(set.union(*family) & bad)>=k
print('Exhaustive tiny coverage implication: passed')
