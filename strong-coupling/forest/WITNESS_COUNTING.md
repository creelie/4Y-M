# From local witnesses to joint rarity

**Proved (finite combinatorial reduction).** Let (A_p) be events indexed
by plaquettes and (B_e) events indexed by candidate bad links. Suppose
each B_e implies the union of A_p over a finite set N(e), with
|N(e)|<=M, and each plaquette belongs to at most D of the sets N(e), for
positive integers M,D. Assume
for every finite T that

    P(intersection_{p in T} A_p) <= q^|T|,     0<=q<=1.

Then, for every finite link set S,

    P(intersection_{e in S} B_e) <= min(1, (M q^(1/D))^|S|).

For each configuration in the intersection choose, for each e, a witness
p(e) in N(e). There are at most M^|S| such choices, and any resulting
image has at least |S|/D distinct plaquettes. Union bound and the assumed
joint plaquette estimate give the result. The estimate is useful only if
M q^(1/D)<1. It is annealed: it says nothing about conditioning on an
arbitrary exterior or the remaining labels.

**Open gauge application.** `FLAT_HOLONOMY_OBSTRUCTION.md` shows that the
absolute-link witness assumed by the earlier sketch need not exist.
`RELATIVE_CHART.md` already proves a local witness relative to R(c), and
`COARSE_LABELS.md` proves corresponding annealed finite-intersection
bounds with explicit coverage constants, conditional on its cited
chessboard input. The present generic deduction is an alternative, often
weaker counting route, not a newly established gauge witness. Even with
that witness, `RARE_LABEL_INFLUENCE.md` prevents treating the resulting
annealed rarity as uniform conditional environment contraction.
