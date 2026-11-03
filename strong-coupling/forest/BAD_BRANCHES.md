# Bad independent coordinates can have several path-link branches

**Proved (explicit example).** The unique-path assertion of
`INDEPENDENT_LABELS.md` is a good-regime result. In the bad regime, even
an abelian-embedded SU(2) example has at least three path-link solutions
with the same seven other outgoing private links and the same coarse
average. Thus the good graph cannot be extended globally just by dropping
its radius assumption.

## 1. A scalar family

Freeze all inside-block links to I and all but one pair's private links.
Use the diagonal subgroup diag(exp(iq),exp(-iq)); q denotes the eigenangle,
so its geodesic length from I near zero is sqrt(2)|q|. Prescribe one of
the seven non-path private angles to -a=-1/2, and the other six to
b=49/100. Let q be the path-link angle. The two offset-zero paths give
zero logs, while each of the other seven private links occurs twice.
When the cut-off is nonzero all relative principal angles are within
its log region. The average angle is consequently

    f(q) = q + chi(q)[-a+6b-7q]/8,
    chi(q) = min(1,max(0,2-2m(q))),
    m(q) = max_j 2 sqrt(2) |sin((theta_j-q)/2)|,
    theta_1=-1/2,   theta_2=...=theta_7=49/100.

If chi=0 the average is the path itself, as usual. There is no wrapping
ambiguity for the values used below.

At q=0 the unique extreme relative angle is -a. Also
1/2<2 sqrt(2) sin(1/4)<1, so the maximum and cut-off have smooth formulas
near zero. Set gamma=f(0). Differentiation gives

    f'(0) = -3/4 + sqrt(2)[(7/2)sin(1/4)-(61/100)cos(1/4)] < 0.

For an elementary strict bound, use sin(1/4)<1/4,
cos(1/4)>31/32 and 1<sqrt(2)<3/2 separately on the positive and negative
terms: f'(0)<-91/3200. The transition condition follows from
sqrt(2)>7/5 and sin(1/4)>1/4-(1/4)^3/6 for the lower bound, and the
same simple upper bounds for the upper bound.

## 2. Three solutions (proved)

For small epsilon>0, f(-epsilon)>gamma>f(epsilon). At q=-1 and q=1 the
maximum chordal distance exceeds one, so the cut-off is zero and
f(-1)=-1, f(1)=1. Since 0<gamma<1, continuity gives a preimage of gamma
in (-1,-epsilon) and another in (epsilon,1), in addition to q=0.
All three represent distinct SU(2) path links. The average of every other
pair is unchanged, since these links are private.

There is also a smooth turning point: throughout 0<=q<=1/5 the unique
maximiser is the -1/2 path and the cut-off stays strictly in its transition
region. At q=1/5, m>9/10 and m<1, using the same sine bounds and
sqrt(2)<10/7. Thus chi<1/5 and

    f'(1/5) > 1-7/40-3*(13/100) > 0.

Together with f'(0)<0, continuity gives a point with f'=0 strictly between
0 and 1/5. The full SU(2) path derivative has a zero axial direction there.
Thus one cannot obtain a global positive lower bound for the path-link
Jacobian from this parametrisation. This does not contradict the paper's
uniform **eight-private-link** submersion bound: that bound allows all
eight links to move, while here seven are frozen.

## 3. Scope

**Proved:** a bad-side description using independent non-path coordinates
must retain all admissible path branches, or use different coordinates.
The single-graph determinant from the good regime does not suffice.

**Sketch / open:** after integrating over the allowed bad non-path links,
these branch effects may regularise. The example does not prove that a
bad-side Wasserstein coefficient is unbounded, nor rule out an island
estimate. It identifies a concrete term that a valid such estimate must
include. The environment comparison and D8 remain open.

**Computation:** `d5/bad_branch.py` checks the rational sign inequalities
and finds the three roots numerically. The intermediate-value argument,
not numerical root finding, proves their existence. Output:
`d5/bad_branch.out`.
