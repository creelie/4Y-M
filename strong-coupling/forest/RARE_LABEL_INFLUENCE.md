# Rare labels do not imply uniform environment contraction

**Proved (finite counterexample).** For 0<e<1/4 put a strictly positive
law on two binary labels by

    P(1,1)=e,  P(1,0)=P(0,1)=e²/2,
    P(0,0)=1-e-e².

Each label has probability e+e²/2 of being one, tending to zero. Yet the
conditional influence of label 2 on label 1 is

    P(B1=1 | B2=1)-P(B1=1 | B2=0)
      = 1/(1+e/2) - (e²/2)/(1-e-e²/2) -> 1.

Thus arbitrarily rare bad labels do not provide a uniform single-label
Dobrushin contraction, even for a strictly positive finite law. This is
compatible with unconditional covariance tending to zero: conditioning on
a rare neighbour can reveal a highly correlated island.

**Sketch (consequence for D/I).** An environment method based only on
unconditional rarity needs a cluster-level or graded treatment. A
single-label criterion such as `ENVIRONMENT_COMPARISON.md` instead needs
uniform conditional odds for all frozen environments. One possible
intermediate target is a bound on the conditional probability that a bad
cluster reaches a prescribed boundary, under a *specified calm collar*;
the collar probability and its correlation with physical observables must
then be controlled separately. This example neither proves nor disproves
such an estimate for constrained SU(2) fibres.
