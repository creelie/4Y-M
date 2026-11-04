# Passing a uniform physical gap to a reconstructed limit

**Proved (abstract Hilbert-space interface).** Let H_a >= 0 be self-adjoint
on Hilbert spaces with normalized vacuum Omega_a, H_a Omega_a=0. Suppose a
set D of labels carries vectors v_a(F) orthogonal to Omega_a, and that for
all F,G in D and t>=0 the inner products

    <v_a(F),v_a(G)>,   <v_a(F), exp(-t H_a)v_a(G)>

converge to the corresponding inner products of vectors v(F) and a
strongly continuous self-adjoint contraction semigroup exp(-tH) on a limit
Hilbert space. Suppose the span of the v(F) is dense in Omega-perp,
where H Omega=0, and for one m>0 independent of a,

    |<v_a(F), exp(-t H_a)v_a(G)>|
        <= exp(-m t) ||v_a(F)|| ||v_a(G)||.                        (1)

Then H restricted to Omega-perp has spectrum in [m,infinity).

*Proof.* Take limits in (1), using convergence of diagonal inner products.
The same hypothesis must hold for linear combinations, or D must already be
a vector space of labels; take it to be so. By density the resulting bound
extends to all vectors in Omega-perp and implies
||exp(-tH)|_{Omega-perp}|| <= exp(-mt). Omega-perp is invariant because the
semigroup is self-adjoint and fixes Omega. The spectral theorem then yields
spectrum(H|Omega-perp) subset [m,infinity). This also excludes any second
zero-energy vector. The bound m is in physical inverse-length units.

## A nontrivial-state witness

**Proved conditional deduction.** If for some F the centered norms obey
liminf_a ||v_a(F)||^2 >= c>0 and their inner products converge, then
v(F) is nonzero. Thus the preceding spectral conclusion is not merely a
statement about a zero-dimensional orthogonal complement. A lower bound on
a normalized correlation ratio without a lower bound on the chosen
renormalized norm does not alone give this conclusion.

## Remaining construction

**Open.** Nothing above proves convergence of vectors, inner products or
semigroups, OS reconstruction, Euclidean invariance, field locality, the
Yang--Mills identification, or (1). It also does not show that one chosen
observable is cyclic. The density assumption requires the full chosen
observable algebra, and the gap estimate must hold for its centered vectors,
not merely one two-point function. Positive-time flow observables can be
useful approximants, but their relationship to a dense local algebra and
reflection-positive time supports must be proved. This interface separates
spectral passage from those indispensable analytic tasks.
