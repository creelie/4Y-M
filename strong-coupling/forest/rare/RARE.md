# Rare bad fibres (direction C2)

`../DIRECTIONS.md` (C2) suggested that "bad fibres are rare and harmless" might follow from disagreement percolation
and the paper's Propositions `anyexterior` and `calmrare`. The two halves come out differently:

- **Rare: proved here, by chessboard estimates rather than disagreement percolation.** Under the SU(2) Wilson measure
  on an even torus, the sites with a large plaquette are dominated, in the sense of all finite intersections, by a
  Bernoulli field of density 6z. Here z ≤ exp(−β_Wρ²/16 + 6 log β_W + 25). So the bad clusters have exponential
  tails above an explicit β_W (table below). This is stronger than Corollary `anyexterior`, which bounds only the
  density of bad cubes.
- **Harmless: not obtainable from disagreement percolation at weak coupling.** The single-link and single-window
  conditional laws under two neighbouring configurations, even good ones, are at total-variation distance of order one
  (§3). The harmless half needs a contraction in a Wasserstein or chart metric with bad regions split off (§4). That is
  the open composition problem of the cube estimate (D), not a combinatorial lemma.

## 1. The lemma

Notation follows `../main.tex` §"large fields":
- G = SU(2);
- s_p = 1 − ½ Tr U_p = δ(U_p, I)²/4, where δ is the Hilbert–Schmidt distance;
- the Wilson measure is ∝ exp(−β_W Σ_p s_p) ∏ dU on the periodic lattice Λ = (ℤ/Lℤ)⁴ with L even and N = L⁴ sites.

Fix ρ ∈ (0, 1]:
- A plaquette is **bad** if δ(U_p, I) > ρ.
- A site x is **bad** if one of the six plaquettes with base point x (x, x+e_μ, x+e_ν, x+e_μ+e_ν) is bad.

Every bad plaquette lies in a unit cube y + [0,1]⁴ whose base y is at distance ≤ 1, in ℓ∞, from a bad site. So a
statement about bad sites controls bad unit cubes up to a ∗-neighbourhood.

**Lemma (rare large plaquettes).** Let

    log z = −β_W ρ²/16 + 6 − 4 log h(t),   t = (8β_W)^{−1/2},   h(t) = (t − sin t cos t)/π,

where h(t) is the Haar mass of {U : angle(U) ≤ t}. Then:

(i) For every finite set S of sites, P(every x ∈ S is bad) ≤ (6z)^{|S|}.

(ii) Let 𝒢 be any adjacency on sites of maximal degree D, e.g. D = 8 for nearest neighbours or D = 80 for ℓ∞
neighbours, and put q = 6eDz. If q < 1, the 𝒢-cluster of bad sites containing a given site has
P(size ≥ n) ≤ qⁿ/(1 − q).

(iii) Uniformly: −4 log h(t) ≤ 6 log β_W + 4 log(3π/2) + 12 log(2√2) + o(1), so log z ≤ −β_Wρ²/16 + 6 log β_W + 24.7.

**Proof.**

*Reflection positivity.* Take a reflection θ of the torus in a pair of hyperplanes {x_μ = k} ∪ {x_μ = k + L/2}
through sites. It splits the torus into two closed halves that share the links in the planes. Every plaquette lies in
one closed half: a plaquette with a link crossing a plane has all four vertices on that link's side or on the plane.
Plaquettes inside a plane depend only on links that θ fixes, so they can be split with weight ½ to each side.
So the action splits as S = S₊ + θS₊ with no cross term. For F measurable with respect to the closed half,

    E[F · θF] = ∫ ∏_{plane links} dU · |∫ F e^{−β_W S₊} ∏_{+ links} dU|² ≥ 0.

This is site-reflection positivity, the same argument as Osterwalder–Seiler [OSei78] for the plane through sites. δ(U_p, I)
is invariant under orientation reversal U_p ↦ U_p^{−1}, so the reflected events are again "plaquette bad".

*Chessboard estimate.* Take the block B₀ = [0,1]⁴ ∩ ℤ⁴ (side 1) and its images B_y = y + B₀ under compositions of the
reflections above, one for each site y. The torus side L is a multiple of 2. For events A_y measurable with respect
to the links of B₀, the chessboard estimate (Fröhlich–Israel–Lieb–Simon [FILS78]; Biskup [Bis09, §5] in this form, with
site reflections and blocks of side B dividing L/2) gives

    P(⋂_y θ_y A_y) ≤ ∏_y z(A_y),   z(A) = P(⋂_y θ_y A)^{1/N},

with z(Ω) = 1. Here θ_y is the composition of reflections mapping B₀ to B_y.

*Bounding z.* Take A = "the plaquette p of B₀ is bad", for any of the 24 plaquettes p of B₀; all are faces of B₀.
- The disseminated event ⋂_y θ_y A says that every image of p is bad. A plaquette lies in exactly four blocks (its two
  transverse coordinates each have two choices), so there are at least N/4 distinct bad images.
- Each bad plaquette has s_p > ρ²/4, so the unnormalised weight of the event is at most exp(−β_W ρ² N/16), using
  total Haar mass 1 and S ≥ 0.
- For the normalisation, restrict every link to angle ≤ t. Angles are subadditive, since the angle is the geodesic
  distance to I on S³ and the metric is bi-invariant. So every plaquette has angle ≤ 4t and s_p = 1 − cos(angle) ≤ 8t².
  With 4N links and 6N plaquettes, Z ≥ h(t)^{4N} e^{−48 N β_W t²}.
- At t² = 1/(8β_W) this is Z ≥ h(t)^{4N} e^{−6N}. Hence z(A) ≤ z.

*Proof of (i).* "x bad" is the union over the six plaquettes with base x. Expand ⋂_{x∈S} ⋃_{p∈P_x} {p bad} into at most
6^{|S|} intersections that pick one plaquette p_x for each x. The plaquette p_x lies in the block B_x, so each
intersection has the form ⋂_{x∈S} θ_x A_x and probability at most z^{|S|} by the chessboard estimate.

*Proof of (ii).* A cluster of size ≥ n containing x₀ contains a 𝒢-connected set of exactly n sites through x₀, all
bad. There are at most (eD)ⁿ such sets (the lattice-animal bound; e.g. Grimmett, *Percolation*, §4.2). Each is all bad
with probability ≤ (6z)ⁿ by (i). Sum over n′ ≥ n.

*Proof of (iii).* h(t) ≥ (2/3π) t³ (1 − t²/5), and t³ = (8β_W)^{−3/2}. ∎

## 2. Constants (`rare.py`, `rare.out`)

Smallest β_W at which q < 1, i.e. at which the bad clusters have exponential tails:

| bad threshold ρ | nearest-neighbour clusters (D = 8) | ∗-clusters (D = 80) |
|---|---|---|
| β_W^{−1/4} (Corollary `anyexterior`) | 3.7·10⁶ | 3.9·10⁶ |
| 0.1 | 1.6·10⁵ | 1.7·10⁵ |
| 0.01 | 2.1·10⁷ | 2.1·10⁷ |
| 3.25·10⁻¹² (Corollary `localcalm`, m = 1) | 6.0·10²⁶ | 6.1·10²⁶ |

Above the threshold the tails improve very fast. At ρ = β_W^{−1/4} and D = 80, q = 10⁻³⁰ at β_W = 10⁷ and 10⁻²¹⁰ at
β_W = 10⁸.

**Along the tuning** β_W ∼ b log(1/a) with ρ = β_W^{−1/4}, log(1/q) ≈ β_W^{1/2}/16. So on a cube of physical side r,
which has (r/a)⁴ sites, the largest ∗-cluster of bad sites has at most

    n_max ≈ 64 log(r/a)/β_W^{1/2} ≈ 64 √(log(1/a)/b)

sites, with probability → 1. The bad regions are therefore isolated islands of size O(√log(1/a)) in lattice units:
the large-field regions are small and separated, which is the geometry Bałaban's analysis [Bal89b] assumes. Compare
Proposition `calmrare`(ii): the fibre estimates of Sec. `smallfibres` hold on windows of |W| blocks once
β_W > 10²⁴|W|², i.e. for |W| ≲ 10⁻¹² β_W^{1/2}. The bad islands and the windows the fibre estimates can handle grow at
the same rate, β_W^{1/2}, but the constants differ by about 10¹³. So with the present constants the islands cannot be
absorbed into calm windows.

## 3. Why disagreement percolation does not give "harmless"

Van den Berg–Maes disagreement percolation [vdBM94] bounds the influence of a boundary condition by a percolation
of disagreement. It needs a site-wise coupling: for any two configurations of the neighbours, the two conditional laws
of a site, or window, must agree with probability ≥ 1 − p, with p below the percolation threshold.

Here the conditional law of one link is the von Mises–Fisher law ∝ exp(β_W⟨u, h⟩) on S³. The staple sum h has |h| ≈ 6
in a calm region, so the concentration is κ ≈ 6β_W. Two neighbour configurations whose staples differ by an angle φ
give conditional laws at total-variation distance ≈ c·φ√κ.
- Typical fluctuations in a good region are φ ∼ β_W^{−1/2}, so the distance is of order one, even with no bad
  plaquette anywhere.
- At the bad threshold φ ∼ ρ = β_W^{−1/4}, the distance is β_W^{1/4} ≫ 1, i.e. close to 1.
- The same holds for windows and for block averages: a shift of the boundary of order the fluctuation scale moves the
  window law by an amount of order one in total variation.

So the per-site disagreement probability is not small anywhere at weak coupling. Disagreement percolation applies
only at strong coupling, where it reproduces Dobrushin-type uniqueness and gives nothing new beyond the blocks of
track A.

## 4. What "harmless" would need

The right replacement is the Bassalygo–Dobrushin scheme for random bad sites [BD86], in a metric where good regions
contract:

- (G) On good regions there is a contraction of influences with exponential decay. The candidate is the chart decay
  of Theorem `chartdecay`, rate α_G = 2.5·10⁻⁴ per unit distance, but it is proved only on a fixed fibre and needs
  the passage of Remark `chartdecayopen`(b).
- (B) Bad islands amplify influence by at most a constant C^{|island|}, where |island| counts the sites of the island.
- (P) Islands have exponential tails, with rate exceeding log C plus the entropy. This is the lemma above.

(P) is now proved with explicit constants on the torus. (G) and (B) are the open content of (D). The lemma also does
not cover a box with an arbitrary exterior, which is what (D)'s supremum over boundary conditions asks for. Chessboard
estimates need the torus, and under a general exterior only the density bound of Corollary `anyexterior` is available.
Extending (i) to boxes near the boundary would need a local version of the energy–entropy count. That is plausible
because the count in Proposition `anyexterior` is already local up to the boundary layer, but it is not done here.

**Status of C2:** rarity is proved with explicit thresholds, on the torus. Harmlessness by disagreement percolation is
ruled out at weak coupling, and the remaining route, (G) + (B) + (P), is the open part of (D).

## References

- [Bal89b] T. Bałaban, *Large field renormalization. II. Localization, exponentiation, and bounds for the R operation*,
  Comm. Math. Phys. 122 (1989) 355–392.
- [BD86] L. A. Bassalygo, R. L. Dobrushin, *Uniqueness of a Gibbs field with a random potential — an elementary
  approach*, Theory Probab. Appl. 31 (1986) 572–589.
- [Bis09] M. Biskup, *Reflection positivity and phase transitions in lattice spin models*, in: Methods of Contemporary
  Mathematical Statistical Physics, Lecture Notes in Math. 1970 (2009) 1–86, arXiv:math-ph/0610025.
- [FILS78] J. Fröhlich, R. Israel, E. H. Lieb, B. Simon, *Phase transitions and reflection positivity. I*,
  Comm. Math. Phys. 62 (1978) 1–34.
- [OSei78] K. Osterwalder, E. Seiler, *Gauge field theories on a lattice*, Ann. Phys. 110 (1978) 440–471.
- [vdBM94] J. van den Berg, C. Maes, *Disagreement percolation in the study of Markov fields*, Ann. Probab. 22 (1994)
  749–763.
