# The large-field wall of (D): what can be moved (November 2026)

This note is a research pass on the two blockers of (D) named in `CALM_LINK.md` §4 and `ANNEALED.md` §5:
- **(R1)/(I1):** one-point calm rarity under every exterior, and quenched coefficients for the laws μ^ω conditioned
  on the large-field set;
- **D8:** composition across scales, with its crossover base case at β_top ≈ 10–700.

Every item is labelled **proved**, **sketch**, **heuristic** or **computation**. Nothing here closes (D).

The main outcome is §3. It changes what is conditioned on: a **fixed product of link balls in fibre gauge**, instead of
plaquette-level calm conditions. This removes the thin-lens problem on both sides of the mixture of `ANNEALED.md` §5,
and it moves the rarity threshold from K/(8n_P) ≈ 10⁻⁷⁹ to a fixed radius. In the island ledger the threshold drops
from 8.3·10¹⁶³ to 4.7·10⁴¹ at basin radius K, and lower for a larger basin (§4). It needs rarity only on the torus and
averaged over the constraint v, which avoids (R1) at a single scale. It does nothing for D8 beyond about 19 scales (§4).

## 1. Where the huge thresholds come from (computation)

The thresholds β ≥ 10⁶⁴ (`CALM_LINK.md` §3) and β ≥ 10¹⁶⁴ (`ISLANDS.md` §3) have two sources.

1. **L² calm conditions.** Proposition `calmfibre` needs every plaquette ≤ K/(4√n_P), and `calmrare`(ii) needs
   ≤ K/(8n_P). Both are L² (sum-over-plaquettes) conditions turned into pointwise ones, so the plaquette threshold
   shrinks with the number n_P ≤ 576m of plaquettes in the window. With m = 2.9·10²⁴, ρ_c ≈ 10⁻⁷⁹, and rarity at ρ_c
   costs β ≳ ρ_c⁻².
2. **Entropy factors with exponent ∝ m.** ε_loc = (3·10⁸ n_P^{3/2}/K³)^{73m} e^{−3λK²/4}. The prefactor is a volume
   ratio of the window's fibre, so its log is linear in m, and it must be beaten by λK².

Both come from treating the window as one chart. A condition that is **L^∞ in the link variables** (each link in a
fixed ball) has a threshold that does not depend on m. That is the motivation for §3.

## 2. A box-local Peierls lemma for excess action (proved), and why it is marginal (heuristic)

**Lemma (excess action).** Let Λ be a box of N links, with any exterior configuration. Let S_Λ be the Wilson action of
the plaquettes that meet Λ, and S_min the minimum of S_Λ over the links of Λ (exterior fixed). Let A be the event
S_Λ ≥ S_min + E. Then

  μ(A | exterior) ≤ η(N, E)⁻¹ e^{−βE/2},   η(N, E) = Haar(geodesic ball of radius E/(2C_N))^N,

where C_N is the Lipschitz constant of S_Λ in the product geodesic metric (C_N ≤ 6·const·N, box-only).

*Proof.* The conditional density is Z⁻¹e^{−βS_Λ}. Then μ(A) ≤ Z⁻¹e^{−β(S_min+E)}. The product ball of radius
r = E/(2C_N) around a minimiser has S_Λ ≤ S_min + E/2, so Z ≥ η e^{−β(S_min+E/2)}. Divide. ∎

The constants depend on N and E only, and the exterior is arbitrary. So **excess** action is rare under every
exterior, with no chessboard and no torus.

**Why it does not give (R1).** Write the curvature in Λ as forced + excess: the forced part is the curvature of the
minimiser, which is fixed by the exterior. (R1) needs the *absolute* curvature to be small. The forced part is, at the
linearised level, a harmonic (Poisson) average of the boundary data, and its sup over the interior is comparable to
the boundary sup: **attenuation 1**. A constant abelian-embedded field (`ROADMAP.md`, lead (b)) realises this. So a
recursion "calm boundary ⇒ calm interior with a gain" gets no gain from the lemma, and the scheme is marginal.
(Heuristic: at the nonlinear level there may be a nonabelian gain, but no mechanism is identified.)

## 3. Fixed link-ball product conditioning in fibre gauge (sketch; the main idea)

**The change.** Work in the fibre gauge of the block constraint (the coordinates of `calmfibre`). Fix a radius r₁ and,
for each link ℓ in that gauge, the geodesic ball B_ℓ = B(1, r₁). The **good set** of a configuration is
G = {ℓ : U_ℓ ∈ B_ℓ}, and the bad set is B = G^c. Mix over ω = B as in `ANNEALED.md` §5.

**Bad side: no thin lens (proved, trivial).** Conditioned on ω, a bad link is restricted to the complement of a
**fixed** ball, which does not move with its neighbours. Its conditional law is e^{−βS} times the indicator of a
fixed set C. For a 1-Lipschitz f, the derivative of E[f] in a neighbour σ is −Cov(f, β∂_σS), because C does not depend
on σ. So it is at most β sup|∂_σS| · diam G, uniformly in the field and in the mass of C. Hence (B) of `ANNEALED.md` §1
holds on ω with Ā = O(β_W). Compare `ANNEALED.md` §4: there, the allowed sets were intersections of
balls centred at the neighbours, hence thin lenses with unbounded response.

**Good side: a convex domain (sketch).** Conditioned on ω, the good links live in the product Π_{ℓ∈G} B_ℓ, which is
geodesically convex for r₁ below the injectivity radius. If the Hessian of βS restricted to the fibre is ≥ βc₀ on
this product (input (a) below), then:
- **Bakry–Émery:** the conditional law satisfies an LSI with constant ≍ (βc₀)⁻¹, uniformly in the shape of G (convex
  restriction does not hurt Bakry–Émery);
- **Herbst:** sub-Gaussian concentration of every Lipschitz function, so the rarity of "a plaquette > r₁/4" inside G
  is e^{−cβr₁²} under the conditioned law, under any configuration of ω;
- **Helffer–Sjöstrand / Brascamp–Lieb:** Cov ≤ ⟨∇F, (Hess)⁻¹∇H⟩, and Combes–Thomas gives exponential decay of the
  inverse Hessian at rate ≍ √c₀ in lattice units, because the Hessian is a local operator;
- **no off-chart term:** the domain is the chart, so the error 320λ|V₁|ε_loc of `ISLANDS.md` §1 does not appear.

So the quenched input (I1) for μ^ω reduces to input (a), a deterministic Hessian bound.

**Rarity of the bad set (sketch).** A bad link forces, within its block closure, a plaquette with δ > r₁/4
(input (b)). The torus chessboard of `rare/RARE.md` then gives P(S ⊂ B) ≤ q(r₁/4)^{|S|/c}, at a **fixed** threshold.
Because the mixture is averaged over v and over ω, this is used only annealed and only on the torus. The one-point
rarity under every exterior, (R1), is **not** needed at this scale.

**Open inputs.**
- **(a) Hessian under L^∞ smallness.** Is Hess(βS|fibre) ≥ βc₀ whenever every link is within r₁ of 1 in fibre gauge?
  The paper's chartdecay/fibrechart bound is stated under calm (plaquette) hypotheses. The Hessian of S at U is
  the lattice Laplacian-type form plus O(max_p δ_p) corrections, and max_p δ_p ≤ 4r₁ on the product. So the bound
  should hold for r₁ ≤ c with c an absolute constant, *provided* the fibre gauge kills the zero modes uniformly. That
  last point is the check to do in `main.tex`.
- **(b) Bad link ⇒ large plaquette.** In fibre gauge a link far from 1 should force curvature nearby, since the gauge
  is fixed by the block averages. This is a gauge-fixing estimate (a Poincaré-type inequality for the fibre
  coordinates in a block); it is plausible and unchecked.
- **(c) The Cov_ω term.** The environment covariance of `ANNEALED.md` §5 is still closed by the bootstrap Lemma 3, with
  ε ≍ q·|N|²/r₁² now, which is smaller than before.
- **(d) D8.** See §4.

## 4. Ledger (computation; `d5/basin.py`, `d5/basin.out`)

SU(2), C_G = 3·10¹³, island cost (4.5β C_G)^{324} per block (`ISLANDS.md` §2), torus rarity
log q(ρ) = −βρ²/16 + 6 log β + 24.7 + log(6e·80). Condition: q · cost < ½, with ρ = r₁/4.

| basin radius r₁ | ρ = r₁/4 | β needed |
|---|---|---|
| 4.8·10⁻¹⁸ (= K) | 1.2·10⁻¹⁸ | 4.7·10⁴¹ |
| 10⁻⁸ | 2.5·10⁻⁹ | 7.2·10²² |
| 10⁻⁴ | 2.5·10⁻⁵ | 5.6·10¹⁴ |
| 10⁻² | 2.5·10⁻³ | 4.8·10¹⁰ |
| 10⁻¹ | 2.5·10⁻² | 4.4·10⁸ |

Compare 8.3·10¹⁶³ (`ISLANDS.md` §3) and 1.1·10⁶⁴ (`CALM_LINK.md` §3). The island cost still carries C_G from the
calm-data route. If input (a) holds, C_G is replaced by the Combes–Thomas constant and the numbers fall further.

**Coarse scales by chessboard alone.** A coarse block of side 2^k with holonomy defect r₁ forces a fine plaquette
> r₁/4^{k+1} (Stokes over ≤ 4^{k+1} plaquettes). Counting the scales k at which the same island condition holds:

| r₁ | β | scales usable |
|---|---|---|
| 4.8·10⁻¹⁸ | 10³⁸ | 0 |
| 4.8·10⁻¹⁸ | 10⁶⁴ | 19 |
| 10⁻² | 10⁸ | 0 |
| 10⁻² | 10²⁰ | 8 |

The fine-scale threshold r₁/4^{k+1} decays geometrically, so chessboard rarity runs out after O(log β) scales. Along
β ∼ b log(1/a) the number of scales needed is ∼ log(1/a), which is exponentially more. **D8 needs RG-improved rarity:
coarse curvature must be rare at the coarse coupling, which is Bałaban's large-field estimate.**

## 5. Negatives recorded

- **Concentration cannot give (R1) cheaply (proved, by citation).** Stroock–Zegarlinski: a uniform LSI for all
  finite-volume conditional laws is equivalent to complete analyticity (Dobrushin–Shlosman). So deriving (R1) from
  concentration under every exterior is as hard as D5 itself. §3 escapes this only because the conditioning set is
  convex, and that is input (a).
- **Partial reflection positivity (heuristic).** A chessboard estimate with Dirichlet boundary on part of the box has
  constants that depend on the volume through the partition-function ratio. It does not give exterior-uniform rarity.
- **No contour structure (proved, elementary).** The Bianchi identity only gives δ_p ≤ Σ of the other five faces of a
  cube. Large plaquettes do not form closed contours, so Peierls contour counting does not apply directly.
- **Dynamics (heuristic).** Path coupling (Hayes–Vigoda type) on Glauber dynamics gives mixing bounds whose constants
  grow with the volume at weak coupling; no volume-uniform statement comes out.
- **Chessboard is not RG-improved (computation, §4).**

## 6. What remains, in order

1. **Input (a)** in `main.tex`: does the fibre-gauge Hessian bound use only L^∞ smallness of the links, with the zero
   modes killed uniformly? If yes, (I1) at one scale is a sketch with a deterministic input, and the islands threshold
   is ≤ 4.7·10⁴¹.
2. **Input (b):** the Poincaré-type estimate "bad link ⇒ large plaquette" in fibre gauge.
3. **D8:** RG-improved coarse rarity. Nothing here helps; this remains Bałaban's large-field expansion.

Net effect on the blockers: (R1)/(I1) at a single scale moves from "a piece of Bałaban's expansion" to "two
deterministic checks, (a) and (b)" (sketch). D8 is unchanged.
