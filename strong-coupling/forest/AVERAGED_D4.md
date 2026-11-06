# Averaged D4: replacing the supremum over exteriors by a Gibbs average (November 2026)

Labels: **proved**, **sketch**, **open**. Nothing here closes D or E.

## 1. Why try it
Pointwise (R1) is false under uniform-field exteriors (`CALM_LINK.md` §4). Chessboard rarity holds **on the torus,
averaged over exteriors** (λ|L₁|p_nc ≤ 10⁻³ at β_W ≥ 2.1·10⁹⁵, `CALM_LINK.md` §3). Lemma `mixture` and Conjecture B
(`GAUSS_GREEN.md`, now proved for local F) consume only averaged statements. So: can D4 be stated in averaged form?

## 2. The exact decomposition (proved)
For the torus Gibbs measure μ, a cube Λ, its exterior σ-field 𝓔, and F, G with supp F ⊂ Λ:

  Cov_μ(F, G) = E_μ[ Cov(F, G | 𝓔) ] + Cov_μ( E[F|𝓔], E[G|𝓔] ).               (1)

The first term is a quenched covariance averaged over Gibbs exteriors. Averaged (R1) plus the calm cube coefficient
(`D5.md`, `CALM_LINK.md` §§2–3) bound it on the calm-exterior event, which has μ-mass ≥ 1 − 10⁻³/(λ|L₁|). On the
complement it is ≤ ‖F‖∞‖G‖∞ · μ(non-calm), with no decay in d(F, G).

## 3. Where it fails (sketch, sharp)
- **Term 2 is the same problem one scale up.** E[F|𝓔] is a function of the exterior, so its covariance with G is a
  decay statement for the law of 𝓔. That is the Corollary `compose` recursion, i.e. D8.
- **The non-calm complement costs μ(non-calm) ≈ 10⁻³/(λ|L₁|) with no decay.** Summed over the ~d⁴ cubes between F and
  G, this gives no exponential decay. Averaged rarity of a single window is not enough: one needs the large-field set
  to be a **subcritical random set** (exponential tail of connected clusters) under μ. Chessboard gives finite-intersection
  rarity e^{−cβρ²|S|} on the torus, so a Peierls count of clusters **does** give subcriticality at β_W ≥ 2.1·10⁹⁵
  (sketch: the standard lattice-animal count with the chessboard bound per site). With that, term 1 decays like
  e^{−α d} plus the cluster-tail term.
- **Conditioning on a large-field cluster** requires the quenched coefficients (I1) at the cluster's boundary. These
  are bounded (row ≤ 4.5β_W, `D4_D8.md` §1) but not contracting. Disagreement percolation then needs the cluster tail
  to beat the e^{O(β_W)} cost per site, and chessboard rarity e^{−cβ_Wρ²} does (ρ fixed). **This step is the live
  candidate**: it would give single-scale averaged D (term 1 of (1)) at β_W ≳ 10⁹⁵. Constants are not checked. It
  must be reconciled with `rare/RARE.md`, which reports disagreement percolation failing at weak coupling; that note
  uses pointwise rarity under arbitrary exteriors, not torus-averaged chessboard rarity, so it may not apply. **Open.**

## 4. Outcome
- Averaged single-scale D: reduced to one checkable step (§3, third bullet). **Open, next to work on.**
- Term 2 of (1) is D8, whose base case sits at β ≈ 10–700 (`D4_D8.md` §2). **Unchanged; out of reach.**
- E: unaffected.
