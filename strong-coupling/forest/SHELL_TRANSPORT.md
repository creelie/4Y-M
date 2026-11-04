# Shell transport: Conjecture B for local observables, and what remains

Status labels: **proved** = checked line by line here; **sketch** = argument complete modulo a stated
standard fact; **open** = not proved. Fronts D and E are **not** closed by anything in this note.

## 0. Review of `AVERAGED_SHELL_MASS.md` and `SHELL_REMAINDER_AUDIT.md`

- SU(2) hazard bound f ≤ e·sec²(1)(κ+1)·P(R ≥ r), its extension to a maximum of dependent radii, and
  the averaging/Fatou step (1)–(2): **proved** (rechecked). The chordal change of variables
  δ = 2√2 sin(R/2) has dδ/dR ≥ 1 for R ≤ π/2, as used. The speed constant 4 is valid; Cauchy–Schwarz
  gives 2 for distinct links (cosmetic).
- The audit's point is **correct**: Ξ differentiates the whole large-path contribution, so a bound on
  the shell *surface* mass alone does not control Ξ; bulk terms must be handled. Its obstruction example
  is right.

## 1. Global transport (proved)

Lemma `privatesubmersion` (`main.tex` l. 8711) holds at **every** differentiability point of V♯, not only
on R = {m < ρ_G/2}. So Y = −J*(JJ*)⁻¹ ∂_t c is defined a.e. on the whole fibre, with |Y| ≤ |∂_t c|/σ₀.

- m > ρ_G: V♯ = ℓ₀ is independent of inner links, ∂_t c = 0, **Y = 0**.
- ρ_G/2 < m < ρ_G: ∂_t c = E(χN)[χ δN + χ′ N], |χ δN| ≤ Λ/16, |χ′| ≤ 2/ρ_G, |ṁ| ≤ 1, |N| ≤ (π/3)ρ_G,
  so |∂_t c| ≤ Λ/16 + 2π/3 and **|Y|₂ ≤ (Λ/16 + 2π/3)/σ₀ ≈ 69.6**.
- m < ρ_G/2: the [CD] transport, unchanged.

Numerical sanity check (heuristic, `d5/shell_transport_check.py`, cluster job 1657627, 60 000 SU(2)
configurations in calm / shell / large regimes): min σ(J) = 0.353 (bound 0.031), max |Y| = 0.84
(bound 69.6), Y ≡ 0 for m > 1. The proved constants are therefore very coarse but not violated.

## 2. Decomposition of Ξ (sketch)

Differentiating E_t F with the global transport (Gauss–Green on the fibre, piecewise smooth):

1. transport term E⟨∇F, Y⟩;
2. bulk terms Cov(F, ℓ̇), ℓ̇ = −λ DQ[Y] + ω, with |DQ[Y]| ≤ 24√r|Y|₁ a.e. on the shell
   (constant in ω not yet computed);
3. singular divergence terms from jumps of the normal component of Y on {m = ρ_G/2}, {m = ρ_G} and on the
   tie sets {δ_a = δ_b} inside the shell.

The singular mass is ≤ 2 sup|Y| × the level-surface mass of `AVERAGED_SHELL_MASS.md` (2a)/(5), i.e.
O(κ+1)·μ(m ≥ ρ_G/2). Tie surfaces fall under the same hazard lemma: condition on everything except the
private link, with measurable centre and level.

**Gauss–Green: now written** in `GAUSS_GREEN.md` (Lemma GG, ambient averaged form, proved modulo Federer).
That note also fixes the localisation (ψ₄ = (Σδ⁸)^{1/4}, giving supp(1−θ) ⊂ {m ≥ ρ/8}). The tie-set
hazard step remains a sketch.

## 3. Consequence and what stays open

**Bulk terms, bounded local F (sketch).** On every regular fibre the bulk score is
λ·(D Q[Y] + ∂_j Q) plus the derivative of the coarea density. Since SU(2) is compact, |∂Q| ≤ C per moved link,
|Y| ≤ 69.6, and |∂ log det JJ*| ≤ C/σ₀² on smooth branches. Everything is supported where 1 − θ ≠ 0, which lies
inside {m ≥ ρ_G/8}. So for |F| ≤ 1,

  |Ξ_j(F)| ≤ C(1 + λ)·π(m_y ≥ ρ_G/8),

which is Conjecture B as stated in [CD], for bounded F supported near j. This is exactly the "O(β) pointwise
bulk score × shell mass" route that `SHELL_REMAINDER_AUDIT.md` names. It does not need H, and it answers the audit's
bulk objection for local F. On m > ρ_G the whole term is −λ Cov(F(1−θ), ∂_j Q), which is bounded the same way.

**What the cube coefficient actually needs: decay (open).** For F supported at distance d from j, the bound must be
C(1+λ)e^{−αd}π(·), not uniform in d. The bulk term is λ Cov(F, ∂_jQ·1_{m≥ρ/2}), a covariance between a far
observable and a local large-field event. Its decay in d is **Conjecture L** (localised off-chart covariance),
which is the quantity called H in the earlier version of this note. Under an arbitrary exterior it needs (R1), and
`CALM_LINK.md` §4 shows that (R1) pointwise in the exterior is false in general (uniform-field boundaries). So:

- Conjecture B, bounded Lipschitz local F, φ-averaged: **sketch** (`GAUSS_GREEN.md` §4; only the tie-set step and C_G remain).
- Conjecture B with separation decay ⇔ Conjecture L ⇐ averaged (R1) + averaged cube coefficient: **open**.
- D4, D5 (global), D7, D8 and E: **open**. The wall is unchanged. It is the large-field part of a multiscale
  expansion (Bałaban/Dimock), where large-field regions are summed in a polymer expansion over the torus measure
  rather than bounded under an adversarial exterior.

**Dead ends checked (not reopened):** disagreement percolation (`rare/RARE.md`: fails at weak coupling);
finite-intersection rarity alone (`ENVIRONMENT_OBSTRUCTION.md`: hidden-bit counterexample); pointwise (R1)
(uniform-field boundary).

## 4. Next steps

1. Compute the constant in ω and write the bulk bound with explicit constants.
2. Write the stratified Gauss–Green lemma (or cite a version for semialgebraic/subanalytic stratifications).
3. The only live route to decay: replace the supremum over exteriors in the cube coefficient by a torus
   polymer expansion of large-field regions (averaged (R1)). This is research-scale (a one-scale Dimock-style large-field expansion for the
   block variables), not a refinement; numerics at the required β_W ≳ 10⁶⁴ are not informative.
