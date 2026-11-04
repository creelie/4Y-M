# Shell transport: averaged Conjecture B from the large-field displacement H

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

**Missing standard fact:** a Gauss–Green/divergence theorem on the stratified, piecewise real-analytic
fibre with Lipschitz Y. Expected routine, not written.

## 3. Consequence and what stays open

- **Sketch:** Conjecture B in averaged form, |E Ξ| ≤ C(1+λ)·μ(m ≥ ρ_G/8), follows from
  H = O(1): the displacement of the reweighted law on {m ≥ ρ_G/2}. This replaces CALM_LINK Step 1's
  unproved "surface-only" assumption and answers the audit's bulk-term objection.
- **Open:** H itself, which is the large-field input (Conjecture L / (R1), D7). The uniform
  abelian-boundary example of `CALM_LINK.md` §4 still forces (R1) to be averaged over exteriors.
- D4, D5 (global), D7, D8 and E remain open. This note moves the obstruction for B onto H; it does
  not remove it.

## 4. Next steps

1. Compute the constant in ω and write the bulk bound with explicit constants.
2. Write the stratified Gauss–Green lemma (or cite a version for semialgebraic/subanalytic stratifications).
3. Attack H: test numerically on small SU(2) blocks whether the reweighted displacement on
   {m ≥ ρ_G/2} stays O(1) uniformly in exterior data (CPU job, multi-seed, within the 96-CPU budget).
