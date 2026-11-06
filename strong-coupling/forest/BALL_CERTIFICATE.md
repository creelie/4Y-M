# Ball certificate: interface statement (direction #1 of CLAY_NEXT_DIRECTIONS.md)

Status: **draft target statement, not a theorem.** Nothing here is proved. It fixes
what a finite-box certificate must say so that it can be plugged into a terminal-scale
effective action produced by an external existence proof (e.g. a Bałaban-type
construction), avoiding the |V| factor that killed Lemma S′ (POLYMER_ARCHITECTURE §11).

## Setting

- Block lattice Λ (periodic or Z⁴), link variables U ∈ SU(2)^{E(Λ)}.
- Reference action Φ₀ = Wilson at coupling β₀ plus finitely many local terms
  (rectangles, chairs) with coefficients c ∈ K, K a compact box in R^m.
- Remainder class 𝓡(δ, κ): gauge-invariant interactions R = Σ_X R_X with
  ‖R‖_κ := sup_e Σ_{X∋e} e^{κ|X|} ‖R_X‖_{C¹} ≤ δ.
  (Norm chosen to match a polymer/Bałaban output: exponential tree decay plus one
  derivative, since couplings are controlled by Lipschitz influence.)
- Family 𝓕 = {Φ₀(c) + R : c ∈ K, R ∈ 𝓡(δ, κ)}.

## Certificate (what a computation must deliver)

Box B of side L, boundary ∂B. For Φ ∈ 𝓕 and boundary configurations ω, ω′ differing
on one link j ∈ ∂B, let μ^Φ_{B,ω} be the conditional measure. The certificate is

  (C) sup_{c∈K} sup_{ω} Σ_{j∈∂B} W₁(μ^{Φ₀(c)}_{B,ω}, μ^{Φ₀(c)}_{B,ω^j})|_{interior} / d(ω_j, ω′_j) ≤ ε₀ < 1,

computed with interval arithmetic over K (dump10 T-1..T-3 engine), restricted to an
interior region at distance ≥ r from ∂B, plus an explicit Lipschitz modulus

  (Lip) |∂_R (influence)| ≤ A(L) ‖R‖_κ,

where A(L) is polynomial in L (from a C¹ perturbation of the conditional measure,
one boundary term per plaquette of B: A(L) ≲ β-dependent const · L^4).

## Target theorem (to prove)

If ε₀ + A(L) δ < 1 then every Φ ∈ 𝓕 satisfies a Dobrushin–Shlosman-type block
condition GDS(B, ε₀ + A(L)δ), hence uniqueness and exponential clustering with rate
≥ −log(ε₀ + A(L)δ)/L in block units.

Key point: the hypothesis is uniform over the family, and A(L)δ involves only the
box volume L⁴, not |Λ|. That is exactly what S′ lacked.

## Weak-form variant (pairs with direction #4)

Replace sup_ω by the reflection-positive periodic state and W₁ by covariance response:
κ(r) = β · Σ_{p at distance r} |Cov(O, P_p)|. This is what
`research/experiments/linresp/` measures numerically (heuristic only). A weak-form
target would need a chessboard estimate transporting κ(r) < 1 to decay. Open.

## Open points

1. Is (Lip) with polynomial A(L) true for sup-over-ω influences? The conditional
   measure depends on R through e^{−R}; the C¹ derivative is a covariance in μ_B,
   bounded by L⁴ ‖R‖ crudely. Need to check the large-field part (|R_X| unbounded?).
   A Bałaban output gives R only on small-field configurations; large fields need a
   separate rarity bound (cf. §10 chessboard density).
2. Which (κ, δ) does a Bałaban terminal-scale action come with? To read off the
   OpenAI text when released.
3. Feasible L: decided by the linresp feasibility run (direction #3).

## Benchmark frozen (2026-11-06)

dump10's 0.348/0.52 vs 0.351/0.77 discrepancy is a version lag, not an error:
commit 1039165 (2026-10-31) took the SU(2) plaquette blocks from 0.345 to 0.351, and
1f7acd4 (2026-11-02) extended to 0.351 / 0.77 and updated the README; dump10 quoted
an older README. Benchmark to beat: SU(2) β_W ≤ 0.351 (α(0.351) ≤ 0.99995, Lean
kernel, fails at 0.352), SU(3) β_W ≤ 0.77.
