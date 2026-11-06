# Discharging `hW` for the SU(2) heat-bath kernel at β_W < 2/9

Labels: **proved**, **sketch**, **open**. The Lean is **not** written; this is the proof plan, with the missing formal piece named.

## Target
`finite_dobrushin_comparison` / `compact_comparison_of_discretizations` (`dump2.txt` l. 3298, 4478) assume `hW`.
For every n, every i ∉ S and every x, y, there must be a finite coupling π of k n i x and k n i y with

  tcost(dist∘e n) π ≤ Σ_j C_ij · dist(e n x_j, e n y_j) + η_n,  η_n → 0,  sup_i Σ_j C_ij < 1.

`finite_dobrushin_of_approx` reduces this to `hWref` for reference kernels (error η) plus `happrox` (the true kernels are
within δ_n of the references), giving total error η + 2δ.

## Proof
1. **Continuum single-link bound (proved, `prd/main.tex` Lemma `rate` l. 822, Thm `quarter` l. 878).** For the SU(2)
   heat bath ν_b ∝ exp(Re tr(b u)) on S³: W₁(ν_b, ν_b′) ≤ |b − b′|·sup L, with L ≤ 1/4.
2. **Staple Lipschitz (proved).** Moving one neighbour link moves the staple b by ≤ 4|β| (chordal), so
   C_ij = |β| for each of the 6(d−1) = 18 neighbours, and Σ_j C_ij = 18|β| < 1 ⇔ β_W = 4β < 2/9.
   This matches `wilson_threshold_d4` (dump2 ~9040–9140).
3. **Reference kernel (sketch).** Let k_ref be the push-forward of ν_{b(e∘x)} under the Voronoi nearest-point map onto the net Y_n.
   Then the Voronoi map moves mass by ≤ mesh_n, so W₁(k_ref x, k_ref y) ≤ W₁(ν_b, ν_b′) + 2 mesh_n, giving η_n = 2 mesh_n.
4. **Approximation (sketch).** The true finite Gibbs single-site kernel on Y_n has weights exp(Re tr(b y))·w_n(y),
   with w_n the cell volumes. It differs from k_ref only through evaluating the density at y rather than averaging it over the cell, so
   its TV distance is ≤ e^{2|b|}·Lip·mesh_n → 0 uniformly (|b| ≤ 24β bounded). So δ_n = O(mesh_n), and with diam(SU(2)) finite,
   W₁ ≤ diam·TV.
5. **Coupling, not just W₁ (open formally).** `hW` asks for an explicit coupling π achieving the cost. On a finite set this is
   finite Kantorovich–Rubinstein duality / existence of an optimal transport plan for a finite LP. It is standard, but
   **I did not find it in Mathlib**. It is the one missing formal ingredient: an optimal coupling exists (compactness of
   the transport polytope) and its cost equals sup over 1-Lipschitz f. Only the "inf attained" half is needed if
   step 1 is itself proved via a coupling (the paper's proof is by a synchronous coupling along the rate ODE,
   which can be discretised directly, avoiding duality).

## Verdict
Mathematically `hW` holds at β_W < 2/9 with C_ij = |β| on the 18 neighbours and η_n = O(mesh_n): **sketch, no gap I
could find**. Formalising it needs (i) the quarter theorem in Lean (not present: `ZonalBound.lean` has
`zonal_final` under 9k² ≤ 88, a related but separate bound), (ii) finite OT existence, and (iii) the Voronoi/mesh
estimates. The chain becomes unconditional in that range only after those are written and checked.

## Update (Nov 2026): what the formal gap in step 5 actually is
- Existence of a minimal coupling is **already in our Lean** (`exists_min_coupling`, used by `finite_dobrushin_comparison`,
  dump2 l. ~3310, by compactness). So "inf attained" is not missing.
- What is missing is **min-cost ≤ dual bound**: finite Kantorovich–Rubinstein, i.e. the minimal coupling cost is at most
  sup over 1-Lipschitz f of ∫f d(k x − k y). The continuum Lemma `rate` is stated in that dual form.
- Mathlib (cluster checkout) has no optimal transport, but it does have the Farkas lemma / cone duality
  (`Mathlib/Analysis/Convex/Cone/Dual.lean`, `InnerDual.lean`). Finite KR is then the LP duality of the transport
  polytope, derived from Farkas. That is standard and a few hundred lines.
- Value check: this raises a *formal strong-coupling* range (β_W < 2/9). It does not touch the weak-coupling wall
  (large-field expansion) that blocks the Clay statement. Priority is low; it would be worth doing as a formalisation
  deliverable, not as Clay progress.
