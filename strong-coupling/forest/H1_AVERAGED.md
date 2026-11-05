# (H1): averaged restatement, and the uniform-field test of the sup-over-V form

Labels: **proved**, **heuristic** (floating-point numerics), **open**.

## 1. Restatement
dump1's (H1) has the form Σ_{e≠j} sup_{V,X,Y} |Cov_{π_V}(a^X_e, a^Y_j)| ≤ (1 − h_{b−2})/3, with a sup over **all** coarse data V.
The averaged form, which is what Lemma `mixture` / `CALM_LINK.md` §2 actually consume, is:

**(H1-avg).** For the torus Gibbs law ρ of the coarse field,
  Σ_{e≠j} ∫ sup_{X,Y} |Cov_{π_V}(a^X_e, a^Y_j)| ρ(dV) ≤ (1 − h_{b−2})/3,
or, more weakly, the same bound with sup_V restricted to a class 𝒱 with ρ(𝒱ᶜ) summable in a polymer expansion.

(H1-avg) is implied by sup-over-V (H1). The converse is false in general.

## 2. Test: is a uniform abelian field a saddle inside its fibre? (heuristic)
Setup: a 2D L = 16 slice. The restriction is exact, because x₃,x₄-constant fluctuations of the (1,2) links with the 3- and 4-links fixed see only the 2D Hessian.
The background is SU(2), U₂ = exp(iφxσ₃), φ = 2πn/L². The Hessian is the Wilson action Hessian (finite differences), restricted to ker J, where J is the linearised straight-path block map:
'corner' uses one line per coarse link; 'all' uses all b parallel lines. Job 1662049 (`uniform-field/`).

Sanity checks pass. At n = 0 there are no negatives. |∇S| ≈ 10⁻¹², so the background is critical. The free negative count is exactly 4n (Nielsen–Olesen).

| b | map | first n with fibre negatives | coarse angle b²φ at onset | last stable b²φ |
|---|---|---|---|---|
| 2 | all/corner | 48 | 4.71 | 3.14 (n=32) |
| 4 | all/corner | 12 | 4.71 | 3.14 (n=8) |
| 8 | all | 12 | 18.8 | 12.6 (n=8) |
| 8 | corner (weak, codim 24) | 4 | 6.28 | 3.14 |

**Finding.** The dimension-count argument in my earlier sketch (negative modes once the flux quanta exceed the constraint count) is **wrong**. The block map kills the lowest-Landau-level modes until the magnetic length is below the block scale. For every tested map, the fibre is stable whenever the coarse plaquette angle b²φ ≤ π. Instability appears only at b²φ ≥ 3π/2, where the coarse field is already large (beyond the angle π, where it aliases), and not a calm datum.

## 3. Verdict
- The uniform-field obstruction does **not** show up at the linearised level for calm coarse data (L = 16, b ≤ 8). It is therefore **not**
  evidence that sup-over-V (H1) fails on the calm class. `CALM_LINK.md` §4's c̄ ≈ 1 attenuation heuristic (no decay transmitted through a
  classical uniform background) is a separate, non-perturbative concern and remains unrefuted. Nothing here proves (H1) in either form.
- Sup-over-V **over all V** is still unreasonable, because large-angle data do have fibre saddles. The sound targets are (H1) on calm V
  (angle < π per coarse plaquette, plus a margin), or (H1-avg). Large V should be handled by the polymer/large-field expansion.
- Limits: linear stability only, one slice size, straight-path maps (not the [CD] V♯ map), floating point.
  Next check: repeat with the V♯ averaging map and L = 32, and scan b²φ ∈ (π, 3π/2) to locate the threshold.
