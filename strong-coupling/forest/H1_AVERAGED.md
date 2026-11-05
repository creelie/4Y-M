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

**Finding.** The dimension-count argument in my earlier sketch (negative modes once the flux quanta exceed the constraint count) is **wrong**. The block map kills the lowest-Landau-level modes until the magnetic length is below the block scale. For every tested map, the fibre is stable whenever the coarse plaquette angle b²φ ≤ π. Instability appears only above π (sharp at L = 32, see below), where the coarse field is already large (beyond the angle π, where it aliases), and not a calm datum.

**L = 32 rerun (job 1662119).** The threshold is sharp. For b = 4 and b = 8, with both block maps, the fibre has **no** negative
directions for coarse angle b²φ ≤ π, and has negatives as soon as b²φ > π (first tested points: 1.12π for b = 4, 1.25π for b = 8).
The L = 16 onsets "3π/2" and "6π" were artefacts of a coarse flux grid and a 2×2 coarse lattice.
*Reading:* past π the principal-log coarse datum aliases to angle b²φ − 2π. The uniform field is then no longer the
low-action point of its own fibre. So the instability is exactly an aliasing effect, not a Nielsen–Olesen
effect inside calm data.

## 3. Verdict
- The uniform-field obstruction does **not** show up at the linearised level for calm coarse data (L = 16, b ≤ 8). It is therefore **not**
  evidence that sup-over-V (H1) fails on the calm class. `CALM_LINK.md` §4's c̄ ≈ 1 attenuation heuristic (no decay transmitted through a
  classical uniform background) is a separate, non-perturbative concern and remains unrefuted. Nothing here proves (H1) in either form.
- Sup-over-V **over all V** is still unreasonable, because large-angle data do have fibre saddles. The sound targets are (H1) on calm V
  (angle < π per coarse plaquette, plus a margin), or (H1-avg). Large V should be handled by the polymer/large-field expansion.
- Limits: linear stability only, one slice size, straight-path maps (not the [CD] V♯ map), floating point.
  The threshold is now located (π, sharp at L = 32). Next check: repeat with the [CD] V♯ averaging map.

## Mechanism of the θ = π threshold: Nielsen–Olesen lowest Landau level (heuristic, checked against existing data, 2026-11-05)
Background diag(e^{iφx}, e^{−iφx}). The adjoint W± components carry charge 2, so they see flux 2φ per plaquette and total flux 2·2πn.
Predictions from the continuum Nielsen–Olesen picture:
- **Unconstrained:** the lowest Landau level has 2n complex modes, i.e. 4n real negative directions, with eigenvalue ≈ −φ/2.
  Data: free_neg = 4n exactly for every n in both runs (L=16 and L=32), and free_min = −0.0248 at φ = 0.0491. ✓
- **Constrained (block b):** an LLL mode can avoid the block constraint only if it is localised inside one block. That requires
  each block to carry more than one charged flux quantum: 2·b²φ > 2π ⇔ θ = b²φ > π. Just above threshold the predicted count
  is one complex mode per block, i.e. 2·(L/b)² real.
  Data: 32 = 2·16 (L=32, b=8, θ=1.25π); 128 = 2·64 (L=32, b=4, θ=1.12π); 32 = 2·16 (L=16, b=4, θ=1.5π, 'corner');
  8 = 2·4 (L=16, b=8, θ=2π, 'corner'). ✓ At θ = π exactly there are no negatives (b2 n=32, b4 n=8 and n=32, b8 n=2). ✓
- **Crossing:** the min eigenvalue at L=32, b=4, θ = 1.12, 1.25, 1.37, 1.5 π is −.034, −.062, −.088, −.115. That is linear in θ − π
  and extrapolates to ≈ −0.006 at π. So the crossing is consistent with continuous closing, not a jump. Job uniform-field-gapscan
  measures the θ < π side, with the prediction gap ∝ (π − θ)/b².

**Consequence.** The instability is the standard Nielsen–Olesen one. The block constraint removes it exactly while the coarse
plaquette angle is below π, i.e. while the coarse plaquette stays away from −I. Bałaban's small-field condition (coarse plaquette
near I) already imposes this margin. So uniform fields are **not** an obstruction to idea D. They do show that the small-field
convexity constant must degrade like (π − θ), and cannot be uniform in V. Still heuristic: these are floating-point numbers,
and the LLL count is matched, not proved.
*Adversarial caveat.* The 'all' map (b parallel lines per block) is stronger. At L=16, b=8 (only 2×2 blocks) it shows no negatives
at θ = 2π, and negatives only appear by θ = 6π. At L=32 it does follow the θ > π threshold (b=8: 32 negatives at 1.25π), but with
smaller |min|. So the clean threshold is established for 'corner'. For 'all' it holds at L=32, while the L=16, b=8 case looks
like a finite-size effect (too few blocks for a block-localised LLL mode to fit). This is not settled.
