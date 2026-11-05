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

## Gap scan result (job 1663501, L=16, heuristic floating point)
**'corner' map:** gap(θ) ≈ gap₀ · (1 − θ/π), with the ratio to that law in [1.01, 1.21] over every θ < π
(b = 2, 4, 8). gap₀ · b² = 2.0, 2.34, 2.44. So **gap ≈ (2–2.4)/b² · (1 − θ/π)**: the gap closes linearly and continuously as θ → π.
At θ = π exactly, the soft modes become exact zero modes (zero count 576 → 704 = 576 + 4n with n = 32, the Landau-level count).
So the "gap" reported at θ = π is the next level up, not a reopening. This confirms the prediction made from the lowest-Landau-level
mechanism before the job was read.

**'all' map: inconclusive.** The zero-mode count depends on n (240 when 8 | n, 208 otherwise, and 504/456/472 for b = 4). At odd n,
32 near-zero eigenvalues (~7e-4) appear even at θ ≈ 0.03π. These are accidental zero modes, flat at n = 0, that commensurability
splits. They are not a convexity gap, and the 1e-6 zero cut cannot tell them apart from soft modes. To measure the 'all' gap,
we need the exact gauge-orbit dimension (to separate gauge modes from soft modes) rather than a numerical cut.

**Consequence for idea D.** On uniform data, small-field fibre convexity holds with constant ≍ (π − θ)/b². A calm class
θ ≤ π − δ gives a gap ≳ δ/b². This is compatible with Bałaban's small-field condition. It does not touch the large-field wall.

## Exact gauge removal (gaugegap v2, request f0c061083e8c714deb103b2d, L=8) — heuristic

The v2 analysis works in ker J coordinates, with the residual gauge computed exactly. Validity: leakage ≤ 5e-10, and the singular-value cut is clean (gap ~1e10).
- 'corner' map: reproduces gapscan exactly (b=2: 0.5 → 0.0745 at θ=0.88π; b=4: 0.146 → 0.090 at 0.5π).
- 'all' map: there are **no zero modes and no negatives below π** after exact removal, so the earlier 1e-6-cut
  "zeros" were gauge. Genuine soft modes remain. Their gap depends on flux parity (b=2: n=0 0.157, n=1 0.0105,
  n=2 0.0358, n=4 0.0872, n=6 0.0141, n=7 0.0022), and the residual-gauge count depends on n (72 vs 56).
- Hypothesis under test: the odd-n soft modes are torus/toron artefacts with gap ~ L^-4.
  Jobs: L=12 (4e528eba147bb9e8d8dc3b1f), L=16 (848a2647c6a2448929299ffe).
- At θ=π exact zeros appear (LLL), as before. Above π, negatives appear, except at b=4 'all', θ=2π, n=4.

### L=12 result (job 1663616) and the 'all'-map soft mode — heuristic, with mechanism
- 'corner' map: the gaps are L-independent (n=0: 0.5 / 0.146 at b=2/4, the same as at L=8). At b=4, n=5 (θ=1.11π)
  there are 18 negatives, as the θ>π rule predicts.
- 'all' map at θ=0: the gap goes 0.157 (L=8) → 0.0692 (L=12), a ratio of 0.44 = (8/12)². Fit: gap = 0.254·(2π/L)².
  **Prediction for L=16: 0.039** (job 848a2647…).
- **Mechanism (to quadratic order, derived here):** 'all' fixes every line holonomy, so each link lies on one
  constrained line. The fibre is then spanned by half-gauge moves: g_x acts on horizontal links only, g_y on
  vertical links only. The curvature change is ∂ₓ∂ᵧ(g_y − g_x), so the Hessian is |∂ₓ∂ᵧψ|², with ψ pinned only at
  block corners. The modes ψ = h(x)s(y) (s vanishing on block rows, h slowly varying) cost ~|h′|² ~ L⁻².
  The combinations f(x)+g(y) are exact gauge, consistent with the measured zero count of 0.
- **Consequence:** a block map that constrains every parallel line (axial-type) has **no volume-uniform fibre
  convexity**, even at zero field. This is the lattice form of axial gauge's known lack of decay. For idea D,
  only averaging maps of corner/[CD] type are viable. The (1−θ/π) law and the θ<π rule are about those maps.
- Odd n: the 'all' gap falls faster (n=1: 0.0105 → 0.0022, about L⁻⁴). It is not needed for the conclusion above.
- *Adversarial note:* "ker J consists of half-gauge moves" is automatic for line constraints: a perturbation with a
  fixed line holonomy is a discrete derivative along the line of a function vanishing at its ends. So an
  eigenvector test of that structure would not discriminate, and it was dropped. The content of the mechanism is
  the L⁻² cost of slowly modulated ψ = h(x)s(y). It is tested by the L=16 prediction (0.039), and analytically by
  the 1D quotient: minimise Σ|∂ₓh|²|∂ᵧs|² with s pinned on block rows, which gives c·(2π/L)².
- *Analytic constant (b=2, back of envelope):* the trial mode ψ = cos(2πx/L)·s(y), with s alternating 0,1,
  has Rayleigh quotient (1/4)(2 sin(π/L))². That gives 0.146 (L=8) and 0.067 (L=12) against the measured 0.157 and 0.069.
  So the scaling and the order of the constant agree; the 4–7% excess is unexplained (a trial-mode
  computation is only a heuristic: ψ is not exactly an eigenmode).
