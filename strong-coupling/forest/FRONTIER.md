# Frontier log (running)

This is a working log for the three dump1 targets: (H1) on every torus, (H2) at weak coupling, and (N_u) / a nontrivial limit.
Each is at the level of the Clay problem. Nothing in this file is a proof unless it is labelled **proved**.

## Literature status (searched 2026-11-05)
- **Verified results.**
  - No peer-reviewed proof of the 4D Yang–Mills mass gap exists.
  - The standard UV-side references are Bałaban's 1980s series and Dimock's exposition, "The Renormalization Group According to Balaban" I–III (small fields, large fields, convergence).
  - Area law at strong coupling: Cao–Nissim–Sheffield, arXiv 2505.16585 (2025), extend it from β < c_d/N to β < c_d for U(N).
  - Large-N surface sums: Cao–Park–Sheffield; Borga–Cao–Shogren-Knaak.
  - None of these touches weak-coupling (H1) or (H2).
- **Unverified claims, do not cite:**
  - Zenodo, philarchive and vixra preprints claiming a full proof.
  - A news-style report (vocal.media) of an "OpenAI existence proof, Oct 2026". I found no primary source for it.

## Idea triage for Target 1, (H1)
| # | Idea | Verdict |
|---|---|---|
| A | Sup over all V | **dropped as target.** Large-angle uniform fields have unstable directions in the fibre (job 1662049). |
| B | Derive averaged (H1) from total-measure decay | **dead.** The absolute value sits inside the V-integral, and a conditional covariance is not the total covariance. |
| C | Bakry–Émery on SU(2)ⁿ, with Ric = 1/2 > 0 | **only strong coupling.** The negative part of the Hessian of βS is O(β). This is already known, and is the regime of `hW`. |
| D | Calm-V fibre convexity + Helffer–Sjöstrand: H\|ker J ≥ c > 0, uniformly over calm V, gives decay at the Green's-function rate | **alive, but it is Bałaban's small-field programme.** Needs gauge fixing and the large-field polymer expansion. |
| E | Uniform-field obstruction restricted to calm V | **not an obstruction at linear order** (L=16). At L=32 (job 1662119) the onset is exactly coarse angle π, the aliasing threshold. Calm uniform data are linearly stable in their fibre. |

## Adversarial notes on our own outputs
- `FLOOR_1e-4.md`: downgraded to "likely explanation" after the exact rational check (`c7febc5`).
- `HW_SU2.md`: the arithmetic checks out:
  - 6 plaquettes × 3 links gives 18 neighbours;
  - a staple Lipschitz constant of 4|β| with L ≤ 1/4 gives C_ij = |β|;
  - so the row sum is 18|β| < 1 ⇔ β_W = 4β < 2/9.

  Metric caveat resolved: `finite_dobrushin_comparison` (dump2 l. 3287) takes a generic metric `dX`, and the
  1/4 bound is chordal (`zonal_final`, dump2 l. 9044), so choosing chordal `dX` is consistent. Remaining
  caveat: the Lean condition is a weighted column sum, Σ_i w_i C_ij ≤ θ w_j. With uniform w and the symmetric C here,
  that equals the row sum.
