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

## Literature: Bałaban–Dimock method (searched 2026-11-05)
Dimock, "The Renormalization Group According to Balaban", I–III: arXiv 1108.1335 (small fields), 1212.5562 (large fields),
1304.0705 (convergence). These treat φ⁴ in d = 3. Large fields are handled "crudely" there, using the stability of the interaction.
A 2025 Roma Tre talk extends this to U(1)-Higgs in d = 2, 3 and the O(4) σ-model. I found no published 4D non-abelian large-field
treatment beyond Bałaban's original series. That is exactly the wall for (H1) and (H2).

## Queued test (not run: CPU budget full)
`uniform-field/gapscan.py` + `plan-gap.json` (reviewed). It measures the gauge-reduced fibre Hessian gap as a function of coarse angle θ ∈ [0, π]
and block size b. The question: does small-field convexity need an angle margin θ ≤ π − δ?

## External review (James's agent, 2026-11-04): what I accept
- **Floor (item 7):** agreed that it is not a demonstrated defect. `FLOOR_1e-4.md` already says "likely ridge, not proved", and nothing
  shows that it affects Romanov's exact rational/SOS certificate (h_F ≤ 14.420001). Any objection would need to name a specific matrix and
  the theorem that depends on it.
- **"No progress on the core difficulty" (items 10, 16):** too categorical. The thresholds belong to different methods: 2/9 single-link,
  1/3 forest, 0.351 plaquette-block (v6). The compact-coupling G1 theorem and the forest combinatorics are real content. The wall
  (large fields, scale composition, continuum limit) is unchanged.
- **H1 sup (item 12):** never claimed disproved. The uniform-field work answers the reviewer's priority 2. The exterior field
  has a Nielsen–Olesen instability, and the block constraint removes it exactly for coarse angle θ < π, with gap ≈ (2–2.4)/b²·(1−θ/π)
  (job 1663501). So uniform fields are not a counterexample on calm data. Whether that holds near and above θ = π is open.
- **Not yet verified by me (do not cite):** Rajabu, 0.419 (Oct 2026); Douglas et al., Lean lattice gauge theory; Douglas, Nature Rev. Phys. 2026.
- **Priority 3 (scale-composition theorem)** is the same wall as Conjecture L / averaged (R1) here. Agreed that it is where the effort should go.
- **Checked 2026-11-05:** Douglas `lgt` (Lean reservoir @mrdouglasny/LGT) exists. It covers U(n) Wilson theory at β < 1/(32n(d−1)), proved by Dobrushin
  (after Chatterjee, Ch. 16), with 0 sorry and only foundational axioms reported. For n = 2, d = 4 that is β < 1/192. The normalisation still
  has to be matched against β_W = 2/9 before comparing. Rajabu 0.419: not found by search, remains unverified.

## Literature check 2026-11-05 (second pass)
- No 2024–26 paper found that extends Bałaban–Dimock large-field control to 4D non-abelian YM. The Roma Tre reading
  list (Tanimoto) is rereading Bałaban's gauge papers and Dimock's scalar QED d=3 in 2025, so it is review work, not new results.
- Adhikari–Cao (Ann. Probab. 2025): weak-coupling correlation decay, **finite** gauge groups only. Its method
  (rarity of vortices) has no direct SU(2) analogue, because continuous groups have spin waves. Logged, not pursued.
- An unreviewed listing (SynapseSocial) claims a "Bałaban on S³, large-field region empty inside the Gribov horizon" result.
  Unverified. Treat it as noise unless an arXiv record turns up.
- Internal: 'all'-type (every-line) block maps have an L⁻² half-gauge soft mode at zero field (H1_AVERAGED.md).
  Restrict idea D to single-average-per-block-link maps ([CD] V♯, corner).
