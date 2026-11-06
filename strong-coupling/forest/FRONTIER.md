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

## Lead: relativized Dobrushin criterion (Nov 2026, not yet read in full)

**Clarification (adversarial re-read of CALM_LINK §3–4).** Averaged one-point (R1) is *already* available on the torus:
chessboard rarity gives λ|L₁|p_nc ≤ 10⁻³ at β_W ≳ 2·10⁹⁵ (`d5/bootstrap.out`). So the wall is not (R1) itself. It is
that the Dobrushin–Shlosman cube coefficient takes a supremum over boundary conditions, and uniform-field boundaries
break that supremum.

**Candidate replacement.** Baffioni–Kuna–Merola–Presutti, "a relativized Dobrushin uniqueness condition"
(SFB701 preprint 07-042, https://sfb701.math.uni-bielefeld.de/preprints/sfb07042.pdf). It asks for the Dobrushin
contraction only on "good" boundary conditions. Failures are compensated by a priori bounds on the probability of bad
conditioning, taken under the Gibbs measures generated by the admissible class of boundaries.

**Next step.** Check whether its a priori estimate can be supplied by (i) the calm contraction C_* (`D5.md`) and (ii) a
rarity bound under the *admissible-class* measures rather than the torus measure. If (ii) needs a bound uniform over
the class, the uniform-field obstruction reappears unless the class excludes it. Watch for the hidden-bit trap
(`ENVIRONMENT_OBSTRUCTION.md`). Status: **open lead**.

**Assessment after reading §§2–3 of BKMP (same day): reduces to the same wall, so the lead is dropped.**
- *Assumption 1* is the Dobrushin bound for boundaries in a good set G_x.
- *Assumption 2* asks for ∫_{G_x^c} p_{C_i}(dt|s) ≤ ζ ⟨|t_x|²⟩², **for every s ∈ X**, where X is the restricted ensemble
  (a hard-core constraint that allows only small deviations from one stable phase).
- The "relativization" therefore controls bad boundaries near the edge of the constraint. They are rare *uniformly
  over the restricted ensemble*. Configurations outside X are excluded by construction, not averaged away.
- Part 2 of the paper puts everything outside X into a Pirogov–Sinai contour/Peierls analysis.

Translated to our setting:
- X is the all-calm ensemble. Uniqueness there is essentially the calm contraction C_* we already have.
- A uniform-field boundary with B ≫ K lies outside X.
- The large-field contours then need a Peierls bound under arbitrary contour exteriors. That is Bałaban's large-field
  expansion, the same wall as (R1)/D7.

Useful as a clean decomposition ("restricted-ensemble Dobrushin + contour Peierls"), but it does not avoid the
non-abelian large-field Peierls estimate. *(Reading of the extracted text, not a full proof check.)*

## 2026-11-05 — Announced OpenAI existence proof for 4D quantum YM (secondhand, UNVERIFIED)

Source: KITP working-group statement, "Proof by OpenAI of the Existence of Quantum Yang-Mills theory",
shepherding committee Douglas, Gross (Chair), Jaffe, Klebanov, Nekrasov, Oct 2026
(https://www.kitp.ucsb.edu/sites/default/files/YM_working_group.pdf). Commentary: Woit
(https://www.math.columbia.edu/~woit/wordpress/?p=15948). The text I extracted from the PDF says:
- The claim is existence of the continuum limit, "including coupling to fermions". Bałaban's RG was among
  the key ingredients, used "to follow the theory all the way to infinite distances".
- "confinement and the mass gap remain very important open challenges."
- The manuscripts total 1147 pages and come with a Lean formalization. Neither is public yet; the committee
  will release them after it has checked them.

Implication for this programme (my reading only): if the claim holds, it settles our "scale composition" and
"continuum limit" walls (§10.3, §10.6), very likely by full Bałaban multiscale control. Our one-scale results
would then be superseded on existence. The Clay-relevant open piece is the **mass gap**, i.e. Conjecture L
in the continuum. Proof text not available, so this changes nothing proved here.
- Follow-up (2026-11-05): OpenAI's own post (https://openai.com/index/building-scientific-understanding-with-academic-communities/)
  gives these dates: natural-language YM existence proof on 2026-09-09, a QCD proof on 2026-09-18, and the Lean formalization
  finished on 2026-10-16. All YM materials are withheld pending the Gross panel. No arXiv preprint found. Mass gap not claimed.

## 2026-11-05 — Literature bearing on §11 idea 3′ (stability of mixing under perturbation)
- Crawford–De Roeck (arXiv 1612.02620). Uniqueness for ferromagnetic Glauber dynamics is stable under small
  perturbations throughout the *whole* uniqueness phase. The input is weak spatial mixing on finite boxes.
  This is the closest analogue of 3′ gap (a): stability needs mixing of the unperturbed system, not a
  perturbative regime. Caveat: ferromagnetic, discrete spins, monotonicity is used. None of that applies to SU(N).
- Martinelli–Olivieri–Schonmann (mp_arc 93-154). In d=2, weak mixing implies strong mixing for squares, with
  counterexamples in higher d. Relevant warning: in 4D, mixing on cubes does not imply strong mixing, so the
  finite-box certificate must be the strong (all-boundary-condition) kind. Our GDS is sup over ω, so it is already
  of that kind.
- Martinelli–Olivieri SSM definition. These are the standard references for idea 3; I found no paper with a
  volume→boundary improvement of perturbative stability for finite-size conditions.

## 2026-11-05 — Chatterjee, "strong mass gap ⇒ confinement" (arXiv 2006.16229, CMP 385 (2021))
- Exponential decay of correlations *uniformly over arbitrary boundary conditions* ("strong mass gap"), with a
  compact connected G of nontrivial center, implies Wilson's area law.
- Our Conjecture L, GDS, and the Doeblin/shell route are all sup-over-ω statements, so they target the *strong*
  mass gap. That is stronger than what Clay needs (decay in the transfer-matrix / periodic state), and it would
  also imply confinement.
- Adversarial implication: we may be aiming at a harder target than necessary. The weak form only needs decay
  for the infinite-volume periodic/free state. Its natural tools are reflection positivity and transfer-matrix
  spectral methods, which our finite-box criterion does not use.
- No 2025–26 converse or all-β result was found in a standard search.

## 2026-11-05 — Weak-form route via finite-size spectral-gap criteria: dead end as is
- Knabe (1988); Lemm (threshold 3/n); Anshu (~1/ℓ²); Lemm–Xiang (arXiv 2112.07756); Lemm–Mozgunov (arXiv
  1801.08915, open boundaries). All of these require **frustration-freeness**, which is called essential.
- The Kogut–Susskind Hamiltonian, i.e. the transfer matrix of lattice YM, is not frustration-free: the electric
  and magnetic terms do not share ground states. No non-FF finite-size gap criterion was found.
- Verdict: the transfer-matrix finite-size route is a dead end without a new non-FF criterion. Proving such a
  criterion would be a major result in quantum many-body theory by itself. Dropped.

## 2026-11-05 — Recent rigorous lattice gap results (none reaches the window)
- Nissim, "U(N) lattice Yang–Mills in the 't Hooft regime" (arXiv 2510.22788).
  - Claims: a mass gap, a unique infinite-volume limit, and a new proof of the large-N limit, via cluster
    expansion plus Langevin dynamics on a U(1)×SU(N) random-environment model.
  - It notes that the earlier SU(N)/SO(N) gaps relied on Bakry–Émery, which fails for U(N).
  - From the abstract, the regime is presumably small 't Hooft coupling, i.e. strong coupling. Not checked in the full text.
  - Not a window result.
- Chatterjee 2026 (Prob. Math. Phys.): Gaussian scaling limit of SU(2) Yang–Mills–Higgs. Higgs matter,
  not the pure-YM gap.
- A reported Oct-28-2026 preprint: a volume-uniform gap in a strong-coupling window only. Primary source not found.
- Many Zenodo/SSRN "solutions": not credible, ignored.

### 2026-11-05 sweep: intermediate-coupling mass gap literature
- Still no accepted rigorous gap for 4D lattice YM at intermediate β. Chatterjee 2401.10507 (an SU(2) YMH scaling limit, Higgs-mechanism mass, Gaussian limit) is not pure YM.
- Zenodo 18923847-family (self-published, Mar 2026): claims an SU(2) Wilson gap for β ≥ ~9587 via block RG.
  - A fixed-large-β gap would be *exactly* our blocked scale-composition step: it would have to beat the asymptotic-freedom flow to the window.
  - Red flag: for fixed β the RG flow reaches strong coupling after ~β/b₀ log steps, so any proof must handle the window. Check before trusting. Status: unverified, likely flawed. Checked record 18923847: it is the CFUT item, a 12.9 kB markdown "proof" with no RG treatment. Dropped.
- Others (CFUT, "Universal Cascade", PhilArchive): non-rigorous, ignore.
