Directions after Hexagon v7 (upstream main c794c07) and James's message of 2026-11-07. The 2026-11-06 plan (old #1–#8) is exhausted and is summarized at the end of the ranked list.

## Where things stand

Hexagon v7 reduces the mass gap to three inputs (Conclusion; Thm `avgtransfer`, Cor `bpuptuning`). Until all three are proved, the gap is not proved:

- **(H1)** The covariance bound at the coarsest scale, along the tuning. See Prop 130 `prop:bpcoarsezero`, Cor 149 `cor:bpfibH1`, and the sufficient condition in Prop 168 `prop:bppathsuff`.
- **(H2)** The split inequality at weak coupling, uniform in depth and volume. See Thm 114 `thm:bpsplitcond`. Prop 656 `prop:h2volume` shows that the naive version depends on the volume.
- **(E)** A nontrivial continuum limit (`ssec:limitE`).

New in main:
- Commit d90d88b: sharp Gaussian block decay and the Wilson fibre criterion.
- `gaussian-decay-section.tex` and `fibre-heatbath-section.tex`.

Only lean01–lean03 are in main. lean05 is not.

The OpenAI existence claim is still unreleased. If it holds, it bears on (E), not on (H1) or (H2) as Hexagon states them.

## Ranked directions

**A. A 4D SU(2) Gaussian fibre theorem, mapped into Thm 114. Highest priority.**
- Fix the 4D data explicitly: the section, the block map c, the curvature operator and the fibres F_{x'}.
- Prove the geometric hypotheses on the lattice. Get a stable decomposition with constants independent of depth and volume.
- Extend the variance estimate from linear observables to all of Gaussian L².
- Deliverable: explicit constants in the Gaussian model for the three hypotheses of Thm 114: (F) C_F, (X) ρ and (D) R. State the result as "Thm 114 holds with …".
- Build on:
  - d90d88b and `gaussian-decay-section.tex`;
  - `thm:bprowgauss`;
  - `cor:bprowFC` (Gaussian row sums (log b)/b).
- Ask James what lean05 covers (4D aligned covers?) and whether it is coming to main.

**B. Nonlinear (H2), organized around Thm 114 and Thm 181 (`thm:bpdeviation`).**
- **What Thm 114 gives.** It yields S(n',n; C_F/(1−η), (n'/(bn))⁴/(1−ρ²−η)), with η = ρ(RC_F)^{1/2}/n'².
- **Per-scale quantities.** Take n_j = 2n_{j−1}, so n'/(bn) = 1 when b = 2. Set η_{k,j} = ρ√(R C_F)/n_j² and δ_{k,j} = ρ² + η_{k,j}.
- **Target.** C_F ≤ C*, δ_{k,j} ≤ q < 1 and sup_k Σ_j δ_{k,j} ≤ S, all independent of depth and volume.
- **Consequence (checked).** C_A ≤ C*/(1−q), and Π_j C_B ≤ Π_j (1−δ_j)^{−1} ≤ exp(S/(1−q)), using −log(1−x) ≤ x/(1−q) on [0,q].
- **Where the content is.** η is summable automatically because n_j² grows geometrically. So the real condition is Σ_j ρ_j² ≤ S: the exterior constant (X) must decay along the scales.
- **Inputs.** The (X) and (D) bounds have to come from A. Thm 181 handles the deviation term when the chain does not compose exactly.

**C. An (H1) computational campaign on Prop 168, on the cluster with ≤ 64 CPUs. Diagnostics only.**
- **What to measure.** The partial sums Σ_{e≠j} |Cov_{π_V}(a_e^X, a_j^Y)|, against the threshold (1−h_{b−2})/3.
- **Order of runs:**
  1. against coarse volume at fixed depth k;
  2. varying k;
  3. stressing the conditioning V.
- **Where the sum accumulates.** Record the contribution from each of:
  - the 3 links sharing j's base point;
  - nearby shells;
  - distant shells, fitting the decay exponent (Prop 168 needs p > 4, with S_p = 3 + 320(1+1/(p−4)));
  - exceptional configurations.
- **Limits.** Sampled conditioning cannot establish the supremum over V, X and Y.
- **Control.** Prop `prop:bphidden`(ii) is a warning case: a mixture of product laws for which the sum grows with the volume. Run it as a control.

**D. Keep (E) alive with one concrete programme.**
- Fix one gauge-invariant observable, such as a smeared plaquette or a small Wilson loop, and one normalization.
- Targets for its separated correlations:
  - a nonzero limit;
  - uniform moments;
  - rotation restoration;
  - convergence.
- Every use of MIX must cite an independent infrared input, so the gap is never an assumed premise.

**Closed or parked (the 2026-11-06 plan; details in the log):**
- #1 ball certificate: interface drafted; waits on the OpenAI text.
- #2 marginal coupling: parked, gain O(β).
- #3 feasibility: done, L* = 2, heuristic only.
- #4 RP weak gap: no-go.
- #5 Polchinski: dropped.
- #6 Tomboulis: closed.
- #7 3D: folded into the others.

The old #8 deferred (H1) and (H2). That is reversed: they are now the main line. James's dump9 composition results are now the frame for B (Thm 181), not just a no-go.

## Suggested sequencing

| When | What |
|---|---|
| Week 1 | **A:** write out the 4D section, block map, curvature operator and fibres; confirm lean05's scope with James. **C:** first fixed-depth volume scan on the cluster. |
| Weeks 2–3 | **A:** stable decomposition and the Gaussian L² variance bound. **B:** bound ρ_j along the scales using A, and check Σδ ≤ S. |
| Weeks 3–4 | **C:** depth scan and conditioning stress. **D:** pick the observable and list every MIX dependency. |
| On OpenAI release | Read it against D and the #1 interface. |

Sources:
- [KITP working-group statement](https://www.kitp.ucsb.edu/sites/default/files/YM_working_group.pdf)
- [AI Weekly](https://aiweekly.co/alerts/openai-hands-eth-kakeya-yang-mills-proofs-to-reviewers)
- [IBTimes](https://co.ibtimes.com/openai-yang-mills-proof-review-2338)
- [Quantum Research Center](https://quantumresearch.center/news/articles/openai-yang-mills-existence/)
- [Chatterjee, strong mass gap implies area law (arXiv 2006.16229)](https://arxiv.org/pdf/2006.16229v2)
- [BBD Polchinski survey (arXiv 2307.07619)](https://arxiv.org/pdf/2307.07619)
- [sinh-Gordon log-Sobolev (arXiv 2609.28091)](https://arxiv.org/pdf/2609.28091)
- [Langevin dynamics of lattice Yang–Mills–Higgs (arXiv 2401.13299)](https://arxiv.org/pdf/2401.13299)
- [Chevyrev–Shen (arXiv 2503.03060)](https://arxiv.org/abs/2503.03060v2)
- [Dimock, Bałaban RG II (arXiv 1212.5562)](https://arxiv.org/pdf/1212.5562)
- [Wikipedia: Yang–Mills existence and mass gap](https://en.wikipedia.org/wiki/Yang%E2%80%93Mills_existence_and_mass_gap)

## Progress log

- **2026-11-06, #5 dropped.** Falsifier `research/experiments/polchinski2d/` (cluster job 1677160): for the single-plaquette SU(2) measure e^{β cos θ} on S³ (2D axial gauge), the multiscale curvature rate along the heat flow divided by the exact Poincaré gap is 0.59, 0.51, 0.33, 0.11, 0.007 at β = 0.5, 1, 2, 4, 8, which is exponential decay. The flow restores convexity only at t* ≈ 0.4, and the antipodal curvature ≈ −2β is paid for before then. β ≥ 16 rows are numerical underflow. Heuristic proxy for the BBD criterion, not a theorem about it.
- **2026-11-06, benchmark frozen** (see `forest/BALL_CERTIFICATE.md`): SU(2) 0.351, SU(3) 0.77. The 0.348/0.52 figures were a stale README.
- **2026-11-06, #1 interface drafted:** `forest/BALL_CERTIFICATE.md`.
- **#3 running:** `research/experiments/linresp/`.
- 2026-11-06: **#3 done** (job 1676562). L* = 2 for β ∈ [2.0, 2.6] at L=8, T=16. Caveats:
  C(r≥2) is at the noise level, and this is a weak-form periodic-state heuristic, not a
  sup-ω influence. Decision rule (L* ≤ 4) says continue with #1/#2 at small L and favours
  **#4** (weak gap in the RP state) over #6.
- 2026-11-06: **#2 step 5.1 done** (job 1679093). Reduced specification is worse in raw
  Lipschitz form (row 5.40 → 6.54 at L=2, β_W=0.30). See MARGINAL_COUPLING_T1 §6.
  Next is the second-derivative bound for the Z-terms.
- 2026-11-06: **#2 parked.** Variants (a) Hessian bound and (b) one-j marginalization give at most an
  11% drop in total boundary influence (L=2, β_W=0.30), shrinking with L and β (MARGINAL_COUPLING_T1 §6a–b).
  Next: **#4**.
- 2026-11-06: **#4 closed** (`forest/RP_WEAK_GAP_4.md`). The correlation-level RP criterion is a no-go, because
  effective masses are non-increasing (log-convexity), so finite data only bound the gap from above. The chessboard
  version reduces to C2 = (D). This also means #3's L* = 2 has no rigorous force. Next: **#6**.
- 2026-11-06: **#6 audited** (`forest/TOMBOULIS_AUDIT_6.md`). The gap localizes to (★′), monotonicity of the vortex
  ratio Z⁻/Z under c_j ↦ αc_j, or to the global continuation step. (★′) holds exactly in 2D for every group,
  including U(1), so 2D tests are uninformative. The decisive test is 4D U(1) vs SU(2) (small Monte Carlo job).
  This is not a route to the gap.
- 2026-11-06: **#7 started: #2 retested in 3D** (job 1681239, MARGINAL_COUPLING_T1 §6c). The gain is 18% at β_W = 0.30, 11% at 0.45,
  6% at 0.55 (L = 2), shrinking with L. It is O(β) by construction (β²/4 vs β weights). #2 is parked in both 3D and 4D.
  Remaining #7 content: #1 is dimension-free, #5 was dropped. Status: every concrete route is parked or closed except waiting for
  the OpenAI text (#1 interface) and the optional #6 Monte Carlo test.
- 2026-11-06: **#6 closed** (job 1681808, TOMBOULIS_AUDIT_6 §6). The (★′) proxy holds for 4D U(1) in its Coulomb phase (Δ < 0 at 11–16σ),
  so the Tomboulis failure is in the global continuation/eq. 5.21, not in the vortex inequality. Plan exhausted; next is the OpenAI text (#1 interface).
- 2026-11-07: **Plan replaced** after Hexagon v7 (main c794c07) and James's message, with four directions:
  - A: a 4D Gaussian fibre theorem, mapped into Thm 114.
  - B: nonlinear (H2) via Thm 114 and Thm 181.
  - C: diagnostics for the Prop 168 covariance sums.
  - D: a continuum programme.

  Checked that Π C_B ≤ exp(S/(1−q)). The real summability condition is Σ_j ρ_j² ≤ S.
- 2026-11-07: **A week 1 written** (`forest/A_GAUSSIAN_FIBRE_4D.md`). 4D data collected. In the Gaussian model, (F) and (X)
  reduce to the first Wiener chaos with the same constants (Lemmas A1, A2). (D) does not: open item A-D (mixed chaos
  overlaps through non-local coarse factors). Computation plan for λ_F, ρ_1, R_1 set; the lean05 question is drafted, not sent.
- 2026-11-07: **C run 1 done** (job 1695040, `forest/H1_PATHCOV.md`). Negative on method: Σ|Ĉov| equals its time-permuted null
  to < 0.4% at β_W ≤ 1.5, and scales by √2 on half the data, so it is the sampling-noise floor (∝ #pairs/√T), not Prop 168's sum.
  A ~10% excess at β_W = 2.3 is unresolved (possibly autocorrelation). Run 2 needs debiased per-shell estimators before b = 8, random V, or bphidden.
- 2026-11-07: **C run 2 submitted** (ym-pathcov-c2-debiased). Split-block unbiased estimator of per-shell ‖Cov‖_F² with jackknife
  errors, checked first on a synthetic hidden-mixture control with known truth (amp 0 and 0.1, L = 4, 6, 8). Chains: b = 4 at L = 4, 6, 8
  and b = 6 at L = 4, β_W ∈ {1.5, 2.3}, 3000 measurements each.
- 2026-11-07: **A/B: Gaussian ρ_1 job submitted** (ym-gaussrho-a1, job 1696054). Abelian Gaussian model: computes the first-chaos
  (X) constant (exact on all of L² by Lemma A2) via div-free path currents, the Hodge Green's function and the cube Hessian
  pseudo-inverse. Covers b = 3, 4, n'/b = 1, 2, 3, L = 6 vs 8 and cube offset. Selftests: Hodge identity, cube block, 0 ≤ ratio ≤ 1.
  The decision is ρ_1 vs n'/b: decreasing feeds B (Σ ρ_j² < ∞); flat near 1 blocks B as framed.
- 2026-11-07: **A week 2, (D):** Lemma A3 gives the exact chaos expansion w = Σ_{S≠∅} B^{⊗S}⊗A^{⊗rest}⊗I, which is orthogonal in |S|. A-D reduces to one combinatorial bound c_k(s) on symmetrized contractions. It is still open; the naive bound is C(k,s)². See forest/A_GAUSSIAN_FIBRE_4D.md §6. The C run 2 and ρ_1 jobs are still running.
- 2026-11-07: **D part 1:** the observable is fixed (24-plaquette action density at 0, Z_k = 𝒩_k(ε)^{−1/2}) and every MIX use in §limitE is listed (forest/D_E_OBSERVABLE.md). Finding: the *nontriviality* of the limit currently comes from the gap (Prop limitEnorm(ii)), which violates the D rule. Next: D1 (γ_k audit of limitEtwo), D2 (MIX-free nontriviality), D3 (audit limitEexact).
