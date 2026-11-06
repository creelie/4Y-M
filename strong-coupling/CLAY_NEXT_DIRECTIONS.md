I read the git log (the §9–§12 run), dump1, dump9, dump10, `strong-coupling/OPEN.md`, `forest/FRONTIER.md`, `POLYMER_ARCHITECTURE.md` §10–12 and `how-to-use-cluster.md`, then searched the web. Here is a ranked list of what to explore next.

## Where things stand

Our own §10.6 and §11 notes have already reached the main point. **Every route we have tried ends up at the same intermediate-coupling window.** Strong-coupling, curvature, Dobrushin, Higgs-continuation and polymer arguments all stop where the running coupling is O(1). Two outside facts change how we should spend effort:

- **OpenAI's existence claim.** It covers existence only, not the mass gap, and the proof is still withheld by the Gross panel with no release date. If it holds, it covers roughly Target 3 and OPEN.md problems 2–3 (except the window). What remains is a gap at the terminal scale.
- **We are aiming at a harder target than Clay needs.** GDS, D(n,α) and Conjecture L all take a supremum over boundary conditions. That is Chatterjee's "strong mass gap", which would also imply confinement. Clay only needs the gap in the periodic or reflection-positive state.

## Ranked proposals

**1. A box certificate that covers a whole ball of actions (our idea 3″, made concrete). Highest Clay value of anything we can do.**
- The finite-box criterion is the only route we know that avoids the window (§11, idea 3). Lemma S′ failed because the certificate was for one action Φ, so stability brought the |V| factor back.
- The fix is to certify GDS(V,ε) uniformly over a parametrised family: Wilson plus a few loop and adjacent-plaquette terms plus a remainder of norm ≤ δ. This is exactly the object a Bałaban-type construction outputs at the terminal scale.
- First deliverable: a precise statement of the form "effective action in ball B(Φ₀,δ) in norm N, plus a robust certificate on that ball, implies a gap". Write it so that we can plug it in the day OpenAI's proof is released.

**2. A certified finite-box marginal coupling (dump10 T-1 to T-3). Concrete, and it feeds #1.**
- Integrate out the links next to the boundary exactly using SU(2) characters, then maximally couple the far-field marginal. Certify it over all boundary holonomies with interval arithmetic.
- Run it first as a benchmark at d=4, L=2, β_W=0.30, then try to beat 0.351 (the Hexagon v6 SU(2) number).
- This is also the engine #1 needs. Freeze the benchmark first: the README says 0.348/0.52 and v6 says 0.351/0.77.

**3. A feasibility study on the cluster before committing to #1 and #2.**
- With the exact-marginal coupling, not the naive linkwise one, measure the GDS/influence coefficient against box size at β_W ≈ 2.0–2.6, on both Wilson and MCRG-blocked actions.
- Output: the smallest box L* at which contraction appears numerically, and the size of the action family that needs certifying.
- If L* ≥ 6 in 4D, a direct certificate is out of reach and we need a dimension reduction (class functions, symmetry). Better to learn that in a week than in two months. This is cheap CPU work.

**4. Go after the weak mass gap with reflection-positive and transfer-matrix tools.**
- Our infrastructure only partly uses these yet: the slice ratio, Theorem 46, and the tree bounds (T1)–(T3).
- Candidate: a chessboard/reflection-positivity argument that turns a gap at one block scale in the periodic state into decay, avoiding the supremum over boundaries entirely.
- The finite-size criterion for non-frustration-free Hamiltonians is dead, but a correlation-level (not Hamiltonian-level) criterion in the reflection-positive state has not been ruled out.
- Worth one focused session to find either a no-go result or a lemma.

**5. Polchinski-flow log-Sobolev (Bauerschmidt–Bodineau–Dagallier) for lattice gauge fields.**
- Idea 4 killed curvature of a *single* effective action through the antipodal configuration.
- The Polchinski method needs something different: a bound on the Hessian of the renormalised potential *relative to the covariance increment at each scale*. Non-convex potentials are allowed if the flow smooths them, as in their φ⁴ and sinh-Gordon work.
- Cheap falsifier: 2D U(1) and SU(2), which are exactly computable. Check whether the multiscale condition holds beyond the Shen–Zhu–Zhu threshold. If it does, this is a genuinely new non-perturbative tool. If the antipodal obstruction survives, drop it quickly.

**6. Transferring the Migdal–Kadanoff result: audit Tomboulis's 2007 confinement claim.**
- OPEN.md notes that confinement is known at every coupling in hierarchical models, with no transfer to Z⁴.
- Tomboulis's interpolation argument tried exactly this transfer, and the community considers it to have a gap (Ito–Tomboulis follow-ups).
- Pin down precisely which inequality fails, and test it numerically in 2D/3D. If the failure can be localised, that is either a repair route or a sharp obstruction worth publishing. This is a literature task plus small numerics.

**7. 3D as a proving ground.**
- The 3D mass gap is also open, and 3D is super-renormalisable, with a shorter and easier RG.
- We already have the 3D forest (β_W < 4/7). Running #1, #2 and #5 in 3D first gives credible milestones and catches failures at lower cost.

**8. Lower priority or defer:**
- **Target 1, (H1) on every torus, and Target 2, (H2):** both are Bałaban-control statements. If OpenAI's construction is correct, it supersedes our own attempts at them. Wait for their text.
- **(N_u) nontriviality:** existence claims normally include it. Defer, but keep the slice-ratio criterion ready as a checkable interface.
- **More h_F tightening or Lean arithmetic:** no Clay value. James's composition results in dump9 are useful as a theorem that the site-additive approach cannot work, nothing more.
- **Already dropped and should stay dropped:** BKMP relativised Dobrushin, frustration-free finite-size criteria, Higgs continuation, curvature of the blocked action, and the shell-transport chain (§10–12, which is one-scale and already within Bałaban's framework).

## Suggested sequencing

| When | What |
|---|---|
| Week 1 | #3 feasibility runs on the cluster, the #1 interface theorem on paper, and the #5 2D falsifier, in parallel |
| Weeks 2–4 | #2 certification at the benchmark, then #4 or #6, depending on what #3 shows |
| On OpenAI release | Read it against the #1 interface, specifically: what norm and what δ their terminal-scale effective action comes with |

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
