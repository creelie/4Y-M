# One-scale large-field polymer expansion: what it must prove (Nov 2026)

Labels: **proved**, **sketch**, **heuristic**, **open**. Nothing here closes D.

## 0. Why the architecture changes
Every attempt to get (R1) under a fixed exterior failed: localcalm, chessboard, Peierls (`CALM_LINK.md` §4),
Dobrushin–Pechersky, and BKMP (`FRONTIER.md`). The calm-collar variant fails as well, because rarity has to hold given the
coupling's past, and that is the hidden-bit trap. In Bałaban/Dimock, large-field regions are **summed** on the torus, never conditioned
on. So the target is a polymer representation, not a cube condition.
Reference: Dimock arXiv:1212.5562 (φ⁴₃ UV). It covers the scalar case only, where the large-field cost comes from λφ⁴. In
the gauge case the cost is the action itself, e^{−βK²/C} per large plaquette.

## 1. Statement needed (open)
On the torus Λ at one block scale, write Z = Σ_Ω Z(Ω). Here Ω is the closure of the set of non-calm blocks plus a collar of width w,
and Z(Ω) integrates over configurations that are non-calm exactly on Ω and calm on Λ∖Ω.
(P1) **Peierls by summation:** Z(Ω)/Z_calm(Ω) ≤ e^{−c β K² |Ω|/w^d}. Here Z_calm(Ω) is the same integral with Ω also calm,
     and the bound is uniform in the configuration on Λ∖Ω **only through its calm collar**.
(P2) **Analytic calm region:** log Z(calm on X | collar data) has a convergent polymer expansion with
     exponentially decaying activities (the chart / `chartdecay` machinery, already in the paper at sketch level).
(P1)+(P2) together give a standard polymer gas, and a cluster expansion then gives decay of correlations on the torus. That
replaces the cube coefficient and D4/D7 at one scale. D8 (composition across scales) stays open.

## 2. The calm-completion step inside (P1) (heuristic)
(P1) needs a calm configuration on Ω compatible with the calm collar, with action ≤ C|Ω|·(small).
- Holonomy: the collar is calm, but loops around Ω can have O(1) holonomy. Spread over area R², that is plaquettes of size ~θ/R²,
  so it stays calm for large Ω. For Ω of a few blocks, θ ≤ π gives plaquettes of size ~π/b², which is calm only if K ≳ π/b².
  **Potential obstruction for small polymers**, where the threshold K ≈ 4.8·10⁻¹⁸ is far below π/b².
- Topology: collar boundary data can carry instanton winding. A fill of size R has plaquettes ~1/R², so the fill is again
  impossible to make calm for small Ω.
- In Bałaban, small polymers with non-trivial collar holonomy simply *are* large-field, because Ω is enlarged until the
  collar is flat-ish. So the polymer definition has to be "closure under fillability", and its entropy has to be counted.
  **This counting is the concrete open lemma.** It has to beat e^{−cβK²} per block. Fillability-enlargement can multiply
  |Ω| only by a constant if holonomy ≤ C·(collar plaquettes)·(perimeter), which is a non-abelian Stokes bound that holds on
  calm collars (sketch).

## 3. Next step
Prove the non-abelian Stokes/fillability lemma at one scale. Statement: a w-collar with all plaquettes ≤ κ around a region
Ω of diameter R has a calm fill with plaquettes ≤ C κ (R/w)^{?}, in the absence of instanton winding. Then check whether
the exponent in (R/w) is compatible with the threshold K. This is testable numerically on small SU(2) lattices
(minimise the action of the fill given random calm collars).

## 4. Adversarial correction to §2 (heuristic)
- **Holonomy obstruction: spurious for ball-like Ω in d=4.** The collar ≅ S³×[0,w] is simply connected. Every loop in
  it bounds a disk inside the collar, so non-abelian Stokes gives holonomy ≤ C·κ·(area in the collar). A loop does not need
  to pass through Ω to be bounded. The calm data on the collar is therefore close to a pure gauge g: S³→SU(2).
- **Winding is the real obstruction.** g has a degree k (π₃(SU(2)) = ℤ). By Lüscher admissibility, if all plaquettes are
  small, the lattice charge is an integer, so a fill must carry charge k. Its action is ≥ 8π²|k|/g², i.e. ≳ β_W·const per
  unit of winding, and the fill can be calm only if its size is ≳ K^{−1/2}.
- So non-trivial winding is suppressed by the action itself, e^{−cβ_W}: it is just one more large-field cost. Non-ball Ω
  (handles, i.e. non-simply-connected collars) bring the holonomy problem back on non-contractible loops. That needs the
  enlargement rule "fill the holes of Ω", and its entropy is a constant factor per block.
- Revised open lemma: **calm fill for k = 0 with plaquettes ≤ Cκ, C independent of R.** The danger is C growing with R
  (elliptic estimates for the fill are fine for Laplace-type problems, but the fill is non-linear). This can be tested
  numerically.

## 5. Literature lead for the fill lemma: quantitative null-homotopy (sketch, not checked in detail)

Chambers–Manin–Weinberger, *Quantitative nullhomotopy and rational homotopy type*, GAFA 2018 (arXiv:1611.03513):
a null-homotopic L-Lipschitz map S^m → S^n has a null-homotopy of thickness C(m,n)(L+1) and width C(m,n)(L+1)².
Companion paper: Chambers–Dotterrer–Manin–Weinberger, *Quantitative null-cobordism*, JAMS 2018.

How this maps onto the fill. The collar's gauge function is g: S³_R → SU(2) ≅ S³ with degree k = 0. Rescaled to the unit
sphere it has Lipschitz constant L ≈ κR, where κ is the per-link deviation.
- **L ≲ 1 (κR ≲ 1).** g stays near a constant up to O(L). A geodesic radial contraction then fills with link
  deviation O(κ), so C is independent of R. This is elementary.
- **L ≫ 1 (κR ≫ 1).** Put the CMW homotopy in the radial direction. Spatial gradient ~ L/R ~ κ, which is fine. But the
  width is ~ L², which forces a radial Lipschitz of order L²/R ~ κ²R. That is **not** uniform in R.
  Either a better (linear-width) bound is needed for maps S³ → S³, or the polymer counting has to pay for
  κR ≫ 1 collars separately. Gromov's linear conjecture is open in general. For S³ → S³, which is degree theory only,
  linear null-homotopies may be easier: **open, to check**.
- The numerical test (job 2023a16d…, 3D cube) cannot see the L² issue, because π₂(SU(2)) = 0 and its collars have L = O(1).

## 6. Correction to §5: the parameter is κR² (curvature), not the gauge Lipschitz constant (sketch / open)

**Scalability lead.** Berdnikov–Manin, *Scalable spaces* (Invent. Math. 2022, arXiv:1912.00590), and Manin, *Plato's cave*
(arXiv:1801.00335), give O(L)-Lipschitz null-homotopies X → Y for scalable Y. Confirmed from the abstract: the scalable class includes spheres and Lie groups. Precise theorem form (linear null-homotopy
constant, domain hypotheses) still to be read from the paper. If S³ = SU(2) is scalable, the L² width in §5 disappears for the *flat* part of a collar.

**Adversarial point: the gauge function is the wrong object.**
- **Exactly pure-gauge collar.** If U_xy = g_x g_y⁻¹ on the collar, then any extension of g to interior sites gives
  plaquettes exactly 1. Lipschitz bounds are irrelevant, and no winding question arises, because k counts curvature
  and the collar's curvature is zero.
- **Calm collar.** Plaquettes ≤ κ only. In a radial gauge from a pole of S³_R, the collar connection has |A| ≲ κ·dist ≤ κR.
  The cone fill A_t = tA has F = tF_A + (t² − t)A∧A, which gives **|F| ≲ κ + (κR)²**.
  - So a naive fill is calm iff κR ≲ 1, i.e. κR² ≲ R. This regime is not uniform.
  - Uhlenbeck-type gauge needs ∫_{S³_R}|F|^{3/2} ~ κ^{3/2}R³ small, which also fails for large R.
- **Consistency check (no contradiction).** For any fill, ∫_B tr F∧F ≡ CS(collar) mod ℤ. Both sides are O(κ²R⁴), and a
  calm fill has ∫|F|² ≤ κ²R⁴. So there is no dimensional obstruction to a calm k = 0 fill.
- **Revised open lemma (multiscale).** For a collar on S³_R with |F| ≤ κ and k = 0, there is a fill with |F| ≤ Cκ and C
  independent of R. A plausible proof is a dyadic decomposition: on each scale-ℓ sub-block, κℓ² ≪ 1 holds once ℓ ≪ κ^{-1/2},
  which allows local Uhlenbeck gauges. The transition functions then form a Lipschitz map into SU(2), which needs scalability to
  be null-homotopic with linear constants. **Open.**
- **Relevance check for (P1).** Calm regions in the polymer expansion have κ ~ g·polylog, so κR² ≫ 1 holds for large
  regions. The hard regime is generic, not an edge case.

**Numerics.** Job 2023a16d (κ ≤ 0.01, R ≤ 10, so κR² ≤ 1) only tests the easy regime. The follow-up needs κ = 0.1–0.3 at R = 10
(κR² = 10–30).

## 7. Calm-fill numerics, 3D SU(2) cooling (heuristic)

Jobs 2023a16d (κ = 0.001, 0.01) and d2257246 (κ = 0.1, 0.3), with R = 3…10. Collar links are a smooth gauge transform composed
with iid noise of size κ. The table gives interior max plaquette / boundary max plaquette:

| κ | R=3 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|
| 0.001 | 0.57 | 0.53 | 0.43 | 0.44 | 0.36 |
| 0.01 | 0.57 | 0.53 | 0.43 | 0.44 | 0.36 |
| 0.1 | 0.53 | 0.38 | 0.40 | 0.47 | 0.36 |
| 0.3 | 0.57 | 0.42 | 0.44 | 0.53 | 0.39 |

No growth with R: C ≤ 1 in every run.

**Adversarial caveat.** With iid noise, the curvature integrated over the collar is a random walk of size ~κR, not
κR², so these collars are nearly the easy case. Coherent abelian fields are trivially fillable by the same field.
The real test is coherent *noncommuting* curvature. That is the hedgehog collar: each face carries field κ with
generator n̂·σ. It runs as job 0e0d65f7 (κ = 0.05 and 0.2, so κR² reaches 20).

**Hedgehog v1 (job 0e0d65f7): invalid as a test.** The generator jumps where faces meet, so the edge plaquettes are of order
κR. bd_max therefore grew linearly in R, from 0.08 to 0.36 at κ = 0.05, and the collar was not calm. Measured against the
intended κ, the interior grows like ≈ 0.3κR. That is ambiguous: it may only reflect the non-calm edges.

**The real threat, identified.** The test case is a reducible connection with direction field n̂ of degree N, covariantly constant, and
uniform curvature F = κ n̂·σ. Total flux is 4πN, with N ≈ κR²/4π. This collar is calm. Filling it in the ball means unwinding
charge N through non-abelian configurations, and **I know no bound on that cost**: the "monopole charge" is not
topological in SU(2) because π₁ = 0, but non-topological does not mean cheap.
- **3D check: job af45d8e1** (`fill_mono.py`). Valid only if flux_err ≈ 0. Its output includes the depth profile.
- **4D analogue.** On S³ the corresponding object is a calm collar carrying a large, non-quantised Chern–Simons
  density. Same open question.

## 8. Calm fill by distributed lumps (sketch; prediction for job af45d8e1)

**3D, degree-N reducible collar.** Place N unit-charge cores, each of radius ℓ ≈ κ^{-1/2}, uniformly in the ball B_R. Away from the cores,
take a reducible connection with F = f n̂, where n̂ has a degree-1 defect at each core. The core count adds up to degree N on the
boundary, as it must. Gauss's law gives |f(r)| ≈ N r/R³ ≈ κ r/R ≤ κ at radius r. The cores are smoothed abelian monopoles,
capped in SU(2), with |F| ≈ 1/ℓ² ≈ κ. ∇n̂ does not enter F. **This gives a calm fill with C independent of R and N.**
So the monopole collar is *not* a counterexample. Prediction for job af45d8e1: the ratio stays O(1).

**4D.** Topological charge can be stored in lumps of scale ℓ ≈ κ^{-1/2} at density ≤ ℓ^{-4}, so the ball holds up to about κ²R⁴ units.
The collar's Chern–Simons term is |CS| ≲ κ²R⁴ (§6). Dimension counting therefore always lets the charge be
placed calmly, and this holds even for k ≠ 0 with |k| ≲ κ²R⁴. **Correction to §4:** winding is an obstruction only when
|k| ≫ κ²R⁴. The calm-fill cost is action, not topology.

**Remaining gap (open).** This handles reducible or "abelian-like" collars. A general calm collar needs a multiscale
construction: local Uhlenbeck gauges on blocks of size ℓ, with transition data glued by the linear null-homotopies of scalable
targets (§6), and charge deposited in lumps. Next steps: look for an existing "small-energy extension / Chern–Simons filling"
lemma, for example in the Uhlenbeck–Wehrheim boundary-value literature, or in Bałaban's own large-field
fill (*Comm. Math. Phys.* 1985–89), which may already do exactly this on the lattice.

**Literature pointer, to check next.** In Bałaban's scheme (*Large field renormalization I/II*, CMP 122, 1989; exposition: Dimock,
*The renormalization group according to Bałaban*, arXiv:1108.1335 and 1212.5562), large-field regions are treated "crudely".
Their weight is bounded by stability plus the e^{-p(g)/g²} suppression per large cube, and the matching lower bound on Z comes
from the global small-field expansion. **If that is right, no regionwise fill lemma is needed. The fill serves only as one
way to get the denominator.** Read Dimock §§ on large fields to decide whether (P1) can follow that route, which would
make §§4–8 optional.

**Checked, from Dimock's exposition (arXiv:1212.5562, scalar φ⁴₃, §1–2):**
- **Suppression is local.** Every large-field cube carries a factor e^{-O(1)p₀²} from the integrand itself, and regions Ω₁ᶜ
  are summed over the whole lattice. This is (P1) by summation, with no conditioning.
- **No fill appears.** The normalisations are explicit Gaussian constants (eqs. 8, 14–16). In the readable part there is no
  lower bound on Z obtained by extending small-field configurations.
- **The boundary dependence is real, but handled differently.** The small-field part sits on Ω₁ with Dirichlet data on ∂Ω₁,
  and the bonds across ∂Ω₁ are split. This is where our "boundary supremum" wall lives. In Bałaban's approach it is controlled
  analytically, not by a probabilistic sup.
- **For gauge theory,** what plays the fill's role is Bałaban's *R operation*: a stability bound for the effective action in
  large-field regions (CMP 122, 1989, parts I–II).

**Conclusion.** The fill lemma (§§4–8) is a sufficient device, but Bałaban does not use it. The real (P1) input is an R-operation-type
stability bound. Since our model is a one-scale strong-coupling forest rather than a UV multiscale, the next concrete question is:
*does our cube coefficient admit a Bałaban-style split, with large-field cubes suppressed by their own action and small-field
Dirichlet dependence controlled analytically?* The calm-fill numerics remain useful as a sanity check of stability. Proved
status of every piece: unchanged (open).

**Job af45d8e1: invalid.** The covariant-phase solve failed (flux_err ≈ π, bd_max ≈ 2.8). The run was cancelled and **not fixed**, because
§8 settles that collar analytically and the Bałaban route makes the fill lemma optional. Dead end as numerics, dropped.

## 9. Bałaban-style split of Conjecture L: reduction and adversarial check (2026-11-05)

**Observation.** The covariance in Conjecture L, λ·Cov(F, ∂_jQ·1_{m≥ρ/2}), is taken in the **fibre measure π₀(c)**. That measure
conditions on V♯ for *every* block. So its small-field part is a *block-constrained* Gaussian. Its decay therefore comes from
the constraint (a Bałaban fluctuation propagator), not from a mass gap. This is the structure Bałaban proves at one scale:
- Bałaban, *Propagators and renormalization transformations for lattice gauge theories I*, CMP 95 (1984) 17–40;
- *Propagators for lattice gauge theories in a background field*, CMP 99 (1985). It proves regularity and exponential decay of the
  constrained propagators around small, non-flat backgrounds.
- Modern exposition of the decay mechanism (Combes–Thomas): Dybalski–Stottmeister–Tanimoto, arXiv:2303.10754.

**Proposed split** (sketch, not proved):
1. **Small field** (every coarse plaquette angle ≤ θ₀ < π, and fine field calm). The constrained Hessian has a uniform gap.
   Our heuristic data (`H1_AVERAGED.md`, gaugegap v2) give gap ≈ (2–2.4)/b²·(1−θ/π), even for uniform-field backgrounds.
   The constrained propagator should then decay like e^{−γd/b} (Bałaban 1985 for the gauge-covariant version), and the
   covariance with a local event decays by a random-walk/cluster expansion.
2. **Large field** (coarse angle > θ₀, or fine m ≥ ρ/2). Each such cube is suppressed by its own action, about e^{−cλθ₀²},
   and is summed in a polymer expansion. Pointwise in c this fails: c may itself be large-field, which is the uniform-field
   counterexample of `CALM_LINK.md` §4. **Averaged over c** (Lemma `mixture`), large coarse fields are suppressed by the
   marginal, which is the averaged (R1).

**Adversarial check.**
- *Degeneration at θ→π.* It is real: the gap ∝ (1−θ/π) → 0. So θ₀ must stay a fixed distance below π, and the large-field
  activity e^{−cλθ₀²} has to beat the polymer entropy. At λ = O(β) this is fine. It is **not** uniform in the coarse field,
  so the split works only in averaged form. That is consistent with the existing wall diagnosis, not a way round it.
- *Bałaban's small-field condition is perturbative* (|1−U_p| ≤ ε₀, with ε₀ small). Our calm region m < ρ/2 with ρ_G = 1 is
  not small in that sense. Either shrink ρ (this costs the constants in `SHELL_TRANSPORT.md` §1 a factor of 1/ρ) or prove the
  gap directly for the calm region. The numerics only cover uniform-field backgrounds.
- *Gauge fixing.* Bałaban's propagators use axial/Landau-type gauge inside blocks. Our transport Y is gauge-covariant, and
  gaugegap v2 removes the gauge directions exactly. No Gribov issue arises at small field, but at the calm/large interface
  this is unchecked.
- *What it does not give.* It gives decay at **one block scale**, which is what the cube coefficient needs. It gives **no**
  mass gap and does not touch the continuum limit.

**Status.** Conjecture L ⇐ (constrained gap on calm fibres, uniform for coarse angle ≤ θ₀) + (Bałaban 1985 decay) +
(polymer sum of large-field cubes under the averaged measure). Each item is **open** for our specific constraint map V♯. The
first item is being tested numerically by job gfdecay-l16 (constrained propagator profile, L = 16, b = 2, n = 0, 4, 12),
which is queued because the cluster is at its CPU budget.

### 9.1 A non-perturbative calm-fibre gap via Weitzenböck (sketch, continuum heuristic)

This is aimed at the "Bałaban's small field is perturbative" objection. Write the Yang–Mills Hessian at a background A on
gauge-orthogonal fluctuations a as ⟨a, (D_A*D_A) a⟩ + 2⟨a, ad(F)·a⟩ (Weitzenböck/Bochner). The curvature term is ≥ −2|F|·|a|².
On each block of side b, with the block constraint removing the lowest mode, a covariant Poincaré inequality gives
D*D ≥ c₁/b² − O(|F|). The holonomy around a block is θ = |F|b², so

  gap ≥ (c₁ − c₂θ)/b²,

which is positive for θ < θ₀ := c₁/c₂. The bound is linear in θ, like the measured law (2–2.4)/b²·(1−θ/π) in `H1_AVERAGED.md`.
The Weitzenböck constant is cruder, though: it puts the zero near θ ≈ 1 rather than π. So the calm condition is only
**|F|b² ≤ θ₀ for a fixed θ₀ ≈ 1**, which is not small in Bałaban's ε₀ sense. Large-field cubes then mean θ > θ₀, and their
suppression is e^{−cλθ₀²} with an O(1) θ₀.

**Caveats.**
- This is a continuum argument. The lattice version needs a lattice Weitzenböck formula with O(|F|²) corrections, which are
  harmless at small θ.
- Covariant Poincaré with zero *gauge-covariant* block average needs a proof. At θ ≪ 1, parallel-transport to the block
  corner gives it, with error O(θ).
- Status: **sketch.** It is consistent with the numerics, which also show the true threshold (π) is better than this bound.

### 9.2 Lattice version of 9.1 (sketch, with the steps made explicit)

Setting: a block B of side b, with every plaquette deviation |1−U_p| ≤ ε, and θ := b²ε.

1. **Axial (comb) gauge from the block corner.** Tree links become 1. Every other link equals the holonomy of a lattice
   strip of at most b plaquettes, so |U_l − 1| ≤ bε = θ/b. This holds nonabelian, because a loop's holonomy deviation is
   at most the sum of the deviations of the conjugated plaquettes of a spanning surface. **Proved** (standard).
2. **Covariant Poincaré.** Take ψ with zero covariant block average, which is the axial-gauge flat average.
   - Then ‖D_Uψ‖ ≥ ‖∇ψ‖ − (θ/b)‖ψ‖.
   - Flat Poincaré gives ‖∇ψ‖ ≥ √c₁‖ψ‖/b.
   - So ‖D_Uψ‖² ≥ (√c₁ − θ)²‖ψ‖²/b².
   **Proved** modulo the flat constant c₁ for the block map used: corner-type maps, not exactly a mean.
3. **Curvature term.** The Wilson Hessian equals Σ_p |(D a)_p|² plus terms bounded by C·Σ_p|1−U_p|·|a|² over the 6
   plaquettes per link. That is ≥ −6Cε‖a‖² = −6Cθ‖a‖²/b².
4. **Coupling the constrained fluctuations to step 2.** This passes through the Hodge split (gauge directions removed exactly,
   with D*D + DD* on 1-forms). It is the **open** step: a 1-form version of steps 1–2 is needed.

**Consequence (modulo step 4):** gap ≥ ((√c₁−θ)² − 6Cθ)/b² > 0 for θ ≤ θ₀ = θ₀(c₁, C) = O(1). That is the calm condition of
§9, now at lattice level with explicit dependence on the constants.

**Adversarial check.** Step 1 is the place where a naive version fails. Bounding links by the *whole-block* flux θ, instead
of by a strip of bε, gives the condition θ ≲ 1/b, which is useless at large b. The strip bound is what makes the condition
scale-free.

### 9.3 Step 4 of 9.2: reduction and a hidden snag (2026-11-05)

- **Inside one block.** Pick the comb tree so that it contains the corner lines of the constraint. Then J at U equals J
  at 1 exactly, and the kernel ker J does not move. In comb gauge H_U − H_1 is relatively form-bounded: for any δ > 0,
  H_U ≥ (1−δ)H_1 − (C²θ²/δ + Cθ)/b². The block gap therefore follows from the **flat** constrained gap c₁, so step 4
  reduces to computing c₁.
- **Snag: links between blocks.** The comb gauge is per block. A plaquette that straddles two blocks mixes two gauges, and
  the crossing links are close to the *coarse* link, not to 1. A global axial gauge would give |U − 1| ≲ Lε, which is
  useless. This is exactly what Bałaban's **complete axial gauge** handles (CMP 95, 1984): the background is replaced by a
  coarse field that is smooth on scale b, and the fluctuation is measured relative to that. So step 4 = Bałaban's
  complete-axial-gauge reduction plus the flat constant c₁. **Open** for our V♯, but now a literature-shaped task rather
  than a new idea.
- **Computing c₁.** It is translation-invariant over the coarse lattice, so a Bloch decomposition reduces it to small
  per-momentum matrices. That can be certified later with interval arithmetic. Job `bloch-flat-gap` computes it:
  2D b = 2, 3, 4 to validate against gaugegap n = 0 (2.0–2.44), and **4D** b = 2, 3, which earlier runs never covered.

### 9.4 Literature for c₁ and step 4 (2026-11-05)

Source: Dimock, *Covariant axial gauge*, arXiv:1410.2798 (abelian, single scale, d = 3 explicitly; the paper says the methods
are "not specific to any dimension").
- **Prop. 6.** For the block-*averaged* map Q (an average over the L³ block with weights L⁻⁴), under QA = 0 and covariant axial
  gauge τA = 0, ‖dA‖² ≥ C‖A‖² with C independent of N and k. The proof adapts Bałaban's lemma 2.4 (CMP 95). This is the
  flat-c₁ statement, **but for averaged maps, not our corner map**.
- **Sec. 3.8.** The constrained operator is strictly positive with an exponentially decaying kernel, and by Bałaban's lemma so
  is its inverse. This is the small-field decay of §9, in the abelian single-scale case.
- **Gap for us.** `H1_AVERAGED.md` found that only corner/[CD]-type maps are viable for our fibres, which rules out
  averaged maps. So Prop. 6 does not transfer directly. The proof pattern (Bałaban lemma 2.4) should, but this is
  **unchecked**. Non-abelian 4D remains Bałaban's CMP 95/99 only.
- Torons: on the torus Dimock imposes Q•A = 0 to suppress them. Our fibre conditions on all of V♯, which includes the
  analogue, so no extra issue is expected.

### 9.5 c₁ > 0 for the corner map, any d (sketch, elementary)

This is flat and abelian, with each colour component treated separately. Gauge: in each block B pick a comb tree T_B that
contains the corner lines starting at B's corner, and set A = 0 on T_B. The only residual gauge left is constant on each
block, i.e. the coarse gauge.
- **Constraint.** In this gauge, the corner-line constraint fixes each coarse link to 0. The remaining coarse gauge is then
  global and constant, which on the torus is trivial (torons are excluded by fixing V♯).
- **Links inside a block that are not tree links.** Such an A_l equals ∮ over the loop tree–l–tree, which is the sum of dA
  over a comb strip of ≤ b plaquettes. Hence |A_l|² ≤ b·Σ_strip|dA|², and each plaquette lies in ≤ b strips.
- **Crossing links** (x ∈ B, x+e ∈ B′). The loop runs from x along T_B to the corner, then along the coarse line (sum 0),
  then along T_{B′} to x+e. It has length O(b) and bounds a surface of O(b²) plaquettes. Each plaquette lies in O(b)
  surfaces.
- **Result.** ‖A‖² ≤ C_d b³‖dA‖², so **c₁ ≥ 1/(C_d b)** in units where the gap is c₁/b², and c₁ > 0 for every fixed b.
  That is all the one-scale cube coefficient needs, since b is a fixed constant. The b³ comes from using a single crude
  surface per crossing link. Averaging over surfaces, as in Bałaban's lemma 2.4, should improve it to b². The numerics
  (2.0–2.44 in 2D) suggest b² is the truth.

**Adversarial check.**
1. A ∈ ker J must be orthogonal to the gauge directions, not merely gauge-fixed. But the gap of the quotient form equals the
   minimum over gauge-fixed representatives only up to the norm distortion of the gauge-fixing projection. Axial-gauge
   representatives have norm ≥ that of the orthogonal representative? **No, the inequality goes the other way.** The
   orthogonal representative has the *smallest* norm, so ‖A_⊥‖ ≤ ‖A_axial‖ ≤ C b^{3/2}‖dA‖. The bound transfers. ✓
2. The corner-line links belong to the tree, but the last link of a corner line ends in B′. Assign the whole line to T_B
   and make T_{B′} avoid it. ✓ (bookkeeping only)
3. Status: **sketch, believed correct.** It removes the "unchecked transfer" item of §9.4 for fixed b. The Bloch job checks
   the constant.

### 9.6 Bypassing the complete axial gauge: IMS localisation + Combes–Thomas (sketch)

The §9.3 snag (between blocks, gauges don't match) matters only if one wants a *global* gauge. A gap is a local property:
- **IMS localisation.** Take a partition of unity Σχ_i² = 1 made of functions that are **constant on each block** and
  vary on a scale of R blocks, with jumps ≤ 1/R across block faces.
  - Multiplying by χ_i preserves ker J, because the constraint is linear and block-local.
  - H has range 1, so H ≥ Σ_i χ_i H χ_i − (C/R²).
  - On the support of each χ_i (diameter ≈ Rb), the axial gauge gives |U − 1| ≤ Rb·ε = Rθ/b. Combined with §9.5 and the
    relative bound of §9.3, this gives χ_iHχ_i ≥ (c₁ − C R θ)/b²·χ_i².
  - Result: gap ≥ (c₁ − CRθ)/b² − C/R². With b fixed, take R = R(b) large, then θ ≤ θ₀(b) > 0.
- **Combes–Thomas.** A finite-range operator with gap g has an inverse that decays like e^{−c√g·d}. So the **small-field
  decay of §9 follows from the gap, with no Bałaban propagator theory needed** at one scale.
- **Gauge directions.** Work with H + D(D*D)⁻¹... or more simply with the gauge-fixed form H + DD*-type term, which is also
  finite range. Then check that IMS commutes with it up to C/R². **Unchecked detail.**

**Adversarial check.**
- θ₀ now depends on b, whereas §9.1 hoped for a b-free θ₀. That is fine for the one-scale cube coefficient, where b is fixed,
  but it means the large-field suppression is only e^{−cλθ₀(b)²}.
- IMS needs χ to be block-constant to keep ker J. Then |∇χ| jumps at block faces, and the error is C/R² in fine-lattice units,
  not C/(Rb)². Accounted for above.
- Combes–Thomas on a *constrained subspace*: the inverse is P(PHP)⁻¹P, with P the projector onto ker J ∩ (gauge)^⊥. P is not
  finite range (the gauge projector is nonlocal). Using the gauge-fixed operator H + κ D D* restricted to ker J instead keeps
  locality, because ker J is block-local. This must be done that way. **Item to verify.**
- Status: **sketch.** It replaces "complete axial gauge (open)" by two routine-looking checks.

### 9.7 The §9.6 items (sketch, checked on paper)

**Gauge-fixed operator.** Let C₀ = {gauge parameters c vanishing on the corner lattice}, Π the projection onto C₀
(it zeroes the corner sites), and D the linearised gauge action. Then D(C₀) ⊂ ker J, because the corner map only
conjugates by c at block corners. Put H_κ = H + κ DΠD*. This operator is finite range (range 2) and preserves
ker J. Write ψ ∈ ker J as ψ = ψ_⊥ + DΠc with ΠD*ψ_⊥ = 0. If the background is a critical point, HD = 0, the cross
terms vanish, and

  ⟨ψ,H_κψ⟩ = ⟨ψ_⊥,Hψ_⊥⟩ + κ‖ΠD*DΠc‖² ≥ gap_phys‖ψ_⊥‖² + κλ₁^D‖DΠc‖²,

where λ₁^D ≳ 1/b² is the Dirichlet Laplacian gap with the corners pinned. So gap(H_κ|ker J) =
min(gap_phys, κλ₁^D), and κ = O(1) gives a gap of order 1/b². This settles the first item.
Caveat: at a non-critical background HD = O(|∇S|), which adds an error term of size small-field·ε.

**IMS.** Take χ_i constant on every corner line, which a block-constant χ is. Then J(χ_iψ) = χ_iJψ = 0, so the
localised pieces stay in ker J. The double commutator Σ[χ_i,[χ_i,H_κ]] is supported within range 2 of the
jumps of χ, and its norm is ≤ C‖H_κ‖/R², which is the claimed −C/R². This settles the second item.

**What remains.** The local gap on an R-block patch with free boundary. §9.5 is per-strip, so it localises, but it
has not been written for patches. The large-field polymer half is untouched. None of this has been tested
numerically yet; gfdecay is still queued.

### 9.8 Patch-local gap: reduction (sketch)

The §9.5 inequality ‖A‖² ≤ C_d b³‖dA‖² is assembled one block at a time. Each strip or surface sum uses only
the links of one block and that block's corner lines. On an R-block patch with free boundary it therefore
holds with the same constant, **provided** every link sits in the gauge of a single block. Links that cross
between blocks are exactly the §9.3 snag: the complete-axial-gauge issue. So after §9.6–9.7, the whole
small-field half reduces to one lemma.

**Lemma X (open).** In tree + corner-line gauge, the links crossing between blocks satisfy
|A_l|² ≤ C b^k (‖dA‖² on the two adjacent blocks) for some k.

For the corner map, a crossing link together with the two corner lines closes a plaquette chain of length O(b).
That suggests k = 3, matching §9.5. Lemma X has not been checked. Also note that b³ instead of b² still gives
θ₀(b) > 0, only with worse b-dependence, which is harmless at a fixed scale.

### 9.9 Lemma X: abelian 2D proof sketch, k = 3

**Setup.** Blocks B₁ = [bx, bx+b) × [by, by+b) and B₂ is its right neighbour. The crossing link l is horizontal
in row y. Inside each block use the axial gauge:
- the vertical links on the block's corner column are 0 below the top link;
- the top link is then 0 too, because the coarse-link constraint (sum along the corner line = 0) forces it;
- the horizontal links in each row are 0 except one per row (row y keeps only l).
This is a legitimate gauge on the patch: the block gauge functions match along the shared corner column
because both are fixed by the corner value.

**Stokes.** Take the rectangle with corners (bx, y), (bx+b, y), (bx+b, by+b), (bx, by+b).
- Its vertical sides lie on corner columns, so they are 0 in this gauge.
- Its top side is the corner line of the block above. Its sum is 0 by the constraint, and in that block's gauge
  its links are individually 0 except one, which is controlled by §9.5 for that block.
- Its bottom side is row y, whose sum is A_l.
Hence A_l = Σ_{rectangle} dA + (§9.5-controlled terms). Cauchy–Schwarz over an area ≤ b² gives

  |A_l|² ≤ 2b²‖dA‖²_{B₁∪B₂∪above} + 2C b³‖dA‖²_{local},

so k = 3, matching §9.5. **Status: sketch**, for abelian 2D only.

**Nonabelian, small field.** The same loop with the parallel transport inserted gives corrections of size O(θ)
relative to the main term (§9.2). It is not yet written in d dimensions, where the loop becomes a b × b strip
in each transverse plane; the argument is the same plane by plane.

**Adversarial check.** My first attempt took the loop with vertical sides on *partial* corner segments and
assumed those segments were constrained. They are not; only the full corner lines are. The fix is the gauge
choice above, which makes the partial segments vanish identically. Still to double-check: whether the axial
gauge in B₂ can fix B₂'s corner column given that B₁'s gauge already fixed the shared corner site. I believe it
can, since each block's gauge is free on its own interior sites, but this is unverified.
