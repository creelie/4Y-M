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

## 10. Redirect: the large-field half (2026-11-05)

**Adversarial verdict on §9–9.9.** The small-field half is Bałaban's propagator theory (CMP 95, 99, 1984–85),
which is already proved in greater generality. §9.1–9.9 re-derive it at sketch level. From now on that half is
**cited, not re-proved**, and the bloch-flat-gap job is dropped. gfdecay stays: it is cheap, and it checks the
corner map directly.

### 10.1 Lead: averaged large-field density via chessboard estimates

Large-field regions are bad under a supremum over exteriors, but not under the torus measure. Averaged
over c, the joint measure is the torus Wilson measure, which is reflection positive (Osterwalder–Seiler).
Chessboard estimates (Fröhlich–Israel–Lieb–Simon) then give, for a set X of plaquettes of one orientation class,

  μ(|1 − U_p| > ε for all p ∈ X) ≤ z_ε^{|X|},
  z_ε = μ(the event at every translate)^{1/N}.

A disseminated bad event costs action ≥ βε²·N/2. The partition function satisfies Z ≥ e^{−C N log β} by a Gaussian
lower bound, so

  z_ε ≤ exp(−βε²/2 + C log β),

which is small once βε² ≫ log β. Summing over orientation classes costs a factor 2^{|X|}. The result is a
**Peierls bound for large-field sets in averaged form, with no exterior supremum**, which is exactly the form §3
said is needed. **Status: sketch.** The technique is standard for spin systems, and its use for lattice gauge theory
is probably known (Borgs; Fröhlich–Israel–Lieb–Simon), so the literature needs checking before this is claimed as new.

**What it does not give.** Conjecture L needs *decay* of Cov(F, ∂_jQ·1_{large near j}), not only small density.
- Density alone gives |Cov| ≤ 2‖F‖ μ(large near j), uniform in d.
- Decay needs either the small-field cluster expansion with the large-field set as a random polymer, which is
  Bałaban's R-operation, or a coupling argument.

**Next concrete test.** Does chessboard density plus Bałaban small-field decay plus an annulus argument give decay?
- Split by whether a large-field path connects j to supp F.
- Connected: probability ≤ (2z_ε)^{d/size} by chessboard, summed over paths.
- Not connected: condition on the large-field cluster of j. The obstacle is that the conditioning breaks reflection
  positivity, and the remaining measure is the fibre measure with a random exterior; that is the original wall
  again. **To check:** whether the conditioning can be avoided by working with the event "no large-field
  crossing of the annulus" and using small-field decay *inside* the annulus only.

### 10.2 Reduction attempt for averaged Conjecture L (sketch, with a flaw found and a repair route)

**Literature check (2026-11-05).** Chessboard and contour bounds in a gauge-theory setting do exist:
- arXiv:2501.10065 (Z₂ π-flux, J. Stat. Phys. 2026) bounds the probability of bad plaquette configurations;
- Biskup–Kotecký, math-ph/0610060, is the general reference for contour suppression via reflection positivity.

I found no explicit continuous-group large-field Peierls bound of the form in §10.1. Treat it as plausible but not
yet located: it may be new, or folklore.

**Choice of ε.** Take ε = β^{−1/4}. Then βε² = β^{1/2} ≫ log β, so chessboard suppression holds. At the same
time the coarse angle θ ≲ εb² → 0, so Bałaban's small-field regime holds. The two halves are compatible at weak coupling.

**Attempt v1 (wrong).** Let N = {no large-field path from j to supp F}, and let S be the outermost separating
small-field surface. S is a stopping set, so the domain Markov property holds for the fibre measure if S is
block-aligned. This gives Cov_c(F,G) ≤ C π_c(N^c) + ‖F‖·osc_{b_S} E[G | b_S], and the plan was to bound the
oscillation by small-field decay across the annulus. **Flaw:** the annulus between S and j's cluster is not
all small-field; it has a density z of large-field islands. A fully small-field annulus of radius d has probability
about (1 − z)^{d⁴} → 0. So v1 secretly assumes what Bałaban's R-operation provides.

**Repair v2: nested shells and a Doeblin condition.**
- On N, explore inward from F's side to get many disjoint small-field shells around j: k ≳ d/t shells of
  thickness t, except on an event whose μ-probability is bounded by chessboard plus path counting.
- The Markov property makes the shell configurations X₁, X₂, … a Markov chain going inward.
- If every transition satisfies a uniform minorisation P(X_{i+1} ∈ · | X_i = x) ≥ η ν(·), then coupling gives
  osc ≤ (1 − η)^k.

**Open Lemma D (the real content).** For a small-field shell of thickness t, with c fixed and θ < θ₀, the law of
the inner boundary of the shell given the outer data has total-variation oscillation ≤ 1 − η, uniformly.
- Heuristically, the Gaussian part shifts the mean by ≲ εe^{−αt} against fluctuations of size β^{−1/2}, with α
  the constraint mass from gfdecay.
- The TV distance also picks up a factor √(area). Hence t ≳ log β + log(area), which gives decay
  e^{−c d/log d}: subexponential but summable. This is enough for the cube coefficient only if it is used with
  polynomial losses; check this.

**What is genuinely new versus Bałaban.** The large-field sets are controlled by an averaged chessboard bound,
not by the R-operation, and the composition is a Markov-chain coupling across shells, not a polymer expansion.
If Lemma D holds, the one-scale large-field half closes in averaged form.

**Not addressed.** Scale composition, i.e. iterating over all scales with uniform constants, and the continuum
limit. Chessboard needs the torus measure at each scale; after one renormalisation step the effective action
is not obviously reflection positive. This is a probable obstruction to iterating.

### 10.3 Scale composition: does the chessboard bound iterate? (quick test, negative for now)

**Idea.** For a block-symmetric map, large-field events at scale k are events of the *original* torus measure,
measurable in cells of side L = b^k. So chessboard with cell-sized reflections applies at every scale, with no need
for the effective action to be reflection positive. Our corner map is not reflection symmetric; an averaged map
would be.

**Energy.** In 4D the classical action is scale invariant. A holonomy angle θ around the L × L loop of a cell costs
action ≳ θ², independent of L, so the disseminated event costs e^{−cβθ²} per cell.

**Entropy, where it fails.** The crude bound Z ≥ e^{−C L⁴ log β} per cell swamps the energy cost once
L⁴ log β ≳ β. A sharper Peierls bound needs an energy–entropy comparison map at scale L. That comparison is
exactly the running-coupling (asymptotic freedom) problem, β_k ≈ β − b₀ k log b. **Verdict:** chessboard does not
bypass scale composition; it relocates it into the per-cell large-field probability at scale L. Iterating is
not available cheaply, so the effort stays at one scale (Lemma D).

### 10.4 Lemma D, Gaussian part (sketch) and adversarial pass

**Gaussian part.** Expand around the background minimiser with fixed outer data x (Bałaban's background field).
The inner-boundary law is N(Mx, Σ), with Σ independent of x. Then:

- TV(N(Mx,Σ), N(Mx′,Σ)) ≤ ‖Σ^{−1/2}M(x − x′)‖.
- **Bound on Σ^{−1/2}.** Σ ≥ (β‖H_κ‖)^{−1}, so ‖Σ^{−1/2}‖ ≤ (Cβ)^{1/2}. No gap is needed for this bound;
  the gap enters only through M.
- **Bound on M.** M = G_f(inner, shell)·H(shell, outer) satisfies ‖M‖ ≤ Ce^{−αt} by Combes–Thomas from the
  §9.6 gap, which is Bałaban's.
- **Bound on x − x′.** With c fixed, two small-field data sets differ by ≲ εb per link, so
  ‖x − x′‖ ≲ εb√A, where A is the shell area.
- **Result:** TV ≲ β^{1/2}·εb·e^{−αt}√A = b·β^{1/4}e^{−αt}√A. So TV ≤ 1/2 once
  t ≥ α^{−1}(¼ log β + ½ log A + log 2b). This matches §10.2.

**Nonlinear part (sketch).** Fluctuations around the background have size β^{−1/2}, so the cubic terms are
β·β^{−3/2} = β^{−1/2}, which is small. The background itself depends nonlinearly on x, and its decay into the
interior is Bałaban's background-field estimate (CMP 1984–85, proved). Passing from KL to TV via Pinsker plus a
perturbative KL bound should finish it.

**Adversarial pass.**
1. *Do small-field shells exist?* I worried this repeats the v1 flaw: a full shell of area r³ is entirely
   small-field with probability (1 − z)^{r³ t} → 0. It does not repeat it, because the shells need not be round:
   they can be block-aligned separating sets that bend around the bad islands. Such sets exist unless a bad
   cluster of diameter ≳ t appears. By chessboard and path counting, that has probability ≤ r⁴(Cz)^{t/b},
   which is negligible. So the shells exist w.h.p., but they are irregular. Lemma D is therefore needed
   uniformly over irregular block domains; §9.8 localisation is per block, so this looks fine.
2. *Data adjacent to large fields.* The inner and outer data may sit next to large-field islands. Lemma D
   conditions on them, and the Gaussian bound used only the small-field property of the shell itself. OK.
3. *Uniform in c?* Only θ < θ₀ on the shell is needed. If the coarse angle is large somewhere, that is a
   coarse-level large field, which belongs to scale composition and is not handled here. **This is a genuine gap.**
4. *Does subexponential decay suffice?* e^{−cd/log d} times polynomial counting factors is summable over the
   cube coefficient's sum over d. It still has to be checked against the specific sum in the cube coefficient.

**Numerical test (to write).** Compute ‖Σ^{−1/2}M‖ exactly for the constrained abelian field on a
shell geometry as a function of t, using gfdecay's assembly code. The prediction is a slope −α equal to the
gfdecay rate. This test is cheap: 2D, 4 CPUs.

### 10.5 Corrections to the §10.4 adversarial pass

**Item 3 (coarse angle) is not a gap; it was self-inflicted.**
- Inside a small-field shell, every plaquette satisfies |1 − U_p| ≤ ε.
- Composing the transported plaquettes gives a block holonomy θ ≤ Cb²ε. This bound holds for nonabelian groups
  too, for small ε.
- With ε = β^{−1/4}, this is below θ₀ once β ≫ (Cb²/θ₀)⁴.
- So θ < θ₀ holds automatically on every shell. Coarse-large blocks are just fine-large blocks, and they are
  already inside the chessboard-controlled bad set.
- §10.3's entropy problem appears only when the cell side L grows with scale. At fixed b the chessboard ratio is
  exp(−cβε² + C log β), with no factor of b⁴.

**Item 4 (summability): fine.**
- The cube coefficient sums Ξ_j over j (`GAUSS_GREEN.md` §5).
- At distance d there are ≲ d³ observables. Hence Σ_d d³ e^{−cd/log d} < ∞, and the prefactor carries the
  chessboard density.
- So subexponential decay suffices for the cube coefficient. Whether the *downstream* polymer expansion (`OPEN.md`,
  "compose") needs exponential weights e^{a|X|} has **not** been checked.

**Status of averaged Conjecture L at one scale.**
- Chessboard density (sketch; literature check incomplete).
- Shell existence (sketch).
- Lemma D: Gaussian part sketched, shelltv job running; nonlinear part by citing Bałaban.
- Markov-chain coupling across shells (standard).
- No ingredient is proved here. The wall is now **Lemma D's nonlinear part**, plus the unchanged multi-scale and
  continuum issues.

### 10.6 Downstream check: what D(n, α) actually needs (adversarial)

`OPEN.md` l. 55 states the consumer: the averaged-cube condition D(n, α) with **n·a ≤ ℓ uniformly along the
continuum tuning**, so n ~ ℓ/a → ∞ in lattice units. §10 delivers Conjecture L only at a *fixed* block size b. To
reach n one must compose about log_b n scales, and that is §10.3, which is blocked.

**Shortcut tested: take one step with b = n (dead end, quantitative).**
- Small field on the shells requires θ ≤ Cn²ε < θ₀, so ε ≲ θ₀/n².
- Chessboard requires βε² ≫ log β, hence **β ≫ n⁴ log β**.
- Along the tuning, asymptotic freedom gives n ~ ℓ/a ~ e^{cβ}, so β ~ log n ≪ n⁴.
- So the single big step fails, and it fails for the physical reason: at the scale of the correlation length the
  coarse field is O(1), not small.
- **No version of the one-scale argument avoids running-coupling control across ~β scales.** This is the Clay
  core, and it is unchanged.

**Net.** §10 gives a plausible *new route for the one-scale large-field step*: averaged chessboard plus shell
coupling instead of Bałaban's R-operation. It does nothing for composition or for the continuum limit.
A one-scale step is already within Bałaban's (proved) framework. So its value is a possibly simpler proof of a
known-type step, **not progress on the Clay problem itself**.

### §10.7 Lemma D, nonlinear part: reduction to Bałaban small-field analyticity (SKETCH, not verified)

Claim. Suppose a shell Σ of thickness t and cross-section A is small-field everywhere. Then for any two outer
boundary data ω, ω', the inner conditional laws satisfy
  ‖P_ω − P_ω'‖_TV ≤ 1 − δ with a fixed δ > 0,  provided t ≥ C log(A β).

Idea.
- Doeblin needs only a constant overlap, not closeness. So it suffices to bound the Radon–Nikodym derivative
  dP_ω/dP_ω' on an event of probability ≥ 1/2.
- In the small-field region, Bałaban's effective action (CMP 95/99, cited) is analytic in the fluctuation.
  Its kernels decay like e^{-α|x−y|}.
- So ω enters the inner region only through couplings that cross the shell. Those are O(e^{-αt}) per site.
- This gives |log dP_ω/dP_ω'| ≤ C A β^{1/4} e^{-αt}, which is the Gaussian §10.4 bound with the
  cubic/quartic terms absorbed into the analytic kernels. Taking t ≥ C log(Aβ) makes it O(1).

Adversarial check.
(i) The hypothesis "shell is small-field" is itself random. Its failure probability is ≤ A t z_ε by §10.1. This
    is small only when β ≫ log(A t), so the argument works at a fixed scale only. It does not touch the §10.6
    wall, and the result is still one-scale.
(ii) The bound on the log-ratio must hold uniformly, not just in expectation. Bałaban gives deterministic
     analyticity bounds on the small-field domain, so this looks fine. I have not checked his domain conventions
     against our nested shells.
(iii) Gauge fixing: the analyticity is stated in axial/Landau-type gauges on blocks. Our shell needs gauge fixing
      that is compatible across the inner boundary. Same issue as §9.7, where block-constant χ handled it. Not
      checked here.
Status: Lemma D, small-field part, is plausibly reducible to cited results at one scale. The scale-composition
wall is unchanged.
