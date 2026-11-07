# Progress since 3795815 (directions after Hexagon v7)

Branch `research/bad-regions-multiscale-continuum`, 8 commits. Nothing here is Lean-checked.

## A. Gaussian fibre model → Thm 114 (`forest/A_GAUSSIAN_FIBRE_4D.md`)
- The 4D data (section, block map, curvature operator, fibres, the three hypotheses (F), (X), (D)) are written out in one place.
- **Lemma A1:** (F) holds on all of Gaussian L² with the first-chaos constant, C_F = λ_F, with no loss.
- **Lemma A2:** (X) holds on all of L² with ρ = ρ_1, the first-chaos constant.
- **Lemma A3:** the (D) vectors have an exact chaos expansion, w = Σ_{S≠∅} B^{⊗S}⊗A^{⊗rest}⊗I, orthogonal in |S|.
  - **Open item A-D:** this reduces the uniform bound on (D) to one combinatorial bound c_k(s) on symmetrized contractions. The naive bound C(k,s)² is not uniform in k.
- The question for James about lean05 is drafted and not sent.

## B. Scale summability
- Waiting on the ρ_1 computation (job `ym-gaussrho-a1`, queued). ρ_1 as a function of n'/b decides whether Σ_j ρ_j² ≤ S is plausible.

## C. (H1) path-covariance campaign (`forest/H1_PATHCOV.md`)
- **Run 1 result (negative, on method):** the raw Σ|Cov| estimator measures sampling noise. It matches the time-permuted null to < 0.4% at β_W ≤ 1.5 and scales as 1/√(samples). At β_W = 2.3 there is an unexplained excess of about 10%.
- **Run 2 is running:**
  - a block estimator of Cov², unbiased, with jackknife errors and a per-shell Cauchy–Schwarz bound;
  - a synthetic planted-signal control.

## D. Keeping (E) alive (`forest/D_E_OBSERVABLE.md`)
- **Observable and normalization:** the 24-plaquette action density, with Z_k = 𝒩_k(ε)^{−1/2}.
- **Audit of every use of the gap (MIX) in §limitE, complete:**
  - The gap is **never a premise** of (E). It only yields conclusions: clustering and the continuum gap.
  - Nontriviality under this normalization is gap-free (Cor limitEtwo(b)). The first pass said otherwise and was corrected.
- **Remaining (E) inputs:**
  - (SD)⇔(UV) and (NE) for n = 2;
  - (M_s) for n ≥ 3;
  - (U_s) and (W_s).

## Next
- Process C run 2. Check the control first, then plan the depth scan, the V stress test and the bphidden control.
- Process ρ_1 and feed it into B.
- Attack c_k(s) for A-D.
