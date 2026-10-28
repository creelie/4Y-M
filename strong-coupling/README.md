# Mass gap of lattice Yang–Mills theory at strong coupling

Deep Bhattacharjee and Ushashi Bhattacharya.

## The current paper (`prd/`)

`prd/main.tex` is the paper "Mass Gap of Lattice Yang–Mills Theory up to the
Dobrushin Ceiling", in REVTeX 4-2 (Physical Review D, two columns, 29 pages),
compiled as `prd/YangMillsDobrushinCeiling.pdf` (also at
`preprintYangMillsDobrushinCeiling/YangMillsDobrushinCeiling.pdf` in the
repository root). Build with `latexmk -pdf -bibtex- main.tex`; the
bibliography is written by hand and every DOI was checked. The figures are
200 dpi PNG files written by `prd/figures/make_figures.py` from the data in
`prd/figures/data/`.

What it proves, all by hand and without conditions:

- **Theorem A.** For SU(N), SO(N) and U(N) in dimension d the Hessian of the
  Wilson action is bounded by 4d, attained with both signs.
- **Theorem B.** Bakry–Émery, log-Sobolev and Poincaré inequalities uniform in
  the volume for |β| < 1/(8d) (SU(N)); β_W < 1/8 for SU(2) at d = 4.
- **Theorem C.** No curvature argument uniform over configurations goes
  further.
- **Theorem D.** A weighted coupling distance removes the factor two of Shen,
  Zhu and Zhu, so uniqueness and decay hold on the whole curvature range.
- **The sharp quarter.** For every von Mises–Fisher law on the three-sphere
  and every unit vector h, the covariance of a chordally 1-Lipschitz function
  with ⟨h,x⟩ is at most 1/4; equivalently W₁(ν_b, ν_b′) ≤ |b − b′|/4, with
  equality at b = 0. Proof: McShane extension and Green's formula in the ball,
  a spectral split into zonal and transverse parts with closed Bessel sums,
  positivity of every Taylor coefficient of one entire function, and a
  Riccati comparison.
- **Theorem E.** Four-dimensional SU(2) has a unique infinite-volume Gibbs
  measure for β_W < 2/9, which is the ceiling of every single-link
  Wasserstein–Dobrushin argument.
- **Theorem F.** In the same ranges the lattice theory has a positive
  transfer-matrix gap, m ≥ −log(9β_W/2) for SU(2).
- **Theorem G.** A gauge-invariant Dobrushin–Shlosman finite-box criterion,
  stable in every 't Hooft twist sector; two-dimensional SU(2) meets it at
  every coupling.

It does not prove the Clay mass gap. Section VIII lists the gaps: the
intermediate fixed coupling, the crossover between strong coupling and
Bałaban's weak-coupling control, the continuum limit, and a gap in physical
units remain open. `OPEN.md` records what is known of each.

## Forest gauge fixing (`forest/`, `gauge-fixed/`)

`forest/main.tex` is the paper "Mass Gap of Four-Dimensional Lattice
Yang–Mills Theory for Every Compact Simple Gauge Group" by Deep
Bhattacharjee, Ushashi Bhattacharya and James Montgomery Lockwood (REVTeX
4-2, two columns, 38 pages), compiled as
`preprintYangMillsForestGauge/YangMillsForestGauge.pdf` in the repository
root; `forest/cmp/main.tex` is the same paper in the AMS article class (58
pages). Its figures are written by `forest/figures/make_figures.py`. The
paper is self-contained: the transfer-matrix reconstruction and the sharp
quarter for the single-link law of SU(2) are proved in its appendices.
The working notes and the scripts behind it are in `gauge-fixed/`.

Fixing the gauge on one link of every plaquette, along a field of straight
lines, lowers the number of Dobrushin partners of a link from 18 to 12. Three
results follow, proved by hand and not yet refereed (see
`gauge-fixed/README.md`):

- four-dimensional SU(2): uniqueness of the free-boundary state, decay and a
  transfer-matrix gap m ≥ −log(3β_W) for β_W < 1/3, the ceiling of
  single-link Dobrushin after any forest gauge fixing;
- every SU(N): the gauge-fixed Hessian is bounded by 29/2 (against 16), so
  Bakry–Émery, decay and a gap hold for |β| < 1/29 (SU(3): β_W < 9/29);
- four-dimensional SU(3): single-link Dobrushin in the same gauge, with an
  analytic lower bound for the spectral gap of the one-link measure and a
  corrector of degree two evaluated with the exact moments of the one-link
  measure at every staple sum, gives uniqueness, decay and a gap for
  β_W ≤ 3/4; blocks of the three free links of a plaquette extend it to
  β_W ≤ 0.76;
- every compact simple group: with κ_G = h^∨/(2ℓ_ρ) the least Ricci
  curvature, uniqueness, clustering and a gap m_G(β) for N|β| < 2κ_G/29;
- conditional: a Dobrushin–Shlosman condition on one cube of fixed physical
  size, uniform along the continuum tuning, gives a mass gap for every
  continuum limit; the cube estimate and the existence of the limit are open.
- what (D) requires: it holds at strong coupling, forces the gauge-fixed links
  of a cube with trivial exterior to stay a fixed mean distance from the
  identity (so it fails at every fixed cube size at large β), and cannot come
  from single-line bounds uniform in the exterior (an SU(2) line with gap
  e^(−b/1280) up to powers of b); for SU(2) it needs
  n⁴ ≥ 2(1−α) β_W/(5 L(β_W)), and the n² lost in that bound is of the right
  order (a stationary flux configuration);
- the quadratic approximation of the forest gauge: the Gaussian field of the
  Wilson action with the forest links fixed has class sums of cube
  coefficients at least the cube size, so (D) fails there at every n;
- a Gaussian benchmark: the averaged cube coefficient of a Gaussian field is
  exactly 1 when it is massless and bounded below 1 at a fixed physical scale
  when it has a mass;
- smoothed observables: the continuum step holds after the spatial gradient
  flow, which keeps reflection positivity and the transfer matrix.

Checks in `prd/checks/` (the paper does not rely on them):

| Script | Checks |
|---|---|
| `pairblock_derive.py` | the pair block M of Lemma 2 from the trace expansion |
| `sos_delta.py` | the reflection M(−1,θ) = M(1,θ+π), the real form Δ*MΔ and its sum of squares |
| `barta_riccati.py` | the Barta computation and the first two Riccati barriers of Appendix B |
| `quarter/check_quarter.py` | the coefficients of the function 𝒟, the Bessel sums, both parts of the norm below 1/4, the third barrier and the factorisation of the zonal majorant (49 checks) |

The Lean files `YangMills/StrongCoupling/{QxCoefficients,ZonalBound,Thresholds}.lean`
certify the finite and algebraic steps of the sharp quarter and of
Theorems E and F (25 theorems, see the root README);
`YangMills/StrongCoupling/GaugeFix.lean` certifies the combinatorics of the
forest gauge (7 theorems).

## The earlier amsart version (`paper/`)

`paper/` holds the source (`main.tex`, figures as TikZ sources and 200 dpi PNG)
and the compiled `YangMillsStrongCoupling.pdf` (17 pages). Build with
`latexmk -pdf -bibtex- main.tex`; the bibliography is written by hand.

## What the paper proves

- **Theorem A.** For SU(N), SO(N) and U(N) on any lattice of dimension d, the
  Hessian of the Wilson action is bounded by 4d (in place of 8(d−1)), and 4d is
  attained with both signs.
- **Theorem B.** Bakry–Émery, log-Sobolev and Poincaré inequalities uniform in
  the volume, and decay of correlations, for |β| < 1/(8d) (SU(N)) and
  |β| < (N−2)/(16dN) (SO(N)). At d = 4: β_W < 1/8 for SU(2), 9/32 for SU(3).
- **Theorem C.** No curvature argument uniform over configurations goes further.
- **Theorem D** (superseded by Theorem E of the PRD version). SU(2): a unique
  infinite-volume Gibbs measure for |β| < 1/(8(d−1)), β_W < 1/6 at d = 4; single-link Wasserstein–Dobrushin
  arguments stop at β_W = 2/9.
- **Theorem E.** A gauge-invariant Dobrushin–Shlosman finite-box criterion,
  stable in every 't Hooft twist sector, whose hereditary form gives the area
  law; two-dimensional SU(2) meets it at every coupling.
- Section 7: a gap on every finite torus, gap ⇔ clustering, and the reduction
  of a gap at fixed coupling to a bound uniform in the torus size.

It does not prove a mass gap. Section 8 states the one missing step (the
crossover between the strong-coupling regime and Bałaban's weak-coupling
control), which is open.

## Checks

`checks/` contains the independent checks of the hand proofs (Python with
NumPy, SciPy, SymPy; `opt_su2.py`, `opt_sun.py` also use PyTorch).

| Script | Checks |
|---|---|
| `sos.py`, `verify_sos.py`, `verify_pair_exact.py`, `mat.py`, `check_pair.py`, `pairblock.py` | Cartan and pair blocks of Lemma 3.2, the SOS form, exactly |
| `check_identity.py`, `verify_lemma.py`, `verify_blocks.py` | Lemma 3.1 (divergence identity) and Lemma 3.2 on random plaquettes |
| `hess.py N d L`, `dd.py N d L`, `climb.py N d`, `hess_su2.py`, `opt_su2.py`, `opt_sun.py` | Theorem A against the numerical Hessian on random, identity and frustrated fields, and a hill-climb on the top eigenvalue |
| `barta.py`, `vmf_w1.py` | Lemma 5.1 (Barta bound on S³) and the Wasserstein rates of Section 5 |
| `su2_2d.py` | the disk integration of Section 6.4 |
| `torus_gap_check.py` | Proposition 7.1 on small tori |
