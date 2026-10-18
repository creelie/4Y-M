# Sharp curvature and finite-box criteria for lattice Yang–Mills theory

Deep Bhattacharjee and Ushashi Bhattacharya.

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
- **Theorem D.** SU(2): a unique infinite-volume Gibbs measure for
  |β| < 1/(8(d−1)), β_W < 1/6 at d = 4; single-link Wasserstein–Dobrushin
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
