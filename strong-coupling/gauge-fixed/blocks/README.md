# Plaquette blocks for SU(2)

Checks for Section V of the forest-gauge paper (`strong-coupling/forest/main.tex`,
Theorem 3, "SU(2) with plaquette blocks"; numbers of the REVTeX version). A block is the set of the three free
links of one plaquette of the line-field gauge; the blocks used are invariant
under translations by 2Z^4.

| Program | What it does | Output |
|---|---|---|
| `cover.py` | finds the sixteen block plaquettes modulo two by an exact-cover search and checks Lemma 10 (partition, three distinct blocks in every other plaquette, ten outside partners, influence within distance 2) on the tori of side 4 and 6 | `cover.out` |
| `constants.py` | bounds, in ball arithmetic (python-flint, 120 bits), the Hessian eigenvalues and the norm of the third derivative of log Z0 on a ball | `constants.out` (`python3 constants.py 2.07 20000`) |
| `alpha.py` | the constants of the block bound and the column sum alpha = 40 beta K at given beta_W, read exactly | `alpha.out` (`python3 alpha.py 1/3 0.34 0.345 0.346`) |

Result: alpha(1/3) <= 0.95602 and alpha(0.345) <= 0.99729, so the block
Dobrushin condition holds for beta_W <= 0.345 and the transfer-matrix gap is at
least -(1/2) log alpha(beta_W); at 0.346 the bound exceeds one.

`float/` holds the floating-point search of Remark 16 (not certified, not used in
any proof): `tri.py` (quadrature on S^3 and the block law), `cert2.py` (a bound
for the block sensitivity that keeps the law of the middle link), `red.py` (the
six parameters left after the SU(2)^3 symmetry), `opt4.py` (random search and
Nelder-Mead, `python3 opt4.py <beta> <samples> <seed>` with beta = beta_W/4), and the
logs at beta_W = 0.34, 0.35, 0.36 (largest values 0.278284, 0.279506, 0.280754).
It needs numpy and scipy.
