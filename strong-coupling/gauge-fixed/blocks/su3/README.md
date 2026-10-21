# Plaquette blocks for SU(3)

Checks for Section VIII of the forest-gauge paper (`strong-coupling/forest/main.tex`,
Theorem 6, "SU(3) with plaquette blocks"; numbers of the REVTeX version). The
blocks are those of `../` (Lemma 11); the one-link quantities come from the
SU(3) cover in `../../checks/su3/deg2/`.

| Program | What it does | Output |
|---|---|---|
| `consts.py` | bounds, on every box of side r/8 that meets the reduced domain, the eigenvalues of the Hessian of log Z (the covariance of Re Tr(U h) under the one-link measure) and its third derivative, from the exact moments of degree at most three at the corners and the centre, with the interpolation of Lemma 28 (Lemma 32) | `r3_2/consts_boxes.txt`, `r3_2/consts_log.txt` (`python3 consts.py --r 3/2 --L 3 --out r3_2`) |
| `alpha3.py` | reads the cover log (Λ(r)) and the constants, and evaluates w, s, C, K_3 and α_3 = 10 Nβ K_3 of equation (15) in ball arithmetic at 120 bits, at given β_W read exactly | `r3_2/alpha3_out.txt` (`python3 alpha3.py ../../checks/su3/deg2/r3_2 r3_2 0.72 3/4`) |

Result at r = 3/2: 841 boxes, eigenvalues in [0.056015, 0.247316] and third
cumulant at most 0.182126; α_3(0.72) ≤ 0.91815 and α_3(3/4) ≤ 0.96521, so the
transfer-matrix gap at β_W = 3/4 is at least 0.0177.

Result at r = 38/25 (`python3 consts.py --r 38/25 --L 3 --out r38_25`, with the cover
`../../checks/su3/deg2/r38_25` from `python3 cover2.py --r 38/25 --thr 1.03 --out r38_25`,
1063 boxes, 2rΛ(r) ≤ 1.029967): 784 boxes, eigenvalues in [0.054978, 0.248280] and third
cumulant at most 0.183546; α_3(0.76) ≤ 0.99950 (`python3 alpha3.py
../../checks/su3/deg2/r38_25 r38_25 0.76`, output in `r38_25/alpha3_out.txt`), so the
blocks give a gap for β_W ≤ 0.76, at least 2.5·10⁻⁴ at β_W = 0.76.

`consts.py` runs four worker processes and sets one BLAS thread per worker;
with the default multithreaded BLAS the four workers oversubscribe four cores.
It needs numpy and python-flint, and the modules of `../../checks/su3/deg2/`.
