# SU(3) at β_W = 3/4: the corrector of degree two

Programs of Section 6 of `../../../su3.md` and of Proposition 20 of the forest-gauge paper.

| File | Content |
|---|---|
| `ring6.py`, `tsmul.c`, `cmul.py` | truncated Taylor polynomials of degree six in 36 real directions; `cmul.py` compiles `tsmul.c` (gcc) on first use and replaces the product by the C kernel |
| `rs6z.py` | moments of degree at most six of the one-link measure at a diagonal staple sum: character expansion truncated at λ_1 − λ_3 ≤ 34, λ_3 ≤ 11, Cauchy bound for the tail, running error analysis for rounding, products in the subring of torus weight zero |
| `rseries6.py` | the same in the full ring (slower; used to cross-check `rs6z.py`) |
| `polys.py` | the exact polynomials φ_h, the fields X_a of su(3), Δ, ⟨∇φ_{F_k}, ·⟩ and the basis of 𝒫₂; `scaled(C)` returns the Gaussian-integer polynomials of the corrector and asserts that every partial sum stays below 2⁵⁰ |
| `box2.py` | the bound of Lemma F″ on a box: choice of C, rounding to 2⁻³⁰, LDLᵀ in ball arithmetic |
| `cover2.py` | the cover at r = 3/2 (`python3 cover2.py --r 3/2 --out r3_2`), output in `r3_2/` |
| `face2.py` | floating-point values of both correctors on the face diag(r, ρ, ρ), output `face2.txt` |

`r3_2/cover_boxes.txt` lists the 2341 certified boxes (lower and upper corners in (a_1, a_2, x_3, y_3), the weight θ of the choice of C in the column headed "kappa", and the bound on 2r sup_h L_h). Requirements: Python 3 with numpy and python-flint, and gcc.
