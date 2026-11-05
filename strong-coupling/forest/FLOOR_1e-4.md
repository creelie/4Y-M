# Where the 10⁻⁴ eigenvalue floor comes from

Status: **likely explanation, not proved.** Revised 2026-11-05 after an adversarial recheck. The first version overstated this.

## Source
`dump2.txt` ~l. 12176–12576, namespace `YMOpenGluonCert`, has 16 theorems `romanov_G_0000` … `romanov_G_1111`.
Each certifies (1/16384)|x|² ≤ xᵀAx for an 8×8 rational A, using an LDLᵀ/SOS certificate.

## Finding
Every one of the 16 matrices has the same spectrum, to the printed precision:

  {1.00·10⁻⁴ (×2), 1.838 (×3), 5.514, 204.9 (×2)}.

The two smallest eigenvectors are antisymmetric, v_{k+4} = −v_k. **Exact recheck (rational arithmetic, `floor_exact.py`).**
- A − 10⁻⁴·I has full rank 8 for all 16 matrices. So the data are not *exactly* "singular matrix + 10⁻⁴·I".
- λ_min is 9.9711·10⁻⁵ (G_0000) and 9.99993·10⁻⁵ (G_1011). These are within 3·10⁻⁷ of 10⁻⁴.
- The entries have denominators of at most 10⁶, so they are rounded to 6 decimals. Rounding a "rank-6 + 10⁻⁴·I" matrix would shift λ_min by about 10⁻⁷ to 10⁻⁶.
- So the data are *consistent with* that construction (a near-null mode plus a ridge of 10⁻⁴). An independent near-null eigenvalue that happens to sit at 10⁻⁴ in all 16 cases is the less likely alternative. This is still an inference.
- I have **not** identified the zero mode physically. The earlier "gauge/translation" label was a guess.
- The certified constant 1/16384 ≈ 0.61·λ_min is then the ridge, not a property of the model.

## Consequence
- Any downstream use of "λ_min ≥ 1/16384" as a physical coercivity/gap constant is **unsupported**. If the
  ridge is removed, the bound becomes 0 on a 2-dimensional subspace.
- The bound is legitimate only if the consumer either (a) works on the quotient by the antisymmetric zero modes
  (they are gauge-trivial there), or (b) adds the same ridge to its own functional.
- The consumer (OpenGluon calm-data note) is **not in the local repo**, so I could not check which of these it does.
  Before anyone relies on it: restate the cert on the 6-dim complement, with a floor read off the next eigenvalue (≈1.84; not certified), or document the ridge.
- Decisive check (not done): find the code that generated these matrices and see whether a ridge is added.
- The `exact4_diag_*` 3×3 certs (floors 1/4, 1/128) are not affected.
