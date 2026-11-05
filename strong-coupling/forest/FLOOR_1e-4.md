# Where the 10⁻⁴ eigenvalue floor comes from

Status: **diagnosed** (exact rational data re-read from `dump2.txt`; spectra computed in floating point, heuristic only).

## Source
`dump2.txt` ~l. 12176–12576, namespace `YMOpenGluonCert`, has 16 theorems `romanov_G_0000` … `romanov_G_1111`.
Each certifies (1/16384)|x|² ≤ xᵀAx for an 8×8 rational A, using an LDLᵀ/SOS certificate.

## Finding
Every one of the 16 matrices has the same spectrum, to the printed precision:

  {1.00·10⁻⁴ (×2), 1.838 (×3), 5.514, 204.9 (×2)}.

The two smallest eigenvectors are antisymmetric, v_{k+4} = −v_k. That is the signature of an exact zero mode,
a gauge/translation direction of the 4+4 split. So the underlying matrix is **singular (rank 6)**. The floor is
an added ridge of about 10⁻⁴·I, rounded to 6 digits. The certified constant 1/16384 ≈ 0.61·λ_min certifies
the ridge, not a property of the model.

## Consequence
- Any downstream use of "λ_min ≥ 1/16384" as a physical coercivity/gap constant is **unsupported**. If the
  ridge is removed, the bound becomes 0 on a 2-dimensional subspace.
- The bound is legitimate only if the consumer either (a) works on the quotient by the antisymmetric zero modes
  (they are gauge-trivial there), or (b) adds the same ridge to its own functional.
- The consumer (OpenGluon calm-data note) is **not in the local repo**, so I could not check which of these it does.
  Before anyone relies on it: restate the cert on the 6-dim complement, with an honest floor (≈1.838 before the
  ridge, i.e. something like 1.8 − O(10⁻⁶) once rounding is accounted for), or document the ridge.
- The `exact4_diag_*` 3×3 certs (floors 1/4, 1/128) are not affected.
