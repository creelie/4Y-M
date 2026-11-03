# Slice ratio checks (open input of Proposition "slice ratio")

`composite_free.py` computes the slice correlations D(r) for composite, gauge-invariant slice fields of the free
lattice gauge field. Its output is in `composite_free.out`, and the results are Remark "composite free slice" in `../main.tex`.

- **Magnetic energy.** This is the free analogue of the paper's slice field, a function of the spatial links of one slice.
  - Ratios at fixed r > 0 converge to the continuum as a → 0, with or without the spatial flow.
  - Box side 6, t_f = 0.05: D(½)/D(¼) = 0.27341, 0.27239, 0.27201, 0.27174 at 1/a = 8, 12, 16, 24, against 0.27152 in the continuum.
  - Unflowed, the ratio tends to (r/r′)⁵ = 1/32.
  - With the flow the first-step energy −log(D(a)/D(0))/a stays bounded: 6.37, 6.50, 6.56, 6.62.
  - Without the flow it diverges, since D(a)/D(0) = 0.0198 at every spacing. So the r > 0 condition is strictly weaker than a first-step bound.
- **Full energy density (E + B).**
  - The electric part has a contact term in time, so D(0) ~ a⁻² even with the flow, and D(a)/D(0) = 0.042, 0.023, 0.0145 tends to 0.
  - The r > 0 ratios still converge to the same continuum values.

The next-to-leading-order leakage analysis (one-loop logarithm of the first-step energy, coefficient 8/11 of 2b₀,
UV-finite norm of the leaking component) is summarised in the last paragraph of the same remark. The full analysis
and its Monte Carlo checks are in `big-bang-monotile/yang-mills-continuum/LEAKAGE.md`.

    python3 composite_free.py      # about 20 minutes on one CPU core
