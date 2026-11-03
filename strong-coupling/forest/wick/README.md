# Finite-box interacting coefficients and nonlinear flow checks

These standalone mathematical programs reproduce the computations in
`../B_WICK.md` and `../NONLINEAR_FLOW_CHECK.md`. They contain no cluster
submission code or machine-specific configuration. Generated result files
are not versioned. Run each program in a separate temporary working
directory because each writes `result.json` in the current directory.

| Program | Inputs and dependencies | Scope |
|---|---|---|
| `exact.py` | Python standard library; no arguments | exact rational 3x2x2x2 Wick coefficients and equality in two maximal trees |
| `wick.py` | NumPy and JAX; `--device cpu --geometries 3x2` | double-precision full finite-box correction, with optional `--reference` quaternion/Wick diagnostics |
| `jet.py` | NumPy and JAX; same arguments plus `--cpu-reference PATH` | independent explicit covariance jets checked against a locally generated `wick.py` result |
| `flow_calibration.py` | NumPy; no arguments | nonlinear single-variable SU(2) flow, exact moment factors and direct Haar quadrature |
| `lattice_flow_jets.py` | NumPy and JAX; no arguments | full-link local spatial-flow vector field on a 3x2x2x2 box |

The JAX programs explicitly enable double precision. Begin with the
smallest box: tensor storage grows cubically in the number of independent
edges, and automatic differentiation can allocate additional intermediates.
The CPU explicit-jet method is preferable for the small boxes tested.
The optional `--device gpu` selects a numerical backend, not a scheduler
allocation; resource acquisition is external to these programs.

For example, from the repository root:

```sh
project_root="$PWD"
scratch_dir=$(mktemp -d)
cd "$scratch_dir"
python3 "$project_root/strong-coupling/forest/wick/exact.py"
```

The adjacent-slice exact coefficient is
`R2 = 26843325164/2795053172625 > 0`. This is a finite-box unflowed
coefficient, not the infinite-volume flowed-vacuum correction.

To cross-check floating methods, generate a reference in one scratch
directory with `wick.py --device cpu --geometries 3x2`, then run
`jet.py --device cpu --geometries 3x2 --cpu-reference /path/to/reference/result.json`
in a different scratch directory. `--tree-order last` changes the maximal
tree; complete coefficients agree but individual sectors can change.

The full flowed Wick covariance, uniform remainder estimates, infrared
limits and B2 remain open. No program here establishes a continuum limit
or a mass-gap theorem. Historical operational records and raw outputs are
maintained separately from these scientific sources.
