# Numerical experiments behind forest/A_GAUSSIAN_FIBRE_4D.md

Each directory holds the code, its launcher (`run*.sh`) and, where available, the raw output it produced (`results*.jsonl`,
one JSON object per line). Launchers call `python` and `gcc`; they were written for a 16-core node, so reduce the
`xargs -P` / thread counts on smaller machines. Python needs numpy and scipy.

| Directory | Launcher → output | Backs |
|---|---|---|
| `gaussrho/` | `run.sh` → `results_rho.jsonl` | §8 ρ_j of the Gaussian fibre |
| | `run_angle.sh` → `results_angle_L4.jsonl`, `results_angle_L6.jsonl` | §9, §9.1 R_∠ (L = 4, 6) |
| | `run_cb.sh`, `run_cb4.sh` → `results_cb.jsonl`, `results_cb4.jsonl` | §8.3 first-chaos bound C_B^lin |
| | `run_massU.sh` → `results_massU_t.jsonl` (transverse); `results_massU_first.jsonl` is the first, longitudinal-dominated run | §10 massive-Gaussian (U) test |
| | `run_k0.sh` → `results_k0.jsonl` | §10 k → 0 limit, n ≤ 24 |
| | `run_kscan.sh` → `results_kscan.jsonl` | §10.2 k-infimum scan |
| | `run_twoscale.sh` → `results_twoscale.jsonl` | §10.3 free-field two-scale constant |
| `u1bp/` | `run.sh` → `results.jsonl` (via `analyze.py`) | §10.4 compact U(1)₃ block-Poincaré |
| `su2bp/` | `run.sh` → `results.jsonl` (via `analyze.py`) | §10.5 SU(2)₃ block-Poincaré |
| `pathcov/` | `run.sh`, `run2.sh`, `run3.sh` | Direction C path-covariance runs (Prop 168) |
| `linresp/`, `reduced_dobrushin/`, `vortex_mc/`, `polchinski2d/` | `run*.sh` | earlier exploratory checks |

Most scripts and `pathcov` have a `selftest` mode that the launcher runs first. The block-Poincaré codes are self-tested by
the n = L run in their launcher, which must give ϖ = 1.

Argument order for the block-Poincaré programs: `u1bp|su2bp L beta n nconf sep R s seed` (lattice size, coupling, cube side,
configurations, sweeps between configurations, inner conditional replicas per configuration, conditional sweeps between
replicas, RNG seed); output goes to stdout, and `analyze.py out/*.dat` turns it into ϖ with jackknife errors.
