<p align="center"><img src="assets/banner.svg" alt="The Yang–Mills Mass Gap from One Bound at One Scale" width="100%"></p>

<p align="center">
  <a href="https://doi.org/10.5281/zenodo.23444338"><img src="https://img.shields.io/badge/DOI-10.5281%2Fzenodo.23444338-22d3ee?style=for-the-badge" alt="DOI"></a>
  <a href="papers/4Y-MPaper.pdf"><img src="https://img.shields.io/badge/paper-PDF-f472b6?style=for-the-badge" alt="paper"></a>
  <img src="https://img.shields.io/badge/Lean_4-265_theorems-7c5cff?style=for-the-badge" alt="Lean 4">
  <img src="https://img.shields.io/badge/Julia-234_checks-34d399?style=for-the-badge" alt="Julia">
  <img src="https://img.shields.io/badge/license-MIT-f59e0b?style=for-the-badge" alt="MIT">
</p>

**Deep Bhattacharjee** (corresponding author, head of the project), Ushashi Bhattacharya, James Montgomery Lockwood, Tristan Romanov, Ikailo John Sessford, Sanjeevan Singha Roy.

| Proved | Open |
| :--- | :--- |
| Mass gap of 4D lattice Yang–Mills theory at strong coupling, every compact simple group, explicit bounds | The one-scale bound at weak coupling |
| **Theorem F**: a Wightman limit on tori and one bound on partition functions at one scale give the continuum mass gap and a unique vacuum | The Wightman limit itself |
| The unrenormalised loop variables have a trivial limit | The Clay problem |

| Folder | Contents | Run |
| :--- | :--- | :--- |
| [`papers/`](papers) | the paper | |
| [`YangMills/`](YangMills) | Lean 4 project over Mathlib | `lake exe cache get && lake build` |
| [`scripts/`](scripts) | axiom audit | `lake env lean scripts/Axioms.lean` |
| [`julia/`](julia) | constants of the paper | `julia julia/constants.jl` |
