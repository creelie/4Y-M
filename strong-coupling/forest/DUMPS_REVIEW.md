# Review of dump1–dump8 (group work, 4–5 Nov 2026)

I read all of `dump1`, `dump3`, `dump5`, `dump7` and `dump8`, and skimmed `dump2`, `dump4` and `dump6` for structure and key
statements. The dumps total about 1 MB, roughly 250k tokens. `dump2` and `dump4` are near-duplicate copies of one
combined Lean file.

## What's there

- **`dump1` (Deep).** The paper's three open targets:
  - (H1): polynomial covariance decay for block variables, uniform over every torus.
  - (H2): a split inequality at weak coupling at every scale, with constants uniform in k and L.
  - (N_u): a nontrivial limit, plus the Osterwalder–Schrader inputs.
  - Deep says each is Clay-level and none is proved. I agree.
- **`dump2`/`dump4`.** One Lean file of about 18 parts that compiles with no `sorry`. Its content falls into three groups:
  - The Dobrushin comparison and coupling machinery (James's main contribution).
  - Some small library facts: Wilson action bounds, the one-loop running coupling, Holley–Stroock, a cumulant
    counterexample, block and tree counting, and a strong-coupling threshold β_W < 2/9.
  - SpectralNoGo and LimitRatio, which are exact facts about spectral measures and limits.
- **`dump5`–`dump7`.** Translations of the Lean statements into mathematical notation, and Studio status: all 213
  declarations translated, 168 checked as consistent, 45 still unchecked.
- **`dump3`/`dump8`.** A link to the MDPI preprint, and a note confirming that the PDF's Lean compiles. The note also
  flags two real issues: an eigenvalue floor at 10⁻⁴ (the note calls it "the Romanov eigenvalue floor"), and `hW`, the
  per-site Dobrushin coupling hypothesis, which nothing in the file proves for the SU(2) kernel.

## Does any of it help toward the Clay problem?

Very little, though the work is careful and correct.

- **The Lean is honest, but it's the easy layer.** What it checks is abstract: a Dobrushin bound gives decay and
  comparison; one-loop running tends to zero; the strong-coupling regime is β_W < 2/9. Each Yang–Mills-specific input
  enters as a hypothesis such as `hW`. The `hW_example` instance is a two-site Boolean toy. Nothing proves `hW` for
  SU(2) at any coupling where the continuum limit lives.
- **The no-go and ratio theorems are useful for keeping us honest.** They show that mass-gap ratios are independent of
  renormalisation, and that a fixed lattice gap gives no continuum gap. That pins down what a proof must deliver, but
  doesn't deliver it.
- **None of it touches the walls our own notes reached.** Those walls are the large-field expansion (D4/D7, R1), the
  crossover base case at β ≈ 10–700 (D8), and E. Deep's (H1)/(H2) are the same walls in different notation: (H1) on all
  tori is our D5/D7 with decay, and (H2) is D8.

## Is the problem formulated sensibly?

Mostly, with three caveats:
1. **The reduction is coherent.** Proving (H1), (H2) and (N_u) would give a gap via the paper's chain, and the spectral
   no-go shows why the gap must be uniform in the lattice spacing. That part is correct.
2. **(H1) asks for too much.** It takes a supremum over all boundary data V. Our `CALM_LINK.md` §4 indicates that
   uniform-field exteriors break pointwise statements of this kind. (H1) should probably be restated with an average
   over Gibbs-distributed exteriors, as in `AVERAGED_D4.md`, or it may simply be false as written.
3. **(H2) contains the crossover problem.** It asks for constants uniform in k from weak coupling down to the coarsest
   scale. That includes β of order 1–700, where neither strong- nor weak-coupling tools work. Calling it "weak
   coupling" understates it.

## How close are we?

Far. The proved parts are:
- strong coupling (β_W ≲ 0.35);
- abstract comparison and coupling theorems;
- single-scale estimates at absurdly large β (10³⁸ to 10⁹⁵);
- local averaged statements, such as Conjecture B (`GAUSS_GREEN.md`).

Still missing is everything that makes Yang–Mills hard: multiscale large-field control, the crossover, and
nontriviality. Bałaban spent about a decade on the first and never finished it. I'd rate us at the level of a clean,
partly formal framework, with no progress on the core difficulty.

## Should our next efforts be based on it?

Not as a research direction, but in two practical ways:
- **Use the Lean framework as the place to land results.** Our Conjecture B or averaged-D4 lemmas could be stated in its
  Dobrushin and coupling language and formalised there.
- **Fix the formulation before attacking it.** I'd restate (H1) in averaged form and check whether the sup-over-V version
  fails under uniform-field boundaries. That's a concrete, achievable question, and the answer changes which target is
  worth attacking.

Two other things are worth doing:
- **The `hW` gap the note flags.** Proving `hW` for the SU(2) heat-bath kernel in the strong-coupling range
  (β_W < 2/9) is feasible, and would make the Lean chain unconditional there.
- **The 10⁻⁴ eigenvalue floor.** It's worth finding where it comes from before anyone relies on it.
