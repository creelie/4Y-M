# Block shapes: is anything better than a plaquette block? (direction A3)

Before proving a Dobrushin bound for larger SU(2) blocks, this directory estimates how far each shape could
go at best. **Answer: not far.** Blocks built from the plaquette blocks of Lemma 11 reach a first-order
ceiling of about β_W ≈ 0.369, against 0.366 for plaquette blocks. A proof for 48-link blocks would gain under
1%, so it is not worth writing.

## The count

For an active link j, the leading-order influence count is

    n_j = Σ over the plaquettes p ∋ j of the number of blocks, other than the block of j, that meet p.

If the two other active links of p lie in one block, that plaquette counts once: the block sees them through
the single 1-Lipschitz product y_e y_e′. At small β the column sum is α_j ≈ (β_W/4)·n_j, so the
leading-order ceiling is 4/max_j n_j. For plaquette blocks n_j = 10 for every link (Lemma 11(iii)), giving
0.4. With the first-order correction K = (1+β)/4 this becomes 0.366.

A large block pays for its internal plaquettes. A link i in m_i plaquettes that lie inside its block passes
the response on with weight β along each of them. This is the z = y₂y₃ route that gives (1+β)/4 for
plaquette blocks. The first-order estimate used here is

    α_j(β_W) ≈ (β_W/4) Σ_{p ∋ j} Σ_{B ≠ B(j) met in p} (1 + β · max_{i ∈ p∩B} m_i),   β = β_W/4.

It ignores chains of length ≥ 2 inside a block, which only raise K, so these ceilings are optimistic.

## Results (periodic box of side 8, `structured.out`, `ceiling.out`)

The table groups block plaquettes by the coarse cell of their base point.

| cell | links per block | max n_j | mean n_j | 4/max n | first-order ceiling |
|---|---|---|---|---|---|
| 1·1·1·1 (plaquette blocks, paired) | 3, 6 | 10 | 9.50 | 0.400 | 0.366 |
| 2·1·1·1 | 3, 9 | 10 | 8.75 | 0.400 | 0.366 |
| 2·2·1·1 | 6, 18 | 10 | 7.75 | 0.400 | 0.356 |
| 2·1·2·1 | 12 | 10 | 8.06 | 0.400 | 0.353 |
| 2·2·2·1 | 24 | 10 | 6.94 | 0.400 | 0.348 |
| **2·2·2·2** | **48** | **9** | 5.63 | **0.444** | **0.369** |
| 4·2·2·2 | 96 | 9 | 5.00 | 0.444 | 0.369 |
| 4·4·2·2 | 192 | 9 | 4.36 | 0.444 | 0.347 |
| 4·4·4·4 | 768 | 9 | 3.02 | 0.444 | 0.347 |

- **The maximum is set by corner links.** Mean n_j falls quickly with block size. The maximum stays at 10 until
  the blocks are 4-dimensional cells, and then stays at 9 even for 768-link blocks: some link at a corner of
  each cell always sees nine other blocks.
- **Internal plaquettes eat the gain.** At the 2⁴ cell the leading ceiling rises to 0.444, but links with up to
  four internal plaquettes bring the first-order ceiling back to 0.369.
- **Annealing finds nothing better** (`merge.py`; groups of at most m plaquette blocks, P = 4, 3000 sweeps):

  | m | max n_j | mean n_j |
  |---|---|---|
  | 2 | 10 | 9.08 |
  | 3 | 10 | 8.73 |
  | 4 | 10 | 8.54 |
  | 8 | 10 | 8.23 |

  For m = 2 there is a proof: each plaquette block has 15 neighbouring blocks, and merging any adjacent pair
  lowers n_j for only 3 or 4 of its 6 links. `anneal.py` moves single links and never improves the maximum.

## Consequence

Plaquette blocks are essentially the best blocks of this family for the Dobrushin column sum. Getting SU(2)
past about 0.37 needs a different metric or method, not a larger block.

## Files

| File | Contents |
|---|---|
| `geometry.py` | Line field σ, the active links, the plaquettes and the plaquette blocks of Lemma 11 on a periodic box; checks n_j = 10 |
| `anneal.py` | Single-link annealing (negative result) |
| `merge.py` | Annealing over unions of plaquette blocks (`python3 merge.py m 4 3000 1`) |
| `structured.py` | Regular coarse-cell blocks: n_j, block sizes, internal plaquettes (`python3 structured.py 8`) |
| `ceiling.py` | First-order ceiling of each coarse-cell shape (`python3 ceiling.py 8`) |
