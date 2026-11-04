# Flat holonomy obstructs an absolute bad-link curvature witness

**Proved (counterexample to the unconditional formulation of input (b)).**
On a periodic four-dimensional fine lattice with even side 2L, divide sites
into 2^4 blocks and use the block tree gauge of `main.tex`. Fix a direction
mu. Set every fine link in that direction that crosses from one block to the
next equal to -I in SU(2); set every other link equal to I. All block-tree
links are I, since they remain inside their blocks, so this is in the stated
tree gauge. Every plaquette is I: a plaquette involving mu contains two
crossing mu-links or none, and the other links are I. Nevertheless each
crossing mu-link is -I, with Hilbert--Schmidt distance 2 sqrt(2) from I.

When L is odd, the Wilson line around the periodic mu direction is -I.
It is invariant under periodic gauge transformations up to conjugation,
so no periodic gauge transformation can make all links I. This topological
observation is stronger than needed for the displayed tree-gauge example.
The coarse path link across such a block is -I as well; flatness does not
mean that a fixed coarse background or its holonomy is trivial.

Consequently a statement asserting that *every absolute tree-gauge link*
outside a fixed small ball about I forces a nearby large plaquette is false
without additional hypotheses. In particular, it cannot be inferred solely
from local plaquette flatness and the block tree gauge. This counterexample
does not rule out a witness for a gauge-invariant or background-relative
bad variable, or for a specified subset of non-tree links within a block.

## A corrected local target

**Proved (relative path comparison).** For a compact matrix group with
bi-invariant Hilbert--Schmidt distance, if two oriented paths with the same
endpoints have a loop that can be reduced by N elementary plaquette moves,
then the distance between their parallel transports is at most N times the
largest distance of an involved plaquette from I. Indeed, each move inserts
a conjugate of a plaquette holonomy or its inverse; bi-invariance and the
triangle inequality give the bound by telescoping. No small-field
assumption or commutativity is used.

**Sketch / open.** The plausible repair for input (b) is to define badness
relative to a locally chosen coarse parallel transport or to compare paths
with the same endpoints, then prove that this definition is compatible with
the constrained chart, labels and chessboard estimate. The required
reference need not be close to I and may carry global holonomy. In
particular, a relative witness does not automatically give the rarity or
conditioning estimates needed for (I1). The existing path-block estimates
in `main.tex` should be checked for the exact constants and domain before
using this observation in a gauge proof.
