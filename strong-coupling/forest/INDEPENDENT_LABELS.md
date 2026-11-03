# Independent-link labels: matching the good-window law to the graph

**Proved (deduction, with explicit conditioning).** Label the independent
physical **non-path** links, together with coarse curvature, rather than
all physical links. On a good window with a fixed finite collar, the full
conditional fibre law then equals the product graph law estimated in
`CANONICAL_CHART.md`. No restriction on the solved path links is needed.
This establishes a relative good-window estimate for this redesigned split;
it does not establish the bad-side or environment comparison, (I1), or (D).

## 1. A global branch bound from the constraint (proved)

Fix a pair k and its average c_k. Suppose all physical non-path links
entering its paths are within geodesic distance r of their reconstruction
values R_e(c). The eight outgoing private links include the path link ell;
the other seven and all inside-block links are non-path links. Every path
other than the two offset-zero copies of ell then satisfies

    d(P_x,c_k) <= 2r.

For **any** solution V^sharp_k=c_k, not assumed near R(c),

    d(ell,c_k) <= 14r.

*Proof.* If the cut-off chi is zero, V^sharp_k=ell, so ell=c_k. If chi>0,
all v_x=ell^(-1)P_x are in the principal-log region delta(v_x,I)<1.
For SU(2) the Hilbert--Schmidt norm of that principal logarithm equals
the bi-invariant geodesic distance from I. The two offset-zero terms are
zero. Thus, putting D=d(ell,c_k),

    D = d(exp(chi N_k),I)
      <= chi ||N_k||
      <= (14/16) max_x d(ell,P_x)
      <= (7/8)(D+2r).

Rearranging gives D<=14r. This argument covers every branch of the cut-off
map. In particular d(v_x,I)<=16r; for r<=10^(-4) the cut-off is identically
one at every solution. Every solution lies in the 44r contraction ball
of `CANONICAL_CHART.md` and hence coincides with its unique fixed point.
Existence was already supplied there. Thus, for fixed near-reconstruction
non-path links, the entire path-link fibre is a **single** point.

There is no claim of uniqueness when non-path links are large.

## 2. The minimal enlargement of a window (proved)

Let Z_global be the active links other than all path links d_k. For a
finite window W set

    K_W = outgoing pairs of W plus incoming pairs into W,
    L_W = A_W union {d_k : k is incoming into W},
    Z_W = L_W minus {d_k:k in K_W}.

Then Z_W consists of the 45 non-path active links per block of W. Unlike
`TOUCHING_CHART.md`, there is no need here to free the other seven private
links of an incoming pair: they stay in the exterior. The good-regime
path derivative is invertible, as already proved. The eight-link global
submersion lemma is not being applied to a bad incoming pair.

Freeze every physical link off L_W and all averages c. Retain the
constraints of every pair in K_W. Other averages do not involve L_W;
their compatibility with c is assumed, as it is for an exterior sampled
from the full conditional law. Label only the links of Z_global, never
the path links. Assume all labels in W are good and all frozen non-path
links in the supports of K_W are good. This imposes precisely

    d(U_e,R_e(c)) < r,     e in Z_W,

and fixed exterior restrictions. No label restricts any free path link.
Section 1 proves that the full fibre over this product is the single
smooth graph Gamma(D_r). There are no other path-link branches to count
or discard. Therefore Haar disintegration, multiplied by the Wilson
weight, gives **exactly** the graph density of `CANONICAL_CHART.md`.
This is an equality of the specified conditional laws, not a domain
sandwich or a comparison with a smaller event.

The choice of a window and its collar may be made from a fixed label
pattern before conditioning on the exterior. The law asserted here is
conditional on both that pattern and the exterior off L_W. Any additional
nonlocal information used to select a window would have to be retained;
it is not silently allowed by the statement.

## 3. Coarse labels and a sufficient collar (proved)

Let a coarse block be bad if either a non-path active link based there
is outside B(R_e(c),r), or a coarse plaquette based there has defect >eta.
Condition on a fixed pattern omega of these labels and on all c. A label
saying that a block is good restricts **each** of its non-path links;
it is a product condition. Coarse-curvature labels are constant on this
fibre. Bad labels elsewhere impose only exterior restrictions when W
and its collar are good.

A conservative sufficient condition is that all blocks in W^(+5) are
good, with periodic coarse l-infinity neighbourhoods. Pair bases K_W lie
in W^(+1). Plaquettes required by `smallframes` for those bases have
vertices in W^(+2), and their links have bases in W^(+4). For any path
link among those links, all its non-path pair data have bases in W^(+5).
Section 1 therefore puts even the **frozen exterior path links** within
14r of their reconstruction values. All non-path links there are within
r; all coarse plaquettes needed for R are bounded by eta.

Consequently every relevant plaquette along Gamma(D_r) has defect at most
eta+176r (this uses the looser 44r bound to agree with the previous
ledger). The proof of `CANONICAL_CHART.md` only needs small frozen
**non-path** coordinates for its fixed-point equations; frozen path links
enter the action but no other pair's constraint equation. Its estimates
therefore still hold with frozen path links within 44r. They also hold
with the smaller free set L_W: each independent coordinate still belongs
to at most two pair supports.

Under

    r>0,     eta>=0,     eta+352r <= 10^(-10),
    lambda>=2*10^11,

the actual conditional good-window law just described has

    Hess V_eff >= lambda/104,
    |Cov(F,H)| <= (208/lambda)
      ||exp(-alpha phi) grad F||_L2 ||exp(alpha phi) grad H||_L2,
    alpha=10^(-14),

in its independent coordinates, with phi 1-Lipschitz in the coarse base
metric. For physical observables with per-link gradients at most one and
supports X,Y, the corresponding bound is

    |Cov(F,H)| <= (3.3*10^6/lambda) sqrt(|X||Y|)
      exp[-10^(-14) max(dist_coarse(X,Y)-2,0)].

These constants are independent of volume, window shape, and values of
bad links beyond the good collar. Compatibility and goodness of the
specified collar are essential hypotheses. The bounds do not compare two
different exteriors; they are covariance bounds for each conditional law.

## 4. Rarity and remaining work

**Proved (deduction):** the non-path relative-link event is a subset of
the relative-link bad event in `RELATIVE_CHART.md`. The augmented labels
therefore satisfy the second-row witness and finite-intersection estimate
of `COARSE_LABELS.md`, at rho_bad=min(r/50,eta/105), for positive eta in its
stated range. No estimate uniform in c or the exterior is inferred.

**Proved (scope):** the explicit background and domain-matching checks are
now available for good windows of this independent-link split. The
counterexamples for rough backgrounds and the all-physical-link label
domain remain valid; the split was changed to address them.

**Sketch / open:**

- Bad non-path coordinates can have multiple path-link solutions and
  nonsmooth cut-off branches. The uniform bad-side response must include
  those branches and their disintegration density; it is not the old
  fixed-support single-link calculation.
- Good-window covariance must be turned into comparison across the island
  collar, with the transport and changes of the conditioning law included.
- The covariance of conditional expectations over the label environment
  still needs an estimate. Finite-intersection torus rarity alone does not
  supply the bootstrap recursion in `ANNEALED.md`.
- The five-block collar changes island sizes and entropy. No working
  island threshold or theorem (I1) is certified here.
- D8 and continuum existence remain open.

Thus one deterministic part of the single-scale proposal is repaired,
while the large-field mixture and multiscale analysis remain separate
research problems.
