# Relative link labels do not control the coarse background

**Proved (explicit counterexample, with exact finite verification).** The event
that every physical link lies within any positive radius of the canonical
reconstruction R(c) does not imply a positive constrained action Hessian,
even locally at its centre. Thus a relative-chart theorem for this event
needs a separate coarse-background hypothesis. This does not contradict
`chartdecay`, which already assumes small background curvature.

## 1. A centre-valued background

Take the coarse torus (Z/2Z)^4, so the fine torus has side 4. Directions are
numbered 0,1,2,3. Prescribe the SU(2) coarse links

    c_(y,mu) = (-I)^(sum over nu<mu of y_nu).

Every coarse plaquette is -I, including those across the periodic seam.
Let U=R(c): links inside blocks are I and links leaving B(y) in direction
mu are c_(y,mu). This is in tree gauge, and V^sharp(U)=c. All sixteen
paths of any given pair equal c_(y,mu), so the cut-off is inactive, even
though the coarse curvature is large. Every relative link label is good:

    d(U_e,R_e(c)) = 0.

## 2. A four-link curve on the exact fibre

Fix a Hilbert--Schmidt unit generator T in su(2), and let xi_e be zero
except at the following links (base point, direction, scalar value):

| fine base | direction | xi_e |
|---|---:|---:|
| (0,1,1,1) | 1 | -1 |
| (0,1,1,1) | 2 | +1 |
| (0,1,0,0) | 1 | +1 |
| (0,0,1,0) | 2 | -1 |

All four links are active and private; the last two are path links.
Define U_e(t)=U_e exp(t xi_e T). The two nonzero private links of each
changed pair have opposite values, each appearing twice in its sixteen
paths. Thus every linear average of xi is exactly zero.

For sufficiently small |t|, all v_x are near I and the cut-off is one.
All perturbations commute with the central background, so the principal
logarithms in the average are exactly their scalar linear expressions:

    P_x(t) = c_k exp(t a_x T),
    v_x(t) = exp(t(a_x-a_0)T),
    V^sharp_k(U(t)) = c_k exp(t [mean_x a_x]T) = c_k.

This proves that U(t) stays in the full nonlinear fibre, not just its
tangent space. Also d(U_e(t),R_e(c)) <= |t|, so for every r>0 a neighbourhood
of t=0 stays within the all-good physical link-ball event of radius r.

## 3. Negative constrained curvature

Write U_p=s_p I, s_p in {-1,+1}, and let q_p=(d xi)_p be the signed scalar
curl. For a Hilbert--Schmidt unit SU(2) generator, eigenvalues of T are
+/- i/sqrt(2). Therefore

    Q(U(t)) = sum_p 4[1-s_p cos(t q_p/sqrt(2))],
    Q''(0) = 2 sum_p s_p q_p^2.

The only nonzero contributions are:

| s_p | q_p^2 | number of plaquettes |
|---:|---:|---:|
| +1 | 1 | 12 |
| -1 | 1 | 10 |
| -1 | 4 | 1 |

Consequently

    sum s_p q_p^2 = 12-10-4 = -2,
    Q''(0) = -4,     T(xi)=4,
    (beta_W Q/4)''(0) = -beta_W.

The table and all 64 zero-average identities are checked with integers by
`d5/rough_certificate.py`, which uses only the Python standard library.
The table is an exhaustive count on the explicitly specified finite torus;
there is no floating-point sign decision in this certificate.
Every ambient first derivative of Q vanishes at U, since every plaquette
is central. Thus this is also a negative intrinsic Hessian direction of
Q on the fibre at U, independently of the chosen smooth coordinates.
The constrained Haar density is smooth and positive near U: the cut-off
is identically one and the path-link derivative is invertible there. Its
logarithmic derivatives in any fixed smooth chart are finite and do not
depend on beta_W. Hence its contribution cannot restore a positive Hessian
proportional to beta_W for all sufficiently large beta_W in that chart.
This last statement concerns local convexity; it is not a lower bound on
mixing time or a proof that a particular covariance estimate fails.

## 4. Consequence and a repair to investigate

**Proved (scope correction).** Canonical relative-link goodness alone
cannot supply input (a). A rough c can put the centre of the all-good
event at an unstable stationary point. The bad-side fixed-support
calculation and the witness r/50 do not address this missing hypothesis.

**Sketch (repair).** Include coarse-plaquette badness in the environment,
or restrict the relative-chart statement to coarse backgrounds whose
reconstruction is calm. The latter permits a quantitative, inhomogeneous
private-link construction around R(c), with the frozen exterior treated
as parameters. It still requires matching the actual label domain and
controlling the environment term. A coarse-curvature label itself has a
fine-plaquette witness, so it can be included in annealed torus rarity;
this adds losses and does not solve the multiscale problem D8.

**Computation (discovery only).** `d5/rough_hessian.py` forms the full
signed quadratic form on the same torus and finds a negative eigenvalue
numerically. It led to the four-link certificate above. The final result is the exact four-link certificate above.
