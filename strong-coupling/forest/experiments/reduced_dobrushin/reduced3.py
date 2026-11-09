"""Variant (b) stats: per-j column ratios, one j marginalized at a time. Variant (a): Z-term cross influences bounded by sup ||Hess_S log Z|| instead of beta*R.
MARGINAL_COUPLING_T1 step 5.1: Lipschitz (Dobrushin-type) row/column sums of the
reduced specification on F = V \ N vs the unreduced Wilson specification on V.
Weight exp(beta * Re U_p) per plaquette, beta = beta_W (Wilson: beta_W (1 - 1/2 Tr U)).
Raw potential-Lipschitz sums only (no vMF W1 conversion factor), per unit link distance."""
import itertools, json, sys
from scipy.special import iv, ive
import numpy as np
def hess_sup(b):
    # h(S)=log Z_b(|S|), Z=2I1(bs)/(bs). Hessian eigenvalues h''(s) (radial), h'(s)/s (tangential, x3).
    # h'(s)=b*I2(bs)/I1(bs). Sup of |eigs| over s in (0,6] on a fine grid.
    s=np.linspace(1e-6,6,200001); x=b*s
    hp=b*ive(2,x)/ive(1,x)
    hpp=np.gradient(hp,s)
    return float(max(np.abs(hpp).max(),(hp/s).max()))
d = 4
def unit(m): return tuple(1 if i == m else 0 for i in range(d))
def add(x, y, s=1): return tuple(a + s*b for a, b in zip(x, y))
def run(L, betaW):
    b = betaW
    R = hess_sup(b)  # second-derivative (cross-term) bound replaces beta*R(6beta)
    Rlip = iv(2, 6*b) / iv(1, 6*b)
    sites = list(itertools.product(range(L), repeat=d))
    # interior links: both endpoints in box
    V = {(x, m) for x in sites for m in range(d) if all(0 <= c < L for c in add(x, unit(m)))}
    def plaqs_of(link):
        x, m = link; out = []
        for n in range(d):
            if n == m: continue
            for s in (1, -1):
                y = x if s == 1 else add(x, unit(n), -1)
                out.append(frozenset({(y, m), (y, n), (add(y, unit(m)), n), (add(y, unit(n)), m)}) if s == 1 or True else None)
        return out
    P = set()
    for e in V: P.update(plaqs_of(e))
    bd = {l for p in P for l in p} - V
    def staples(e):  # links (other than e) in plaquettes of e
        return {l for p in plaqs_of(e) for l in p if l != e}
    # unreduced
    def wilson_row(f, Vset, excl=frozenset()):
        r_in = r_bd = 0.0
        for p in plaqs_of(f):
            if p & excl: continue
            for g in p - {f}:
                if g in Vset: r_in += b
                else: r_bd += b
        return r_in, r_bd
    un = [wilson_row(f, V) for f in V]
    un_row = max(a + c for a, c in un)
    un_col = {j: b * sum(1 for p in P if j in p for f in p if f in V) for j in bd}
    res = {"L": L, "beta_W": betaW, "R_lip": Rlip, "H_sup": R, "beta2_over_4": b*b/4, "|V|": len(V), "|bd|": len(bd),
           "unreduced_max_row": un_row, "unreduced_max_bd_col": max(un_col.values()),
           "unreduced_total_bd_infl": sum(un_col.values())}
    worst_row = 0.0; cols_red = []
    for j in bd:
        x, m = j
        N = {(add(x, unit(n), s), m) for n in range(d) if n != m for s in (1, -1)} & V
        F = V - N
        # j-column in reduced spec: Wilson plaquettes through j free of N, + beta R per e in N
        colj = 0.0
        for p in plaqs_of(j):
            if p & N: continue
            colj += b * sum(1 for f in p if f in F)
        for e in N: colj += R * len(staples(e) & F)  # j moves S_e; each f in staples of e sees it
        cols_red.append((colj, un_col[j]))
        for f in F:
            ri, rb = wilson_row(f, F, excl=frozenset(N))
            for e in N:
                st = staples(e)
                if f in st:
                    ri += R * len((st & F) - {f}); rb += R * len(st - V)
            worst_row = max(worst_row, ri + rb)
    res.update({"reduced_max_row": worst_row,
                "reduced_max_bd_col": max(c for c, _ in cols_red),
                "col_ratio_red_over_un_max": max(c / u for c, u in cols_red),
                "col_ratio_quantiles": [float(q) for q in np.quantile([c / u for c, u in cols_red], [0, .25, .5, .75, 1])],
                "total_bd_infl_red": sum(c for c, _ in cols_red),
                "n_j_improved": sum(1 for c, u in cols_red if c < u), "n_j": len(cols_red)})
    return res
for L in (2, 3):
    for bw in (0.30, 0.351):
        print(json.dumps(run(L, bw)), flush=True)
