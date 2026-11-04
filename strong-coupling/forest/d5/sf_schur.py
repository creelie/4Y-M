"""Exact rational computation for the obstruction (smallfield.tex).

Scalar linearised model on the W-fibre of Sec. ssec:largefields: real phi on
the active links with base in the blocks of W, phi = 0 on tree links and on all
links off A_W; q(gamma) = min sum_{p in P_W} (dphi)_p^2 subject to
V_{y,mu}(phi) = gamma_{y,mu} for y in W (Definition def:blockgauge).
q(gamma) = gamma^T M gamma.  We compute M exactly (python-flint fmpq_mat)
for W a row of n blocks along e_0 and report M_{AA}, M_{AB} for the coarse
links A = (0, 1), B = ((n-1) e_0, 1).
"""
import itertools, sys
from fractions import Fraction
import flint
D = 4
def off(x): return tuple(c % 2 for c in x)
def blk(x): return tuple(c // 2 for c in x)
def add(a, b): return tuple(i + j for i, j in zip(a, b))
def e(m): return tuple(1 if i == m else 0 for i in range(D))
O16 = [tuple(o) for o in itertools.product((0, 1), repeat=D)]
def is_tree(x, m):
    o = off(x); return o[m] == 0 and all(o[n] == 0 for n in range(m + 1, D))

def model(n):
    W = [tuple([k, 0, 0, 0]) for k in range(n)]
    Wset = set(W)
    var = {}
    for y in W:
        for o in O16:
            x = add(tuple(2 * c for c in y), o)
            for m in range(D):
                if not is_tree(x, m):
                    var[(x, m)] = len(var)
    assert len(var) == 49 * n
    # plaquettes containing a variable link
    plaq = set()
    for (x, m) in var:
        for nn in range(D):
            if nn == m: continue
            a, b = min(m, nn), max(m, nn)
            # plaquettes in plane {a,b} containing link (x,m): corner x or x - e_other
            other = nn
            for corner in (x, tuple(c - d for c, d in zip(x, e(other)))):
                plaq.add((corner, a, b))
    rows = []
    for (x, a, b) in plaq:
        L = [((x, a), 1), ((add(x, e(a)), b), 1), ((add(x, e(b)), a), -1), ((x, b), -1)]
        r = {}
        for l, s in L:
            if l in var: r[var[l]] = r.get(var[l], 0) + s
        rows.append(r)
    nv = len(var)
    # G = D^T D
    G = [[0] * nv for _ in range(nv)]
    for r in rows:
        for i, si in r.items():
            for j, sj in r.items():
                G[i][j] += si * sj
    # constraints, multiplied by 16
    cons = []
    for y in W:
        for m in range(D):
            c = {}
            for o in O16:
                x = add(tuple(2 * cc for cc in y), o)
                for z in (x, add(x, e(m))):
                    if (z, m) in var: c[var[(z, m)]] = c.get(var[(z, m)], 0) + 1
            cons.append(((y, m), c))
    return var, rows, G, cons

def schur(n):
    var, rows, G, cons = model(n)
    nv, nc = len(var), len(cons)
    N = nv + nc
    K = flint.fmpq_mat(N, N)
    for i in range(nv):
        for j in range(nv):
            if G[i][j]: K[i, j] = 2 * G[i][j]
    for k, (_, c) in enumerate(cons):
        for i, s in c.items():
            K[nv + k, i] = s; K[i, nv + k] = s
    # right-hand sides: 16*gamma = unit vectors (constraint rows are 16 V)
    Rhs = flint.fmpq_mat(N, nc)
    for k in range(nc):
        Rhs[nv + k, k] = 16
    Sol = K.solve(Rhs)
    assert K * Sol == Rhs          # exact rational solution of the KKT system
    # q(e_k) etc: phi^T G phi; M_{kl} = phi_k^T G phi_l
    phis = [[Sol[i, k] for i in range(nv)] for k in range(nc)]
    def form(a, b):
        s = flint.fmpq(0)
        for r in rows:
            da = sum((s_ * a[i] for i, s_ in r.items()), flint.fmpq(0))
            db = sum((s_ * b[i] for i, s_ in r.items()), flint.fmpq(0))
            s += da * db
        return s
    labels = [lab for lab, _ in cons]
    kA = labels.index(((0, 0, 0, 0), 1)); kB = labels.index(((n - 1, 0, 0, 0), 1))
    MAA = form(phis[kA], phis[kA]); MAB = form(phis[kA], phis[kB])
    # check the constraint and the stationarity are exact
    return len(var), len(rows), MAA, MAB

def dec(q, digits=8):
    # decimal rounding of an exact rational, sign and 'digits' significant digits
    from fractions import Fraction
    f = Fraction(int(q.p), int(q.q))
    return '%.*e' % (digits - 1, f)

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    print("Row of n blocks along e_0; A=(0,1), B=((n-1)e_0,1); exact rationals, rounded:")
    for n in range(1, nmax + 1):
        nv, npl, MAA, MAB = schur(n)
        print("n=%d: |A_W|=%d, n_P=%d, M_AA=%s, M_AB=%s" % (n, nv, npl, dec(MAA), dec(MAB)),
              flush=True)
        if n == 4:
            print("  n=4 exact M_AB =", MAB)
            assert MAB < 0

