"""Weak-coupling (Gaussian) approximation of the forest gauge: exact cube coefficients.

Linearise U = exp(i A) with A = 0 on the fixed links F; the action becomes
beta sum_p (dA)_p^2, and beta drops out of the response. For a Gaussian the
conditional law of the cube V is a translate whose mean depends linearly on the
boundary links, so the coefficient C_{V,j} of the paper is the l1 norm of the
response column, exactly and for every exterior:
    K_c = sum_{j in c, j not in V} sum_{e in V} |(Q_VV^{-1} Q_Vj)_e|.
Prints max, mean and min of K_c / N_n over the 48 classes. Needs scipy.
Run: python3 gaussian_forest.py 2 4 6 8
Nothing in the paper rests on it."""
import itertools, sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from influence_mc import SIGMA
def sig(x): return SIGMA[tuple(v % 2 for v in x)]
def run(n, off=2):
    V = {}
    for x in itertools.product(range(off, off + n), repeat=4):
        for mu in range(4):
            if sig(x) != mu: V[(x, mu)] = len(V)
    E = lambda x, m: tuple(v + (1 if k == m else 0) for k, v in enumerate(x))
    rows = []  # plaquettes meeting V: list of (link, sign)
    seen = set()
    for (x, mu) in V:
        for nu in range(4):
            if nu == mu: continue
            for base in (x, tuple(v - (1 if k == nu else 0) for k, v in enumerate(x))):
                a, b = min(mu, nu), max(mu, nu)
                key = (base, a, b)
                if key in seen: continue
                seen.add(key)
                links = [((base, a), 1), ((E(base, a), b), 1), ((E(base, b), a), -1), ((base, b), -1)]
                rows.append([(l, s) for l, s in links if sig(l[0]) != l[1]])  # F links are fixed to 0
    B = {}
    for r in rows:
        for l, s in r:
            if l not in V and l not in B: B[l] = len(B)
    Qvv = sp.lil_matrix((len(V), len(V))); Qvb = sp.lil_matrix((len(V), len(B)))
    for r in rows:
        for l1, s1 in r:
            if l1 not in V: continue
            for l2, s2 in r:
                if l2 in V: Qvv[V[l1], V[l2]] += s1 * s2
                else: Qvb[V[l1], B[l2]] += s1 * s2
    lu = spl.splu(sp.csc_matrix(Qvv))
    Qvb = sp.csc_matrix(Qvb); C = np.zeros(len(B))
    for s0 in range(0, len(B), 800):   # response of the cube means to each boundary link
        M = lu.solve(Qvb[:, s0:s0 + 800].toarray())
        C[s0:s0 + 800] = np.abs(M).sum(0)
    K = {}
    for (x, mu), i in B.items():
        c = (tuple(v % 2 for v in x), mu)
        K[c] = K.get(c, 0) + C[i]
    Nn = (n // 2) ** 4
    k = np.array(list(K.values())) / Nn
    print(f"n={n} |V|={len(V)} boundary={len(B)} classes={len(K)}  K_c/N_n: max {k.max():.4f} mean {k.mean():.4f} min {k.min():.4f}", flush=True)
for n in map(int, sys.argv[1:]): run(n)
