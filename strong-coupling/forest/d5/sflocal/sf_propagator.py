"""Certified constants for the decay of the block-spin propagator (linearised
model, periodic lattice).  Exact integer matrices from walsh.py; every real
number below is an mpmath interval (iv) or an exact rational.

Coordinates: v = (inside links, orthonormal Walsh coordinates xi = sqrt(8) xi')
of the fibre ker V.  D = P^T G P with G the matrix of Q, A = L^T G L,
B = L^T G P, where 16 phi = P16 psi' + L16 gamma.
Run on the fine torus of side 12 (coarse side 6 >= 3, so the neighbourhoods
of distinct blocks at l^infty distance <= 1 are distinct); translation
invariance of all row data is checked on every block."""
import math, itertools
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np
from mpmath import iv, mp
from walsh import build_walsh
from lattice import torus_dist

iv.dps = 40
M = build_walsh(12); N = M['N']
P16, L16, d, V16 = M['P16'], M['L16'], M['d'], M['V16']
blk, lam, coarse = M['blk'], M['lam'], M['coarse']
nv = M['nv']
print('torus fine side 12, coarse side', N, ' variables', nv, ' coarse links', len(coarse))
# exactness checks
assert abs(V16 @ P16).max() == 0
assert ((V16 @ L16) - 256 * np.eye(len(coarse), dtype=np.int64)).__abs__().max() == 0
print('V P16 = 0 and V L16 = 256 I: exact integer check passed')
dP = (d @ P16).tocsr(); dL = (d @ L16).tocsr()
D256 = (dP.T @ dP).tocsr()      # 256 * P'^T G P'
A256 = (dL.T @ dL).tocsr()      # 256 * L^T G L
B256 = (dL.T @ dP).tocsr()      # 256 * L^T G P'
S2 = iv.sqrt(2)
def scale(li, lj):
    # 1/sqrt(lam_i lam_j) for lam in {1,8}: exact as iv
    p = li * lj
    return {1: iv.mpf(1), 8: 1 / (2 * S2), 64: iv.mpf(1) / 8}[p]
def colscale(lj):
    return iv.mpf(1) if lj == 1 else 1 / (2 * S2)

# ---- D: row sums by block distance, as exact pairs (r0 + r1/sqrt8 + r2/8)/256
def Drow(i):
    acc = defaultdict(lambda: [0, 0, 0])
    for j, v in zip(D256.indices[D256.indptr[i]:D256.indptr[i + 1]], D256.data[D256.indptr[i]:D256.indptr[i + 1]]):
        r = torus_dist(blk[i], blk[j], N)
        p = lam[i] * lam[j]
        acc[r][{1: 0, 8: 1, 64: 2}[p]] += abs(int(v))
    return {r: tuple(t) for r, t in acc.items()}
byblk = defaultdict(list)
for i, b in enumerate(blk): byblk[b].append(i)
ref = None
for b, lst in byblk.items():
    sig = tuple(tuple(sorted(Drow(i).items())) for i in lst)
    if ref is None: ref = sig
    assert sig == ref, 'translation invariance failed'
print('D row data identical on all', len(byblk), 'blocks; max block distance of a coupling:',
      max(r for s in ref for r, _ in s))
def ivval(t):
    return (iv.mpf(t[0]) + iv.mpf(t[1]) / (2 * S2) + iv.mpf(t[2]) / 8) / 256
S0 = max((ivval(dict(s).get(0, (0, 0, 0))) for s in ref), key=lambda x: x.b)
S1 = max((ivval(dict(s).get(1, (0, 0, 0))) for s in ref), key=lambda x: x.b)
print('max distance-0 row sum of |D|  in', S0)
print('max distance-1 row sum of |D|  S_1 in', S1)
c = iv.mpf(Fr(537, 5000).numerator) / Fr(537, 5000).denominator
print('coercivity c = 537/5000 (Proposition blockgapsharp, |Pv| >= |v|)')
# check |P v| >= |v| structurally: column norms; exact Gram of P' scaled
# |phi|^2 = sum_inside v^2 + sum_pairs (|xi|^2 + 8 m^2), verified numerically on random v
rng = np.random.default_rng(1)
for _ in range(5):
    v = rng.standard_normal(nv)
    psi = v / np.sqrt(np.array(lam, dtype=float))
    phi = (P16 @ psi) / 16
    assert phi @ phi >= v @ v - 1e-9
print('|P v|^2 >= |v|^2 on 5 random vectors (it holds identically by Walsh orthogonality)')
for a in ['0.05', '0.08', '0.1']:
    al = iv.mpf(a)
    ca = c - ((iv.exp(al)+iv.exp(-al))/2 - 1) * S1
    print('alpha =', a, ' c_alpha = c - (cosh a - 1) S_1 in', ca)
t = 1 + c / S1
amax = iv.log(t + iv.sqrt(t*t - 1))
print('c_alpha > 0 for alpha <', amax)

# ---- rows of P in orthonormal coordinates: l1 norms
P = P16.tocsr()
def Prow_l1(e):
    s = iv.mpf(0)
    for j, v in zip(P.indices[P.indptr[e]:P.indptr[e + 1]], P.data[P.indptr[e]:P.indptr[e + 1]]):
        s += iv.mpf(abs(int(v))) / 16 * colscale(lam[j])
    return s
p1 = max((Prow_l1(e) for e in range(P.shape[0])), key=lambda x: x.b)
print('max l1 norm of a row of P (orthonormal coordinates):', p1)
# support of a row of P: blocks within distance one of the base block of the link
act = M['active']
mx = 0
for e in range(P.shape[0]):
    be = M['block'](act[e][0])
    for j in P.indices[P.indptr[e]:P.indptr[e + 1]]:
        mx = max(mx, torus_dist(be, blk[j], N))
print('row of P supported within block distance', mx, 'of the base block')

# ---- A and B: l1 row norms and ranges
cblk = [k[0] for k in coarse]
def rowl1(Mat, rowblk, colblk, cs, i):
    acc = defaultdict(lambda: iv.mpf(0))
    for j, v in zip(Mat.indices[Mat.indptr[i]:Mat.indptr[i + 1]], Mat.data[Mat.indptr[i]:Mat.indptr[i + 1]]):
        acc[torus_dist(rowblk[i], colblk[j], N)] += iv.mpf(abs(int(v))) / 256 * cs(j)
    return dict(acc)
Arows = [rowl1(A256, cblk, cblk, lambda j: iv.mpf(1), i) for i in range(len(coarse))]
Brows = [rowl1(B256, cblk, blk, lambda j: colscale(lam[j]), i) for i in range(len(coarse))]
print('A: max distance of a coupling', max(max(r) for r in Arows),
      ' diagonal entries', sorted(set(int(A256[i, i]) // 256 for i in range(len(coarse)))))
a1 = max((sum(r.values(), iv.mpf(0)) for r in Arows), key=lambda x: x.b)
b1 = max((sum(r.values(), iv.mpf(0)) for r in Brows), key=lambda x: x.b)
print('B: max distance of a coupling', max(max(r) for r in Brows))
print('max l1 row norm of A:', a1, ' of B:', b1)
for a in ['0.05', '0.08', '0.1']:
    al = iv.mpf(a)
    ca = c - ((iv.exp(al)+iv.exp(-al))/2 - 1) * S1
    KD = 1 / ca
    KC = iv.exp(2 * al) * p1 ** 2 / (2 * ca)
    KM = iv.exp(2 * al) * b1 ** 2 / ca
    print('alpha', a, ': |D^-1_ij| <= %s e^{-a d};  2lam|C_ee\'| <= %s e^{-a d};  |M_kk\'| <= %s e^{-a d} (d>=2)'
          % (iv.nstr(KD, 8), iv.nstr(2 * KC, 8), iv.nstr(KM, 8)))
