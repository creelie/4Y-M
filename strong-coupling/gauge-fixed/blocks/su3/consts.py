"""Constants of the SU(3) plaquette-block bound, in rigorous floating point with error bounds.

For b in r conv SU(3) let nu_b(dU) ~ exp(Re Tr(Ub)) dU, phi_a(U) = Re Tr(U E_a) for the orthonormal basis E_a of
C^{3x3} (E_ij, i E_ij), m(b) the mean of phi and H(b) its covariance (the Hessian of log Z).  On every box Q of
side r/2^L of the cover (diagonal representatives b = diag(a1, a2, x3 + i y3), as in ../../checks/su3/deg2) this
program bounds
  tau  >= every eigenvalue of H(b),
  lam  <= every eigenvalue of H(b),
  M    >= sup_{|u|=|v|=1} |D^3 log Z(b)(u, v, v)|,  the third cumulant of phi,
for all b in Q, from the exact moments of degree <= 3 at the sixteen corners and at the centre b0:
  (i)   Var_b(phi_v) <= nu_b[(phi_v - k_v)^2] <= e^{s/8} max_i nu_{b_i}[(phi_v - k_v)^2]   (k fixed, Lemma 21 (ii));
  (ii)  Var_b(phi_v) >= e^{-s/8} min_i Var_{b_i}(phi_v)   (the left inequality of (21) and Lemma 21 (i));
  (iii) with k = m(b0), delta = m(b) - k, |delta| <= tau rho + e_k, rho = max |b - b0|_HS,
        kappa3_b(u,v,v) = T_b(u,v,v) - delta_u Var(phi_v) - 2 delta_v Cov(phi_u, phi_v) - delta_u delta_v^2,
        T_b = nu_b[(phi - k)^{(x)3}], and |T_b(u,v,v)| <= (max_i ||T_{b_i}||_unf + eps max_i nu_{b_i}|g|)/(1 - eps),
        eps = 1 - e^{-s/8}, ||T||_unf the spectral norm of the 18 x 324 unfolding, nu|g| <= 2 sqrt(3) lambda_max(S2).
Here s = w1^2 + w2^2 + max(w3^2, w4^2) as in Lemma 21.  Usage:
  python3 consts.py --r 77/50 --L 4 --cache r77_50/cover_cache.pkl --out r77_50
(four worker processes; corner moments are read from and added to the cache).  Writes consts_boxes.txt and
consts_log.txt in the output directory."""
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')   # one BLAS thread per worker: four workers already fill four cores
import os, sys, time, pickle, argparse, itertools
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'checks', 'su3', 'deg2'))
os.chdir(os.path.join(HERE, '..', '..', 'checks', 'su3', 'deg2'))
from multiprocessing import Pool
import flint
from box2 import corner2, full, lam, gk
from polys import PHI2, N3, R as RING
import cover2
arb = flint.arb

AL = (PHI2 / 2)
LIN = [i for i in range(N3) if RING.degs[i] == 1]
A1 = AL[:, LIN]                      # phi_a = sum_p A1[a, p] x_p over the 18 degree-one monomials x_p
MON1 = [RING.mons[i] for i in LIN]
def _ix(*ps):
    return RING.idx[tuple(sum(c) for c in zip(*[MON1[p] for p in ps]))]
I2 = np.array([[_ix(p, q) for q in range(18)] for p in range(18)])
I3 = np.array([[[_ix(p, q, s) for s in range(18)] for q in range(18)] for p in range(18)])
I1 = np.array(LIN)
AA = np.abs(A1)
EPSR = gk(18 ** 3)                   # relative rounding allowance for every contraction below

def raw(cd):
    """first, second, third raw moments of phi at one corner, with absolute error bounds"""
    w, e = full(cd)
    m1 = (A1 @ w[I1]).real
    e1 = AA @ e[I1] + EPSR * (AA @ np.abs(w[I1]))
    m2 = np.einsum('ap,bq,pq->ab', A1, A1, w[I2]).real
    e2 = np.einsum('ap,bq,pq->ab', AA, AA, e[I2] + EPSR * np.abs(w[I2]))
    m3 = np.einsum('ap,bq,cs,pqs->abc', A1, A1, A1, w[I3]).real
    e3 = np.einsum('ap,bq,cs,pqs->abc', AA, AA, AA, e[I3] + EPSR * np.abs(w[I3]))
    return m1, e1, m2, e2, m3, e3

def centred(m1, e1, m2, e2, m3, e3, k):
    """S2 = nu[(phi-k)(phi-k)^T], T = nu[(phi-k)^{(x)3}] for an exact vector k, with error bounds"""
    ak = np.abs(k)
    S2 = m2 - np.outer(m1, k) - np.outer(k, m1) + np.outer(k, k)
    E2 = e2 + np.outer(e1, ak) + np.outer(ak, e1) + 8 * gk(4) * (np.abs(m2) + 2 * np.outer(np.abs(m1), ak) + np.outer(ak, ak))
    sym = lambda X: X + X.transpose(1, 2, 0) + X.transpose(2, 0, 1)
    T = m3 - sym(np.einsum('ab,c->abc', m2, k)) + sym(np.einsum('a,b,c->abc', m1, k, k)) - np.einsum('a,b,c->abc', k, k, k)
    E3 = (e3 + sym(np.einsum('ab,c->abc', e2, ak)) + sym(np.einsum('a,b,c->abc', e1, ak, ak))
          + 8 * gk(8) * (np.abs(m3) + sym(np.einsum('ab,c->abc', np.abs(m2), ak)) + sym(np.einsum('a,b,c->abc', np.abs(m1), ak, ak))
                         + np.einsum('a,b,c->abc', ak, ak, ak)))
    return S2, E2, T, E3

def lam_min(V, E):
    """certified lower bound of the smallest eigenvalue of the symmetric matrix within E of V"""
    return -lam(-V, E)

def unf(T, E):
    """certified upper bound of the spectral norm of the 18 x 324 unfolding of the tensor within E of T"""
    U = T.reshape(18, 324); EU = E.reshape(18, 324)
    G = U @ U.T
    EG = np.abs(U) @ EU.T + EU @ np.abs(U).T + EU @ EU.T + gk(330) * (np.abs(U) @ np.abs(U).T)
    return float(np.sqrt(lam(G, EG))) * (1 + 1e-15)

def box_consts(args):
    lo, hi, cds, cc = args
    corners = cover2.corners_of(lo, hi)
    wd = hi - lo
    sig = arb(float(wd[0])) ** 2 + arb(float(wd[1])) ** 2 + max(arb(float(wd[2])) ** 2, arb(float(wd[3])) ** 2)
    up = float((sig / 8).exp().upper()); dn = float((-sig / 8).exp().lower()); eps = 1 - dn
    rho = float((arb(float(wd[0])) ** 2 + arb(float(wd[1])) ** 2 + arb(float(wd[2])) ** 2 + arb(float(wd[3])) ** 2).sqrt().upper()) / 2
    m1c, e1c = raw(cc)[:2]
    k = m1c.copy(); ek = float(np.sqrt((e1c ** 2).sum())) * (1 + 1e-14)
    lS = 0.0; lmin = 9.0; un = 0.0
    for c in corners:
        m1, e1, m2, e2, m3, e3 = raw(cds[c])
        S2, E2, T, E3 = centred(m1, e1, m2, e2, m3, e3, k)
        lS = max(lS, lam((S2 + S2.T) / 2, E2))
        Cv = m2 - np.outer(m1, m1)
        Ecv = e2 + np.outer(e1, np.abs(m1)) + np.outer(np.abs(m1), e1) + np.outer(e1, e1) + 4 * gk(4) * (np.abs(m2) + np.outer(np.abs(m1), np.abs(m1)))
        lmin = min(lmin, lam_min((Cv + Cv.T) / 2, Ecv))
        un = max(un, unf(T, E3))
    tau = arb(up) * arb(lS)
    lamb = arb(dn) * arb(lmin)
    D = tau * arb(rho) + arb(ek)
    Tb = (arb(un) + arb(eps) * 2 * arb(3).sqrt() * arb(lS)) / (1 - arb(eps))
    M = Tb + 3 * tau * D + D ** 3
    return lo, hi, float(tau.upper()), float(lamb.lower()), float(M.upper())

def boxes(r, L, L0=3):
    """the boxes of side r/2^L meeting R_r, generated by halving the boxes of side r/2^L0 as in cover2.py, so
    that the corner coordinates coincide with those of the cover and its cache"""
    w = r / 2 ** L0; level = []
    for i1 in range(2 ** L0):
        for i2 in range(2 ** L0):
            for i3 in range(2 ** (L0 + 1)):
                for i4 in range(2 ** L0):
                    lo = np.array([i1 * w, i2 * w, -r + i3 * w, i4 * w]); level.append((lo, lo + w))
    for l in range(L0, L + 1):
        level = [(lo, hi) for lo, hi in level if not (cover2.misses_order(lo, hi) or cover2.separated(lo, hi))]
        if l == L:
            return level
        nxt = []
        for lo, hi in level:
            h = (hi - lo) / 2
            for c in itertools.product((0, 1), repeat=4):
                clo = lo + np.array(c) * h; nxt.append((clo, clo + h))
        level = nxt

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--r', default='77/50'); ap.add_argument('--L', type=int, default=4)
    ap.add_argument('--cache', default=None); ap.add_argument('--out', default=None)
    a = ap.parse_args()
    r = float(Fraction(a.r)); cover2.R = r
    out = os.path.abspath(os.path.join(HERE, a.out or 'r' + a.r.replace('/', '_'))); os.makedirs(out, exist_ok=True)
    cpath = os.path.abspath(os.path.join(HERE, a.cache)) if a.cache else os.path.join(out, 'consts_cache.pkl')
    logf = os.path.join(out, 'consts_log.txt')
    def log(*x):
        s = ' '.join(str(y) for y in x); print(s, flush=True)
        with open(logf, 'a') as f: f.write(s + '\n')
    cache = {}
    if os.path.exists(cpath):
        with open(cpath, 'rb') as f: cache = pickle.load(f)
    t0 = time.time()
    B = boxes(r, a.L)
    centres = [tuple((lo + hi) / 2) for lo, hi in B]
    need = sorted({c for lo, hi in B for c in cover2.corners_of(lo, hi)} | set(centres))
    need = [c for c in need if c not in cache]
    log(f'r = {a.r}, side r/2^{a.L}: {len(B)} boxes meet R_r, {len(need)} new moment points')
    with Pool(4) as pool:
        done = 0
        for c, d in zip(need, pool.imap(corner2, need, chunksize=2)):
            cache[c] = d; done += 1
            if done % 500 == 0:
                log(f'  points {done}/{len(need)} {time.time() - t0:.0f}s')
                with open(cpath, 'wb') as f: pickle.dump(cache, f)
        with open(cpath, 'wb') as f: pickle.dump(cache, f)
        tasks = [(lo, hi, {c: cache[c] for c in cover2.corners_of(lo, hi)}, cache[ce]) for (lo, hi), ce in zip(B, centres)]
        res = list(pool.imap(box_consts, tasks, chunksize=8))
    with open(os.path.join(out, 'consts_boxes.txt'), 'w') as f:
        f.write('# lower corner, upper corner in (a1, a2, x3, y3); tau (upper), lambda (lower), M (upper)\n')
        for lo, hi, t, l, m in res:
            f.write(' '.join(f'{x:.9f}' for x in lo) + '   ' + ' '.join(f'{x:.9f}' for x in hi) + f'   {t:.6f} {l:.6f} {m:.6f}\n')
    tau = max(x[2] for x in res); lamb = min(x[3] for x in res); M = max(x[4] for x in res)
    log(f'{len(res)} boxes, {time.time() - t0:.0f}s: tau <= {tau:.6f}, lambda >= {lamb:.6f}, M <= {M:.6f}')

if __name__ == '__main__':
    main()
