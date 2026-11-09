"""Linear-response influence probe for 4D SU(2) Wilson theory on a periodic torus.

Feasibility study #3 of strong-coupling/CLAY_NEXT_DIRECTIONS.md.
For a perturbation of the action at plaquette p (d/dt of beta*Tr U_p/2), the
response of a gauge-invariant observable O is beta*Cov(O, P_p) exactly.
We measure the slice-summed plaquette correlator C(r) = Cov(S_0, S_r),
S_t = sum of plaquettes on time slice t, and the influence ratio
    kappa(r) = beta * |C(r)| / Var(S_0) * (number of boundary plaquettes ratio 1),
i.e. the response of one slice to unit perturbation of all plaquettes of a slice
at distance r, normalised. Contraction at physical scale needs kappa(r) < 1 with
r small; we report the smallest such r (L*) and the effective mass.
Heuristic: measures the weak (periodic-state) form only.
"""
import argparse, json, math, sys, time
import numpy as np

def qmul(a, b):
    a0,a1,a2,a3 = np.moveaxis(a,-1,0); b0,b1,b2,b3 = np.moveaxis(b,-1,0)
    return np.stack([a0*b0-a1*b1-a2*b2-a3*b3, a0*b1+a1*b0+a2*b3-a3*b2,
                     a0*b2-a1*b3+a2*b0+a3*b1, a0*b3+a1*b2-a2*b1+a3*b0], -1)
def qc(a): return a*np.array([1.,-1.,-1.,-1.])
def sh(x, mu, s):  # x(n + s e_mu); lattice axes are 0..3 of x
    return np.roll(x, -s, axis=mu)

def staple(U, mu):
    S = np.zeros(U.shape[1:])
    for nu in range(4):
        if nu == mu: continue
        a = qmul(qmul(sh(U[nu],mu,1), qc(sh(U[mu],nu,1))), qc(U[nu]))
        b = qmul(qmul(qc(sh(sh(U[nu],mu,1),nu,-1)), qc(sh(U[mu],nu,-1))), sh(U[nu],nu,-1))
        S += a + b
    return S

def kp_heatbath(k, rng):
    """sample a0 with density sqrt(1-a0^2) exp(k a0) (Kennedy-Pendleton)."""
    out = np.empty(k.shape); todo = np.ones(k.shape, bool)
    while todo.any():
        kk = k[todo]
        r1, r2, r3, r4 = (1-rng.random((4, kk.size)))
        lam2 = -(np.log(r1) + np.cos(2*np.pi*r2)**2*np.log(r3))/(2*kk)
        ok = r4**2 <= 1 - lam2
        idx = np.flatnonzero(todo)[ok]
        out[idx] = 1 - 2*lam2[ok]; todo[idx] = False
    return out

def sweep(U, beta, par, rng):
    for mu in range(4):
        for p in (0, 1):
            A = staple(U, mu)                     # sum of staples, quaternion
            nA = np.linalg.norm(A, axis=-1)
            k = beta*nA                            # action beta/2 Tr(U A) = beta*<U,A>
            a0 = kp_heatbath(k[par==p], rng)
            v = rng.normal(size=(a0.size, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
            X = np.concatenate([a0[:, None], np.sqrt(1-a0**2)[:, None]*v], 1)
            Abar = qc(A[par==p]/nA[par==p][:, None])
            U[mu][par==p] = qmul(X, Abar)

def plaq(U):
    P = np.zeros(U.shape[1:5])
    for mu in range(4):
        for nu in range(mu+1, 4):
            P += qmul(qmul(U[mu], sh(U[nu],mu,1)), qmul(qc(sh(U[mu],nu,1)), qc(U[nu])))[..., 0]
    return P / 6

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, required=True)
    ap.add_argument("--L", type=int, default=8); ap.add_argument("--T", type=int, default=16)
    ap.add_argument("--therm", type=int, default=300); ap.add_argument("--meas", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--out", default=None)
    a = ap.parse_args(); rng = np.random.default_rng(a.seed)
    shape = (a.T, a.L, a.L, a.L)
    U = np.zeros((4,)+shape+(4,)); U[..., 0] = 1
    par = np.indices(shape).sum(0) % 2
    t0 = time.time()
    for _ in range(a.therm): sweep(U, a.beta, par, rng)
    S = np.empty((a.meas, a.T)); pl = np.empty(a.meas)
    for i in range(a.meas):
        sweep(U, a.beta, par, rng)
        P = plaq(U); S[i] = P.sum(axis=(1, 2, 3)); pl[i] = P.mean()
    S -= S.mean(0)
    C = np.array([np.mean(S*np.roll(S, -r, 1)) for r in range(a.T//2+1)])
    # jackknife errors over 20 bins
    nb = 20; bs = a.meas//nb; Cj = []
    for b in range(nb):
        m = np.ones(a.meas, bool); m[b*bs:(b+1)*bs] = False
        Sb = S[m]-S[m].mean(0)
        Cj.append([np.mean(Sb*np.roll(Sb, -r, 1)) for r in range(a.T//2+1)])
    Cj = np.array(Cj); dC = np.sqrt((nb-1)*Cj.var(0))
    kappa = 6*a.beta*np.abs(C)/C[0]*1.0   # response per unit slice perturbation (6 plaquettes/site)
    meff = [float(np.log(C[r]/C[r+1])) if C[r] > 0 and C[r+1] > 0 else None for r in range(a.T//2)]
    res = dict(beta=a.beta, L=a.L, T=a.T, meas=a.meas, plaquette=float(pl.mean()),
               C=C.tolist(), dC=dC.tolist(), meff=meff, kappa=kappa.tolist(),
               Lstar=next((r for r in range(1, a.T//2+1) if kappa[r] + 6*a.beta*dC[r]/C[0] < 1), None),
               seconds=time.time()-t0)
    print(json.dumps(res), flush=True)
    if a.out: open(a.out, "w").write(json.dumps(res))

if __name__ == "__main__": main()
