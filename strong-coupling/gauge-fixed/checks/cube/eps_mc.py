"""Monte Carlo estimate of eps_V for SU(2) in the forest gauge.

eps_V(beta) is the mean chordal distance from the identity of the links of the
cube V of side n under the gauge-fixed measure with the identity outside V.
The theorem on what the cube estimate requires says that D(n, alpha) forces
eps_V >= (1 - alpha) D / 3 = 2 (1 - alpha) / 3. This script samples the
measure by heat bath and prints eps_V; a value well above 2/3 means that this
necessary condition excludes nothing at that coupling and size.

Run: python3 eps_mc.py --beta 2.3 --n 6 [--chains 4 --therm 300 --sweeps 300]
Prints one line per run; nothing in the paper rests on it.
"""
import argparse
import numpy as np
from influence_mc import Box, SIGMA


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, required=True, help="beta_W")
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--therm", type=int, default=300)
    ap.add_argument("--sweeps", type=int, default=300)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    n, P = a.n, a.n + 4
    idx = np.indices((P,) * 4)
    incube = np.ones((P,) * 4, bool)
    for d in range(4):
        incube &= (idx[d] >= 2) & (idx[d] < 2 + n)
    V = np.stack([(SIGMA[tuple(idx % 2)] != mu) & incube for mu in range(4)])
    box = Box(P, a.chains, rng)
    one = np.array([1.0, 0.0, 0.0, 0.0])
    series = []
    for t in range(a.therm + a.sweeps):
        box.heatbath(a.beta, rng, mask_extra=V)
        if t >= a.therm:
            d = np.linalg.norm(box.U - one, axis=-1)            # (B, 4, P, P, P, P)
            series.append((d * V[None]).sum(axis=(1, 2, 3, 4, 5)) / V.sum())
    s = np.array(series)                                        # (sweeps, chains)
    per_chain = s.mean(axis=0)
    half = len(s) // 2
    print(f"beta_W={a.beta} n={n} |V|={int(V.sum())} eps_V={per_chain.mean():.4f} "
          f"(chains sd {per_chain.std():.4f}, late-minus-early "
          f"{s[half:].mean() - s[:half].mean():+.4f}) threshold 2/3", flush=True)


if __name__ == "__main__":
    main()
