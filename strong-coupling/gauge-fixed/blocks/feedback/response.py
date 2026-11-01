"""Linear response of the SU(2) plaquette-block law to one outside field (direction A2).

mu_b(dy) ~ exp(sum_i <y_i, b_i> + 4 beta Re(y1 y2 y3)) on (S^3)^3.  For a unit h, the derivative of mu_{b1 + t h}(f)
is Cov_mu(f, <y1, h>).  K_true(b) = sup over f 1-Lipschitz for d = sum |y_i - y_i'| of this covariance.
We bound K_true from below with test families whose best member is linear in a free vector:
  (A) <y1,h1> + <y2,k2> + <y3,k3>          (sum of the three marginal responses, linear part)
  (B) <y1,h1> + <y2 y3, k>                  (the z = y2 y3 route; 1-Lipschitz since |d(y2 y3)| <= |dy2| + |dy3|)
  (C) <y1,h1> + <y1 y2, k>... not allowed (y1 twice);  (C) is <y3,k3> + <y1 y2,k>  (response through y3 and y1y2)
  (D) <y1 y2 y3, k> + 0                     (whole plaquette)
For a family f = sum_a <u_a(y), k_a> the best unit k_a is Cov(u_a, <y1,h>)/|.|, so the bound is sum_a |Cov(u_a, <y1,h>)|.
Sampling: exact heat-bath sweeps over the three links, many chains in parallel, vMF on S^3 by rejection.
Usage: python3 response.py BETA_W NCHAINS SWEEPS
"""
import sys
import numpy as np

rng = np.random.default_rng(1)


def qmul(a, b):
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([a0*b0 - a1*b1 - a2*b2 - a3*b3, a0*b1 + a1*b0 + a2*b3 - a3*b2,
                     a0*b2 - a1*b3 + a2*b0 + a3*b1, a0*b3 + a1*b2 - a2*b1 + a3*b0], -1)


def conj(a):
    return a * np.array([1, -1, -1, -1])


def vmf(c):
    """one sample of nu_c(dx) ~ exp(<c,x>) dx on S^3 for each row of c."""
    n = len(c); k = np.linalg.norm(c, axis=1); out = np.empty((n, 4)); todo = np.arange(n)
    while len(todo):
        x = rng.normal(size=(len(todo), 4)); x /= np.linalg.norm(x, axis=1, keepdims=True)
        acc = rng.random(len(todo)) < np.exp(np.einsum('ij,ij->i', x, c[todo]) - k[todo])
        out[todo[acc]] = x[acc]; todo = todo[~acc]
    return out


def run(bW, b1, b2, b3, h, n, sweeps, burn=30):
    beta = bW / 4
    y = [vmf(np.zeros((n, 4)) + 0.0) for _ in range(3)]
    acc = {k: [] for k in 'y1 y2 y3 z23 z12 P'.split()}
    s1 = []
    for sw in range(burn + sweeps):
        # Re(y1 y2 y3) = <y1, conj(y2 y3)> = <y2, conj(y3 y1)> = <y3, conj(y1 y2)>
        y[0] = vmf(b1 + 4 * beta * conj(qmul(y[1], y[2])))
        y[1] = vmf(b2 + 4 * beta * conj(qmul(y[2], y[0])))
        y[2] = vmf(b3 + 4 * beta * conj(qmul(y[0], y[1])))
        if sw >= burn:
            s1.append(y[0] @ h)
            acc['y1'].append(y[0]); acc['y2'].append(y[1]); acc['y3'].append(y[2])
            acc['z23'].append(qmul(y[1], y[2])); acc['z12'].append(qmul(y[0], y[1]))
            acc['P'].append(qmul(qmul(y[0], y[1]), y[2]))
    s1 = np.concatenate(s1); s1 = s1 - s1.mean()
    cov = {}
    for k, v in acc.items():
        v = np.concatenate(v); cov[k] = (v - v.mean(0)).T @ s1 / len(s1)
    nrm = {k: np.linalg.norm(v) for k, v in cov.items()}
    se = 1 / np.sqrt(len(s1))  # rough: each covariance has sd ~ 0.25/sqrt(N_eff)
    return nrm, se


if __name__ == '__main__':
    bW = float(sys.argv[1]); n = int(sys.argv[2]); sweeps = int(sys.argv[3])
    beta = bW / 4; R = 20 * beta
    e = np.eye(4)
    cases = {
        'all zero': (0 * e[0], 0 * e[0], 0 * e[0]),
        'b2 max': (0 * e[0], R * e[0], 0 * e[0]),
        'b2,b3 max': (0 * e[0], R * e[0], R * e[1]),
        'b2 max, b3 max aligned': (0 * e[0], R * e[0], R * e[0]),
    }
    for name, (b1, b2, b3) in cases.items():
        best = None
        for h in (e[0], e[1], (e[0] + e[1]) / np.sqrt(2)):
            nr, se = run(bW, b1, b2, b3, h, n, sweeps)
            fam = {'A y1+y2+y3': nr['y1'] + nr['y2'] + nr['y3'], 'B y1+z23': nr['y1'] + nr['z23'],
                   'C y3+z12': nr['y3'] + nr['z12'], 'D P': nr['P'], 'y1 alone': nr['y1']}
            line = ' '.join(f'{k}={v:.4f}' for k, v in fam.items())
            print(f'{name:24s} h={np.round(h, 2)}  {line}', flush=True)
    print('reference: paper K(bW) bound; (1+beta)/4 =', (1 + beta) / 4)
