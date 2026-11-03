"""Split the order-g^2 change of Var(F) in leak.py data into its odd (cubic, F3) and even (F2-F4) parts.

Columns of data/leak/d4_L*.npy (big-bang-monotile/yang-mills-continuum, leak.py): F(+g), F(-g), G(+g), G(-g) at
g = 0.4, then F, G at g0 = 1e-3.  F(g) = g^2 F2 + g^3 F3 + g^4 F4 + ...
  odd  = (F(+g) - F(-g))/2 = g^3 F3 + O(g^5)
  even = (F(+g) + F(-g))/2 = g^2 F2 + g^4 F4 + O(g^6)
dVar/Var per g^2 = [Var(F(g))/g^4 - Var F2] / (g^2 Var F2) = [Var F3 + 2 Cov(F2,F4)]/Var F2 + O(g^2),
where Var over the +-g pooled samples equals Var(even) + Var(odd) (sign symmetry).
    python3 split_existing.py DATA_DIR
"""
import sys, glob
import numpy as np

g, g0 = 0.4, 1e-3
def parts(X):
    Fp, Fm, F0 = X[:, 0], X[:, 1], X[:, 4]
    F2 = F0 / g0**2
    v2 = F2.var()
    tot = (np.concatenate([Fp, Fm]).var() / g**4 - v2) / (g**2 * v2)
    odd = ((Fp - Fm) / 2).var() / g**6 / v2
    even = (((Fp + Fm) / 2).var() / g**4 - v2) / (g**2 * v2)
    return np.array([tot, odd, even])

def jack(X, nb=20):
    b = np.array_split(np.arange(len(X)), nb)
    full = parts(X)
    est = np.array([parts(np.delete(X, i, axis=0)) for i in b])
    return full, np.sqrt((nb - 1) / nb * np.sum((est - est.mean(0))**2, 0))

print(f"{'L':>3} {'tau/a^2':>8} {'n':>5}  {'total':>14} {'odd: Var F3/Var F2':>20} {'even: 2Cov(F2,F4)/Var F2':>26}")
for p in sorted(glob.glob(sys.argv[1] + '/d4_L*.npy'), key=lambda p: int(p.split('_L')[-1][:-4])):
    L = int(p.split('_L')[-1][:-4]); X = np.load(p)
    f, s = jack(X)
    print(f"{L:3d} {(L/4)**2:8.2f} {len(X):5d}  {f[0]:7.3f}({s[0]:.3f}) {f[1]:12.3f}({s[1]:.3f}) {f[2]:17.3f}({s[2]:.3f})")
