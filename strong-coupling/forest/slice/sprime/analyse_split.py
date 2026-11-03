"""Analyse leak_split.py output (LATTICE_SPRIME.md, item 3).  Per unit g^2, relative to leading order:
  odd   = Var F3 / Var F2                     (Fock norm^2 of the cubic part, relative to ||F2 Omega||^2)
  hard  = E Var(F3 | theta_s) / Var F2        (= ||w_Lambda||^2/||v||^2 per g^2, Lambda = 2/sqrt(tau))
  soft  = odd - hard
  even  = 2 Cov(F2, F4) / Var F2              (normalisation change of the two-particle sector)
  Eo    = sqrt(tau) * mean energy of F3|Omega> (leading order)
  de    = order-g^2 coefficient of the even-sector energy, relative to E0 = Ebar/sqrt(tau)  (eta - eps)
  E2/E0 = de + odd (Eo/E0 - 1)                (total, comparable with leak_norm.out)
    python3 analyse_split.py data/split_L*.npy
"""
import sys
import numpy as np

g, g0 = 0.4, 1e-3


def Ebar(L):
    tau = (L / 4) ** 2
    f = 2 * np.pi * np.fft.fftfreq(L)
    k = np.meshgrid(f, f, f, indexing="ij")
    kh2 = sum(4 * np.sin(x / 2) ** 2 for x in k)
    kh = np.sqrt(kh2)
    w = np.where(kh2 > 0, kh2 * np.exp(-4 * tau * kh2) / 4, 0)
    return np.sum(w * 2 * kh) / np.sum(w) * np.sqrt(tau)


def NM(L, M=2):
    tau = (L / 4) ** 2
    f = 2 * np.pi * np.fft.fftfreq(L)
    k = np.meshgrid(f, f, f, indexing="ij")
    kh = np.sqrt(sum(4 * np.sin(x / 2) ** 2 for x in k))
    lam = M / np.sqrt(tau)
    H = lam * np.sum(kh[kh >= lam] ** -4.0) / L ** 3
    return H * Ebar(L) / M


def stats(X, L):
    st = np.sqrt((L / 4) ** 2)
    E0 = Ebar(L) / st
    A, B = X[:, :5], X[:, 5:]
    P = np.concatenate([A, B])
    F2 = P[:, 4] / g0 ** 2
    v2 = F2.var()
    o = (P[:, 0] - P[:, 1]) / 2
    e = (P[:, 0] + P[:, 1]) / 2
    odd = o.var() / g ** 6 / v2
    F3a, F3b = (A[:, 0] - A[:, 1]) / (2 * g ** 3), (B[:, 0] - B[:, 1]) / (2 * g ** 3)
    hard = np.mean((F3a - F3b) ** 2) / 2 / v2
    even = ((e / g ** 2).var() - v2) / (g ** 2 * v2)
    Eo = g ** 2 / 2 * P[:, 3].mean() / o.var()
    Ee = g ** 2 / 2 * P[:, 2].mean() / e.var()
    de = (Ee / E0 - 1) / g ** 2
    tot = de + odd * (Eo / E0 - 1)
    return np.array([odd, hard, odd - hard, even, Eo * st, de, tot])


def jack(X, L, nb=20):
    b = np.array_split(np.arange(len(X)), nb)
    full = stats(X, L)
    est = np.array([stats(np.delete(X, i, axis=0), L) for i in b])
    return full, np.sqrt((nb - 1) / nb * np.sum((est - est.mean(0)) ** 2, 0))


if __name__ == "__main__":
    names = ["odd", "hard", "soft", "even", "Eo*st", "dEe", "E2/E0"]
    print(f"{'L':>3} {'tau':>6} {'n':>4} {'N(2)':>6} {'Ebar':>6} " + " ".join(f"{n:>13}" for n in names))
    for p in sorted(sys.argv[1:], key=lambda p: int(p.split("_L")[-1][:-4])):
        L = int(p.split("_L")[-1][:-4])
        X = np.load(p)
        f, s = jack(X, L)
        print(f"{L:3d} {(L/4)**2:6.2f} {len(X):4d} {NM(L):6.3f} {Ebar(L):6.3f} " +
              " ".join(f"{a:7.3f}({b:.3f})" for a, b in zip(f, s)))
