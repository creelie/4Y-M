"""Lattice sums of the order-g^2 S' bound (LATTICE_SPRIME.md), lattice units a = 1, slice dimension 3.

Box side L = C sqrt(tau) (C = 4 as in leak.py), so the soft modes are the same at every spacing.
  H(Lambda) = Lambda * (1/V) sum_{k in BZ, |khat| >= Lambda} khat^-4     -> 1/(2 pi^2) = 0.0507 in the continuum
  Ebar      = (sum_P w_P 2 Phat) / (sum_P w_P) * sqrt(tau),  w_P = Phat^4 exp(-4 tau Phat^2) / (4 Phat^2)
              (Jensen bound of the leading-order one-step energy, times sqrt(tau))
  N(M)      = H(M) * Ebar / M: the hard-norm ratio ||w_Lambda||^2 / (c g^2 ||v||^2) at Lambda = M / sqrt(tau).
              The soft-leg factor [sum_P Phat^4 e^{-4 tau Phat^2}/(2 Phat)] / [sum_P w_P] is the same weighted
              mean of 2 Phat as in Ebar, so it equals Ebar / sqrt(tau).
"""
import numpy as np

C = 4.0
print(f"{'tau/a^2':>8} {'L':>5} {'Ebar=sqrt(tau)E_LO bound':>25} " + " ".join(f"H(M={M})".rjust(10) for M in (2, 4, 8)) + " " + " ".join(f"N(M={M})".rjust(10) for M in (2, 4, 8)))
for st in (1, 2, 4, 8, 16, 32, 64):
    tau = float(st) ** 2
    L = int(round(C * st))
    f = 2 * np.pi * np.fft.fftfreq(L)
    kx, ky, kz = np.meshgrid(f, f, f, indexing="ij")
    kh2 = 4 * (np.sin(kx / 2) ** 2 + np.sin(ky / 2) ** 2 + np.sin(kz / 2) ** 2)
    kh = np.sqrt(kh2)
    nz = kh2 > 0
    V = L ** 3
    w = np.where(nz, kh2 ** 2 * np.exp(-4 * tau * kh2) / (4 * np.where(nz, kh2, 1)), 0)
    Ebar = np.sum(w * 2 * kh) / np.sum(w) * st
    Hs = []
    for M in (2, 4, 8):
        Lam = M / st
        sel = kh >= Lam
        Hs.append(Lam * np.sum(kh[sel] ** -4.0) / V)
    print(f"{tau:8.0f} {L:5d} {Ebar:25.4f} " + " ".join(f"{h:10.4f}" for h in Hs) + " " + " ".join(f"{h*Ebar/M:10.4f}" for h, M in zip(Hs, (2, 4, 8))))
print(f"continuum H = 1/(2 pi^2) = {1/(2*np.pi**2):.4f} (infinite box)")
