"""Split of the order-g^2 norm and energy of the flowed slice state (LATTICE_SPRIME.md, open item 3).

Same set-up as big-bang-monotile/yang-mills-continuum/leak.py (free transverse vacuum, exact nonlinear lattice
spatial Wilson flow, F = zero-momentum spatial plaquette energy, box L = 4 sqrt(tau)), with three additions:
  * the parity split: F = even + odd under theta -> -theta; grad F(+g) = de + do, grad F(-g) = -de + do, so the
    even and odd parts of the energy numerator sum (dF/dphi)^2 are measured separately;
  * the soft/hard split at Lambda = M / sqrt(tau): theta = theta_s + theta_h with independent soft and hard parts
    (Fourier modes with phat < Lambda, >= Lambda).  For each soft draw, two independent hard draws are taken;
    E[Var(F3 | theta_s)] = E[(F3^1 - F3^2)^2]/2 is the squared Fock norm of the components of F3|Omega> with at
    least one hard leg (conditional expectation on theta_s = projection onto all-soft Fock states);
  * only F and its gradient are needed.
Each output row, for hard draw i = 1, 2:  F(+g), F(-g), |de|^2, |do|^2, F(g0)   (g = 0.4, g0 = 1e-3).
    python3 leak_split.py L NSAMP SEED M OUT.npy
"""
import sys, time
import numpy as np
sys.path.insert(0, "/home/user/big-bang-monotile/yang-mills-continuum")
import su2_slice_jax as SJ

G, G0, EPS, d = 0.4, 1e-3, 0.1, 4


def sampler(L, rng):
    n = 3
    axes = (1, 2, 3)
    ps = np.meshgrid(*([2 * np.pi * np.arange(L) / L] * n), indexing="ij")
    k = np.stack([1 - np.exp(-1j * p) for p in ps])
    k2 = np.sum(np.abs(k) ** 2, axis=0)
    zero = k2 < 1e-12
    k2s = np.where(zero, 1.0, k2)
    omega = np.sqrt(k2s)

    def fourier():
        out = np.empty((n,) + (L,) * n + (3,), complex)
        for c in range(3):
            eta = np.fft.fftn(rng.normal(size=(n,) + (L,) * n), axes=axes, norm="ortho")
            proj = eta - np.conj(k) * (np.sum(k * eta, axis=0) / k2s)
            th = proj / np.sqrt(2 * omega)
            th[:, zero] = 0.0
            out[..., c] = th
        return out

    def real(th):
        return np.fft.ifftn(th, axes=axes, norm="ortho").real
    return fourier, real, np.sqrt(k2)


def build(L, tau):
    K = SJ.build(d, True)
    jax, jnp = K["jax"], K["jnp"]
    n = 3
    shape = (1,) + (L,) * n
    nsteps = max(1, int(np.ceil(tau / EPS)))
    eps = tau / nsteps
    step = jax.checkpoint(lambda V: K["rk3_step"](V, eps))

    def F_of(phi):
        U = jnp.zeros((d,) + shape + (4,)).at[..., 0].set(1.0)
        U = U.at[1:].set(K["qexp_vec"](0.5 * phi[:, None]))
        V = jax.lax.scan(lambda V, _: (step(V), None), U, None, length=nsteps)[0]
        return K["slice_energy_raw"](V)[0]
    # gradient w.r.t. the left generator: perturb U -> exp(e/2) U; at first order this is the same as
    # leak.py's convention for the norm; we differentiate through phi with the left-generator chart.
    def F_left(e, phi):
        U = jnp.zeros((d,) + shape + (4,)).at[..., 0].set(1.0)
        U = U.at[1:].set(K["qmul"](K["qexp_vec"](0.5 * e[:, None]), K["qexp_vec"](0.5 * phi[:, None])))
        V = jax.lax.scan(lambda V, _: (step(V), None), U, None, length=nsteps)[0]
        return K["slice_energy_raw"](V)[0]
    FG = jax.jit(lambda phi: jax.value_and_grad(F_left)(jnp.zeros_like(phi), phi))
    Fonly = jax.jit(F_of)
    return FG, Fonly, nsteps


def main():
    L, nsamp, seed, M, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    tau = (L / 4) ** 2
    rng = np.random.default_rng(seed)
    fourier, real, phat = sampler(L, rng)
    hard = (phat >= M / np.sqrt(tau))[None, ..., None]
    FG, Fonly, nsteps = build(L, tau)
    print(f"L={L} tau/a^2={tau:.2f} M={M} hard modes={hard.mean():.3f} flow steps={nsteps}", flush=True)
    rows, t0 = [], time.time()
    for i in range(nsamp):
        A = fourier()
        row = []
        for _ in range(2):
            B = fourier()
            th = real(np.where(hard, B, A))
            fp, gp = FG(G * th)
            fm, gm = FG(-G * th)
            gp, gm = np.asarray(gp), np.asarray(gm)
            de, do = (gp - gm) / 2, (gp + gm) / 2
            row += [float(fp), float(fm), float(np.sum(de * de)), float(np.sum(do * do)), float(Fonly(G0 * th))]
        rows.append(row)
        if i % 10 == 0:
            print(f"  sample {i} t={time.time() - t0:.0f}s", flush=True)
            np.save(out, np.array(rows))
    np.save(out, np.array(rows))


if __name__ == "__main__":
    main()
