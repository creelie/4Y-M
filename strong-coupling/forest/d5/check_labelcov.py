"""Certified arithmetic for the subsection ssec:labelcov (interval arithmetic, mpmath.iv).
Inputs from the paper (Theorem thm:islandwindow, Lemma lem:lambdatwo, Proposition
prop:islandmass): alpha_G = 2.5e-4, rho = 2.17e-16, r = 43 rho, beta* (enclosures),
log10(M q*) (upper bounds), bar Lambda_2 = 99.0001 beta + 4.65e6,
sigma(kappa) = sum_{y in Z^4} exp(-kappa |y|_inf) in closed form.
Every claimed inequality is checked on interval enclosures; a check passes only if it
holds for the whole interval.  All quantities are evaluated at beta = beta*, which
bounds them for beta_W >= beta* by the monotonicity proved in Corollary
cor:fibrecovnumbers (b_* <= (beta_W/beta*) b_*(beta*); b_* eps^{1/2}, b_* eps and g_R
decrease)."""
from mpmath import iv
iv.dps = 50
ln10 = iv.log(10)
def log10(x): return iv.log(x)/ln10
def hi(x): return x.b
def lo(x): return x.a
ok = []
def check(name, cond):
    ok.append((name, bool(cond))); print(("PASS " if cond else "FAIL ") + name)
def show(name, x): print(f"{name} in [{float(lo(x)):.6g}, {float(hi(x)):.6g}]")
def sig(kap):
    x = iv.exp(-kap)
    return 1 + 64*x*(1+4*x+x**2)/(1-x)**4 + 16*x/(1-x)**2

aG = iv.mpf('2.5e-4')
rho = iv.mpf('2.17e-16'); r = 43*rho
s2 = sig(2*aG); s1 = sig(aG)
show("sigma(2 alpha_G)", s2)
check("sigma(2 alpha_G) <= 6.145e15 (value in prop:islandmass)", hi(s2) <= 6.145e15)
show("sigma(alpha_G)", s1)
check("sigma(alpha_G) <= 9.831e16", hi(s1) <= 9.831e16)
check("11^4 = 14641", 11**4 == 14641)

data = {1: dict(beta=iv.mpf(['1.12205e37','1.12207e37']), L10e=iv.mpf('-14102.4')),
        5: dict(beta=iv.mpf(['1.32162e37','1.32164e37']), L10e=iv.mpf('-16652.7'))}
out = {}
for s, D in data.items():
    print(f"\n=== s = {s} ===")
    beta = D['beta']; L10e = D['L10e']; lneps = L10e*ln10
    Lam = iv.mpf('99.0001')*beta + iv.mpf('4.65e6')
    bstar = 1920*Lam**2*r**2/beta
    nu = 45*(2*s+9)**4
    bn = bstar*nu
    show("b_*", bstar); print("nu_s =", nu); show("log10(b_* nu_s)", log10(bn))
    check(f"s={s}: beta* rho^2/16 > 8 (b_* eps^(1/2), b_* eps decrease)", lo(beta*rho**2/16) > 8)
    check(f"s={s}: 2 sqrt(eps) <= 1", hi(L10e) < -1)
    R = 2*s+9
    lam = -(iv.log(iv.sqrt(2)) + 4*iv.log(2*R+1) + lneps/4)   # -ln g_R at R = 2s+9
    show("lambda = -ln g_R at R = 2s+9", lam)
    check(f"s={s}: g_(2s+9) <= 1/2", lo(lam) >= hi(iv.log(2)))
    # Phi_R at R = 2s+9: b nu sigma(alpha) eps^{1/2} e^{alpha(2s+8)}
    lPhi = log10(bn*s1) + L10e/2 + aG*(2*s+8)/ln10
    show("log10 Phi_(2s+9)", lPhi)

    # ---- Corollary (i): labels, R = 2s+9, t = R(ceil(d/(1e4 R)) + 2), d >= 6s+34
    check(f"s={s}: (6s+34)(1-1e-4) >= 6s+33 (d >= t+6)", lo((6*s+34)*(1-iv.mpf('1e-4'))) >= 6*s+33)
    cL = log10(14641*bn) + 2*aG*(8*s+41)/ln10
    show("log10 c_s = log10(11^4 b nu e^{2a(8s+41)})", cL)
    fL = log10(iv.mpf(2)) + lPhi
    show("log10 2 Phi_R (relative floor)", fL)
    lamL = lam/(10**4*R)
    show("chain rate lambda/(1e4 R)", lamL)
    # d at which the main relative term reaches 1, with V = 1
    d0 = cL*ln10/(2*aG*(1-iv.mpf('1e-4')))
    show("d where 10^{c_s} e^{-2a(1-1e-4)d} = 1", d0)
    out[s] = dict(bn=hi(log10(bn)), lamR=lo(lam/R), cL=hi(cL), fL=hi(fL), lamL=lo(lamL), d0=hi(d0))

    # ---- Corollary (ii): fibre, regime E, R = 2s+9, t = R ceil(d/(1e4 R)), d >= 2s+10
    C1 = iv.mpf(9)/2*bn*iv.exp(4*aG*(2*s+9))
    show("log10 C_1 = log10((9/2) b nu e^{4a(2s+9)})", log10(C1))
    check(f"s={s}: lambda/(1e4 R) >= 2 alpha_G", lo(lamL) >= hi(2*aG))
    bnlo = 1920*(iv.mpf('99.0001')*lo(beta)+iv.mpf('4.65e6'))**2*r**2/lo(beta)*nu
    check(f"s={s}: 10 <= (9/4) b nu", lo(iv.mpf(9)/4*bnlo) >= 10)
    check(f"s={s}: (2s+10)(1-1e-4) > 2s+9 (t < d)", lo((2*s+10)*(1-iv.mpf('1e-4'))) > 2*s+9)
    aF = log10(iv.mpf(17)/4*bn*s1) + L10e + aG*(2*s+8)/ln10
    show("log10 (17/4) eps^{1/2} Phi_R (fibre floor)", aF)
    out[s].update(C1=hi(log10(C1)), aF=hi(aF))

    # ---- Corollary (iii): fibre, regime S, R = ceil(2000 sqrt d), t = floor(d/2), 1e8 <= d <= 1e100
    Rmax = 2000*iv.mpf(10)**50 + 1
    lamS = -(iv.log(iv.sqrt(2)) + 4*iv.log(2*Rmax+1) + lneps/4)
    show("-ln g_R at R_max", lamS)
    check(f"s={s}: g(R_max) <= 1/2", lo(lamS) >= hi(iv.log(2)))
    check(f"s={s}: lambda_S/8002 > 1/2", lo(lamS/8002) > 0.5)
    dS = iv.mpf('1e8')
    chainS = log10(iv.mpf(10)) - (lamS/8002 - iv.mpf(1)/2)*iv.sqrt(dS)/ln10
    mainS  = log10(iv.mpf(9)/4*bn) + (2*aG*(2*s+9) - (aG*iv.sqrt(dS) - iv.mpf(1)/2)*iv.sqrt(dS))/ln10
    farS   = log10(iv.mpf(17)/4*bn*s1) + L10e + aG*(4*s+17)/ln10
    show("log10 chain/e^{-sqrt d/2} at d=1e8", chainS)
    show("log10 main/e^{-sqrt d/2} at d=1e8", mainS)
    show("log10 far coefficient", farS)
    tot = iv.mpf(10)**hi(chainS) + iv.mpf(10)**hi(mainS) + iv.mpf(10)**hi(farS)
    show("log10 C_S", log10(tot))
    out[s].update(CS=hi(log10(tot)), chainS=hi(chainS))

print("\n=== regime S, d-independent conditions ===")
dS = iv.mpf('1e8')
check("d=1e8: d/2 - 1 >= 2000 sqrt d + 1 (t >= R)", lo(dS/2 - 1 - 2000*iv.sqrt(dS) - 1) >= 0)
check("d=1e8: d/2 >= sqrt d + 4001 (floor(t/R) >= sqrt d/8002)", lo(dS/2 - iv.sqrt(dS) - 4001) >= 0)
check("d=1e8: alpha_G > 1/(4 sqrt d) (main ratio decreases)", lo(aG - 1/(4*iv.sqrt(dS))) > 0)
import math
check("d<=1e100: ceil(2000 sqrt d) <= 2000*10^50 + 1 = R_max (exact integers)", math.isqrt(10**100) == 10**50)

print("\n=== stated numbers ===")
check("s=1: log10(b_* nu) <= 22.084", out[1]['bn'] <= 22.084)
check("s=5: log10(b_* nu) <= 23.104", out[5]['bn'] <= 23.104)
check("s=1: -ln g per block >= 736", out[1]['lamR'] >= 736)
check("s=5: -ln g per block >= 503", out[5]['lamR'] >= 503)
check("s=1: c_1 <= 26.26", out[1]['cL'] <= 26.26)
check("s=5: c_5 <= 27.29", out[5]['cL'] <= 27.29)
check("s=1: 2 Phi_R <= 10^-7011.8", out[1]['fL'] <= -7011.8)
check("s=5: 2 Phi_R <= 10^-8285.9", out[5]["fL"] <= -8285.9)
check("s=1: lambda_1 >= 0.0736", out[1]['lamL'] >= 0.0736)
check("s=5: lambda_5 >= 0.0503", out[5]['lamL'] >= 0.0503)
check("s=1: d_0 <= 1.21e5", out[1]['d0'] <= 1.21e5)
check("s=5: d_0 <= 1.26e5", out[5]['d0'] <= 1.26e5)
check("s=1: a_1 >= log10 C_1 (22.75)", out[1]['C1'] <= 22.75)
check("s=5: a_5 >= log10 C_1 (23.77)", out[5]['C1'] <= 23.77)
check("s=1: fibre floor <= 10^-14062.6", out[1]['aF'] <= -14062.6)
check("s=5: fibre floor <= 10^-16611.8", out[5]['aF'] <= -16611.8)
check("s=1: C_S <= 10^-1965", out[1]['CS'] <= -1965)
check("s=5: C_S <= 10^-2761", out[5]['CS'] <= -2761)
print("\nALL PASS" if all(v for _, v in ok) else "\nSOME FAIL")
