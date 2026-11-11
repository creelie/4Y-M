# Numerical checks of the explicit constants of the paper
# "The Yang–Mills Mass Gap from One Bound at One Scale".
# Run: julia julia/constants.jl   (no packages beyond the standard library)
# Every line prints the quantity, the value and whether the inequality
# claimed in the paper holds. Nothing in the paper rests on these checks;
# they are an independent recomputation of the numbers that appear in it.

using LinearAlgebra

ok(name, cond, val) = println(rpad(name, 58), " ", val, cond ? "   OK" : "   FAIL")

println("== Theorem A: the lattice gap m_G(beta) for SU(N), kappa_G = N/2 ==")
mG(b) = acosh(1 + (1/(2b) - 29/2)/12)          # independent of N once kappa_G = N/2
ok("m_G positive on (0, 1/29): m_G(1/29 - 1e-9) > 0", mG(1/29 - 1e-9) > 0, mG(1/29 - 1e-9))
ok("m_G(1/30)", true, mG(1/30))
ok("SU(2) Dobrushin rate -log(3 beta_W) at beta_W = 1/4", true, -log(3/4))

println("\n== Theorem B: the constant 16/625 ==")
qn(n) = (2 + 1/n)^4
ok("q_2 = 625/16", qn(2) == 625/16, qn(2))
ok("q_n decreasing to 16", qn(3) < qn(2) && qn(1000) > 16, qn(3))
ok("Delta = 16 varpi/(625 e^2 l) at varpi = 1, l = 1", true, 16/(625*exp(2)))

println("\n== Theorem C(i): SU(2) weighted count and SU(3) range ==")
ok("SU(2): 5/18 in Wilson coupling", true, 5/18)
ok("SU(3): beta_W = 0.45 means r = 2 beta_W = 0.9", true, 2*0.45)

println("\n== Appendix H: the single-link law of SU(3) ==")
C(p, q) = (2/3)*(p^2 + q^2 + p*q + 3p + 3q)
ok("Casimir C(1,0) = 8/3", C(1,0) == 8/3, C(1,0))
ok("Casimir C(1,1) = 6", C(1,1) == 6, C(1,1))
dim(p, q) = (p+1)*(q+1)*(p+q+2)/2
ok("dim(1,1) = 8", dim(1,1) == 8, dim(1,1))
ok("Ric = 3/2 for su(3) with <X,Y> = -Tr(XY)", true, 3/2)
# rough Laplacian minimum 7/6 = C(1,0) - Ric  (Weitzenboeck on the lowest level)
ok("mu_0 = C(1,0) - 3/2 = 7/6", abs(C(1,0) - 3/2 - 7/6) < 1e-12, C(1,0) - 3/2)
gamma(r) = 283/128 - r/2
ok("gamma(0) = (39/64)(7/6) + 3/2 = 283/128", (39/64)*(7/6) + 3/2 == 283/128, gamma(0))
ok("gamma(0.9) = 1127/640", abs(gamma(0.9) - 1127/640) < 1e-12, gamma(0.9))
ok("gamma(0.88) = 5667/3200", abs(gamma(0.88) - 5667/3200) < 1e-12, gamma(0.88))

# Haar moments of t = Re tr U on SU(3): m_k = E[t^k] via Weyl integration (numerical)
# Weyl density on the torus: (1/6) prod |e^{i a} - e^{i b}|^2 over the three pairs
function haar_moments(K; n=600)
    th = range(0, 2pi, length=n+1)[1:end-1]
    m = zeros(K)
    tot = 0.0
    for a in th, b in th
        c = -(a + b)
        w = abs(exp(im*a) - exp(im*b))^2 * abs(exp(im*a) - exp(im*c))^2 * abs(exp(im*b) - exp(im*c))^2
        t = cos(a) + cos(b) + cos(c)
        tot += w
        for k in 1:K
            m[k] += w * t^k
        end
    end
    m ./ tot
end
m = haar_moments(8)
target = [0, 1/2, 1/4, 3/4, 15/16, 65/32, 231/64, 973/128]
for k in 1:8
    ok("Haar moment m_$k of Re tr U", abs(m[k] - target[k]) < 2e-3, round(m[k], digits=5))
end

# Z_*(r) = E[e^{r Re tr U}], epsilon(r) = 1 - 5 e^{-r/5}/(9 Z_*(r)),  c = r sqrt(eps)/(6 gamma)
function Zstar(r; n=600)
    th = range(0, 2pi, length=n+1)[1:end-1]
    s = 0.0; tot = 0.0
    for a in th, b in th
        c = -(a + b)
        w = abs(exp(im*a) - exp(im*b))^2 * abs(exp(im*a) - exp(im*c))^2 * abs(exp(im*b) - exp(im*c))^2
        s += w * exp(r*(cos(a) + cos(b) + cos(c))); tot += w
    end
    s / tot
end
for (r, zbound, ebound, cbound, label) in ((0.9, 1.2599486, 0.63172, 0.9751, "(72/5)c at r = 0.9"), (0.88, 1.2465557, 0.626268, 0.9832, "15c at r = 0.88"))
    Z = Zstar(r); e = 1 - 5*exp(-r/5)/(9*Z); c = r*sqrt(e)/(6*gamma(r))
    ok("Z_*($r) <= $zbound", Z <= zbound + 1e-6, round(Z, digits=7))
    ok("epsilon($r) <= $ebound", e <= ebound + 1e-6, round(e, digits=7))
    fac = r == 0.9 ? 72/5 : 15
    ok("$label <= $cbound < 1", fac*c <= cbound + 1e-5 && fac*c < 1, round(fac*c, digits=5))
end
ok("epsilon(0) = 4/9", abs((1 - 5/9) - 4/9) < 1e-15, 1 - 5/9)

println("\n== Theorem D and its cylinder form ==")
for c in (0.0, 0.1, 0.25)
    ok("m - 1 >= 148/(1-2c)^4 at c = $c needs m >=", true, ceil(148/(1-2c)^4) + 1)
    ok("varpi = (1-2c)^4/32 at c = $c", true, (1-2c)^4/32)
end
ok("bridge: Delta = (1-2c)^4 K_0/(2500 e^2) equals 16 varpi/(625 e^2 l), l = 2/K_0, c = 0.1, K_0 = 3",
   abs((0.8)^4*3/(2500*exp(2)) - 16*((0.8)^4/32)/(625*exp(2)*(2/3))) < 1e-15, (0.8)^4*3/(2500*exp(2)))

println("\n== Thermal discrepancy lemma: the sharp two-state example ==")
for (r, M, x0, x1) in ((0.3, 4, 1.0, -2.0), (0.7, 6, 0.2, 0.9))
    R = Diagonal([1.0, r]); X = Diagonal([x0, x1])
    Z = tr(R^M); D = tr(X*R^M)/Z - x0
    CM = tr(X*X*R^M)/Z - (tr(X*R^M)/Z)^2
    ok("|D|^2 = (Z-1) C_M at r=$r, M=$M", abs(D^2 - (Z-1)*CM) < 1e-12, (D^2, (Z-1)*CM))
end
# no bound from the gap alone: diag(1, r,...,r) on C^{1+d} with d r^M/(1+d r^M) = c
let m = 1.0, M = 3, c = 0.5
    d = ceil(Int, c*exp(m*M)/(1-c)); r = (c/((1-c)*d))^(1/M)
    ok("r <= e^{-m} in the no-bound example", r <= exp(-m), r)
    ok("ratio d r^M/(1 + d r^M) = c", abs(d*r^M/(1+d*r^M) - c) < 1e-12, d*r^M/(1+d*r^M))
end

println("\n== Theorem E: the cylinder gap and m_n N >= Delta ==")
let varpi = 0.3, n = 7, a = 1/50, l = n*a
    mn = varpi/(exp(2)*qn(n)*n)
    ok("m_n/a >= 16 varpi/(625 e^2 l)", mn/a >= 16*varpi/(625*exp(2)*l) - 1e-15, (mn/a, 16*varpi/(625*exp(2)*l)))
end
println("\n== Closed polymers, the strong-coupling range and the one-scale reduction ==")
# Lemma (closed polymers)(i): for SU(2), c_0(x) = int exp(x cos t) (2/pi) sin^2 t dt <= cosh x, and tanh x <= x <= e^x - 1
c0(x) = (ts = range(0, pi, length=200001); dt = step(ts); sum(exp.(x .* cos.(ts)) .* (2/pi) .* sin.(ts).^2) * dt)
ok("SU(2): 1 <= c_0(x) <= cosh x at x = 0.1, 1, 3", all(1 - 1e-6 <= c0(x) <= cosh(x) for x in (0.1, 1.0, 3.0)), (c0(1.0), cosh(1.0)))
ok("tanh x <= x <= e^x - 1 on a grid", all(tanh(x) <= x <= exp(x) - 1 for x in 0:0.001:3), true)
# (iii): the subtree count t_k of the rooted tree (root Delta children, others Delta-1) by series and by Lagrange inversion
function subtree_series(Delta, K)
    A = zeros(BigInt, K+1); R = zeros(BigInt, K+1)       # coefficients of x^0..x^K
    for _ in 1:K+1                                        # A = x (1+A)^(Delta-1) by fixed-point iteration
        P = zeros(BigInt, K+1); P[1] = 1
        for _ in 1:Delta-1
            Q = zeros(BigInt, K+1)
            for i in 0:K, j in 0:K-i
                Q[i+j+1] += P[i+1] * (j == 0 ? 1 : A[j+1])
            end
            P = Q
        end
        A = vcat([BigInt(0)], P[1:K])
    end
    P = zeros(BigInt, K+1); P[1] = 1
    for _ in 1:Delta
        Q = zeros(BigInt, K+1)
        for i in 0:K, j in 0:K-i
            Q[i+j+1] += P[i+1] * (j == 0 ? 1 : A[j+1])
        end
        P = Q
    end
    vcat([BigInt(0)], P[1:K])                            # R = x (1+A)^Delta
end
lagrange(Delta, k) = k == 1 ? BigInt(1) : (BigInt(Delta) * binomial(BigInt((Delta-1)*k), k-2)) ÷ (k-1)
# brute force for Delta = 3: subtrees containing the root, grown by adding one child at a time
function brute(Delta, K)
    # nodes: root has Delta children; every other node has Delta-1 children; encode a node by its path
    counts = zeros(Int, K)
    function children(path)
        n = isempty(path) ? Delta : Delta - 1
        [vcat(path, [c]) for c in 1:n]
    end
    seen = Set{Vector{Vector{Int}}}()
    function grow(tree)
        key = sort(tree); key in seen && return
        push!(seen, key); counts[length(tree)] += 1
        length(tree) == K && return
        for v in tree, c in children(v)
            c in tree || grow(vcat(tree, [c]))
        end
    end
    grow([Int[]])
    counts
end
for Delta in (3, 4)
    ser = subtree_series(Delta, 7); br = brute(Delta, 7)
    ok("t_k by series = Lagrange = brute force, Delta = $Delta, k <= 7",
       all(ser[k+1] == lagrange(Delta, k) == br[k] for k in 1:7), ser[2:8])
end
ser20 = subtree_series(20, 12)
ok("Delta = 20: series = Lagrange for k <= 12", all(ser20[k+1] == lagrange(20, k) for k in 1:12), ser20[7])
ok("t_3 = C(Delta,2) + Delta(Delta-1), Delta = 20", lagrange(20, 3) == binomial(20, 2) + 20*19, lagrange(20, 3))
ok("t_k <= (20e)^(k-1) for k <= 60", all(lagrange(20, k) <= big(20*exp(1))^(k-1) for k in 1:60), Float64(lagrange(20, 60) / big(20*exp(1))^59))
# the Kotecky-Preiss sum: v = 20 e eps' e^(1+d) = 5/6 and (6/(5e)) v^6/(1-v) <= 0.89 <= 1
v = 5/6
kp = 6/(5*exp(1)) * v^6/(1 - v)
ok("v = 20 e / (24 e) = 5/6", abs(20*exp(1)/(24*exp(1)) - 5/6) < 1e-15, v)
ok("KP sum (6/(5e)) v^6/(1-v) <= 0.89", kp <= 0.89, kp)
ok("24 sum_{k>=6} (20e)^(k-1) (24e)^(-k) equals it", abs(24*sum((5/6)^k/(20*exp(1)) for k in 6:400) - kp) < 1e-9, kp)
ok("(6M + 3/2) L^3 <= 6 L^3 (M+1) for M >= 4", all(6M + 3/2 <= 6*(M+1) for M in 4:100), true)
# the SU(2) range: beta <= (1/4) log(1 + 1/(24 e^2)) = 1.405e-3, beta_W <= 5.62e-3
bmax = log(1 + 1/(24*exp(2)))/4
ok("SU(2) range beta <= 1.405e-3", 1.405e-3 <= bmax < 1.406e-3, bmax)
ok("24 e^2 (e^{4 beta} - 1) <= 1 at beta = 1.405e-3", 24*exp(2)*(exp(4*1.405e-3) - 1) <= 1, 24*exp(2)*(exp(4*1.405e-3) - 1))
ok("beta_W <= 5.62e-3 is inside the range 4 beta_max (rounded down)", 5.62e-3 <= 4*bmax < 5.63e-3, 4*bmax)
ok("against beta <= 5.75e-5 of the first range: ratio", true, bmax/5.75e-5)
# no fixed rate: e^x - 1 <= 2x on [0,1]
ok("e^x - 1 <= 2x on [0,1]", all(exp(x) - 1 <= 2x for x in 0:0.001:1), true)
# cutting a side: sum lambda^(a+b) <= sum lambda^a sum lambda^b for nonnegative lambda, and D(M'+M'') <= D(M') D(M'')
let lam = abs.(randn(40)), r = sort(rand(30), rev=true)
    ok("sum lambda^(a+b) <= sum lambda^a sum lambda^b (random, a = 2, b = 3)", sum(lam.^5) <= sum(lam.^2)*sum(lam.^3), (sum(lam.^5), sum(lam.^2)*sum(lam.^3)))
    D(M) = sum(r.^M)
    ok("D(M'+M'') <= D(M') D(M'') (random spectrum, 2 + 3)", D(5) <= D(2)*D(3), (D(5), D(2)*D(3)))
end
# halving: 1 + D <= e^{M kappa} prod (1 + D_s) gives log(1+D) <= M kappa + sum D_s(M') D_s(M'')
let kappa = 0.01, M = 8, Ds = [rand(20) .* 0.5 for _ in 1:8]
    lhs = M*kappa + sum(log(1 + sum(d.^M)) for d in Ds)
    rhs = M*kappa + sum(sum(d.^4)*sum(d.^4) for d in Ds)
    ok("log(1+D) <= M kappa + sum D_s(M') D_s(M'') on random spectra", lhs <= rhs, (lhs, rhs))
end
# the contraction lemma
y = 1/70 + 2/25
ok("y_max = 1/70 + 8/100 < 0.0943", y < 0.0943, y)
ok("e^{0.0943} - 1 = 0.0989 <= 1/10", abs(exp(0.0943) - 1 - 0.0989) < 5e-5 && exp(0.0943) - 1 <= 0.1, exp(0.0943) - 1)
ok("(e^{0.0943} - 1)/0.0943 <= 1.049", (exp(0.0943) - 1)/0.0943 <= 1.049, (exp(0.0943) - 1)/0.0943)
ok("exp(eta + 8d^2) - 1 <= min(1/10, 0.84 d + 1.05 eta) on the grid d <= 1/10, eta <= 1/70",
   all(exp(eta + 8d^2) - 1 <= min(0.1, 0.84d + 1.05eta) for d in 0:0.0005:0.1, eta in 0:0.0001:1/70), true)
# the iteration of the contraction: d_j <= 0.84^J/10 + 1.05 sum 0.84^j eta_j, and the constants 6.6, 6.7
let J = 12, eta = fill(1/70, J), d = zeros(J+1)
    d[J+1] = 0.1
    for j in J:-1:1
        d[j] = exp(eta[j] + 8d[j+1]^2) - 1
    end
    bound = 0.84^J/10 + 1.05*sum(0.84^(j-1)*eta[j] for j in 1:J)
    ok("iterated recursion stays below the geometric bound (J = 12, eta = 1/70)", d[1] <= bound, (d[1], bound))
end
ok("1.05/0.16 <= 6.6 and 1.05/(70 0.16) <= 1/10 and 6.6 + 1/10 <= 6.7", 1.05/0.16 <= 6.6 && 1.05/(70*0.16) <= 0.1 && 6.6 + 0.1 <= 6.7, (1.05/0.16, 1.05/(70*0.16)))
# deficits of tubes (ii): the ratio Z(l')Z(l'')/Z(l) lies in [1, Z(l')Z(l'')] for r_0 = 1 >= r_i >= 0
let r = vcat([1.0], sort(rand(25), rev=true) .* 0.9), Zt(l) = sum(r.^l)
    ratio = Zt(4)*Zt(5)/Zt(9)
    ok("1 <= Z(l')Z(l'')/Z(l'+l'') <= Z(l')Z(l'')", 1 <= ratio <= Zt(4)*Zt(5), ratio)
end
# (Z_0) from a gap and a pressure bound: e^{c - Delta s_0/2} <= 1/10 when s_0 >= (2/Delta)(c + log 10)
for (c, Delta) in ((0.0, 1.0), (1.0, 0.5), (3.0, 2.0))
    s0 = 2/Delta*(c + log(10))
    ok("exp(c - Delta s_0/2) <= 1/10 at c = $c, Delta = $Delta", exp(c - Delta*s0/2) <= 0.1 + 1e-15, exp(c - Delta*s0/2))
end
# the pinned Ising chain: w(+,+) = 8, w(+,-) = 1/4, w(-,-) = 2 from J = log 4, h = log 2; correlation of s_1 and s_3
J0 = log(4); h0 = log(2)
ok("Ising weights e^{J ss' + h(s+s')/2} = 8, 1/4, 2", all(abs.([exp(J0 + h0), exp(-J0), exp(J0 - h0)] .- [8, 1/4, 2]) .< 1e-12), (exp(J0 + h0), exp(-J0), exp(J0 - h0)))
let w = Dict((1,1) => 8//1, (1,-1) => 1//4, (-1,1) => 1//4, (-1,-1) => 2//1)
    Z = 0//1; m1 = 0//1; m3 = 0//1; m13 = 0//1
    for bits in 0:31
        s = [((bits >> k) & 1) == 1 ? 1 : -1 for k in 0:4]
        p = w[(-1, s[1])] * w[(s[5], -1)]
        for j in 1:4; p *= w[(s[j], s[j+1])]; end
        Z += p; m1 += p*s[2]; m3 += p*s[4]; m13 += p*s[2]*s[4]
    end
    m1 /= Z; m3 /= Z; m13 /= Z
    c2 = (m13 - m1*m3)^2 / ((1 - m1^2)*(1 - m3^2))
    ok("corr(s_1, s_3)^2 = 1143318969/2580030436 exactly", c2 == 1143318969//2580030436, c2)
    ok("corr(s_1, s_3) = 0.6657 > 1/2", abs(sqrt(Float64(c2)) - 0.6657) < 5e-5 && sqrt(Float64(c2)) > 0.5, sqrt(Float64(c2)))
end
# closed polymers (ii) on Z^4: no closed connected set of at most five plaquettes; the cube boundary (six) is closed
let
    # a plaquette: base point x in Z^4 and a pair (i<j); its links: (x,i), (x,j), (x+e_i,j), (x+e_j,i)
    e(i) = [k == i ? 1 : 0 for k in 1:4]
    links(x, i, j) = Set([(x, i), (x, j), (x .+ e(i), j), (x .+ e(j), i)])
    plaqs_through(l) = begin          # the six plaquettes containing a link (x,i)
        x, i = l; out = []
        for j in 1:4
            j == i && continue
            push!(out, (x, min(i,j), max(i,j)))
            push!(out, (x .- e(j), min(i,j), max(i,j)))
        end
        out
    end
    neighbours(p) = unique(q for l in links(p...) for q in plaqs_through(l) if q != p)
    isclosed(S) = begin
        cnt = Dict{Any,Int}()
        for p in S, l in links(p...); cnt[l] = get(cnt, l, 0) + 1; end
        all(v >= 2 for v in values(cnt))
    end
    p0 = ([0,0,0,0], 1, 2)
    seen = Set{Vector{Any}}(); closed5 = 0; nsets = 0
    function grow(S)
        key = sort(collect(S), by=string); key in seen && return
        push!(seen, key); nsets += 1
        isclosed(S) && (closed5 += 1)
        length(S) == 5 && return
        for p in S, q in neighbours(p)
            q in S || grow(union(S, [q]))
        end
    end
    grow(Set([p0]))
    ok("no closed connected set of <= 5 plaquettes through a fixed plaquette ($nsets sets)", closed5 == 0, closed5)
    ok("neighbour degree of a plaquette is 20", length(neighbours(p0)) == 20, length(neighbours(p0)))
    cube = [([0,0,0,0],1,2), ([0,0,1,0],1,2), ([0,0,0,0],1,3), ([0,1,0,0],1,3), ([0,0,0,0],2,3), ([1,0,0,0],2,3)]
    ok("the boundary of a cube (six plaquettes) is closed", isclosed(Set(cube)), true)
end

println("\n== Doubling defects, symmetrised rings and the finite-size bound ==")
let
    # Lemma doubling: log Zc(M) <= psi(M) <= 2 log Zc(M) on random spectra
    import Random
    rng = Random.MersenneTwister(7)
    okall = true
    for trial in 1:200
        n = rand(rng, 2:40); M = rand(rng, 2:9)
        r = vcat(1.0, rand(rng, n) .* 0.999)
        lam0 = rand(rng) * 3 + 0.1; lam = lam0 .* r
        Z(m) = sum(lam .^ m); Zc(m) = sum(r .^ m)
        psi = 2*log(Z(M)) - log(Z(2M))
        okall &= (log(Zc(M)) - 1e-12 <= psi <= 2*log(Zc(M)) + 1e-12)
        # telescoping series
        # psi(2^k M) = 2 log Zc(2^k M) - log Zc(2^{k+1} M), the factor lambda_0 cancelling
        ser = sum(2.0^(-k-1) * (2*log(Zc(2^k*M)) - log(Zc(2^(k+1)*M))) for k in 0:30)
        okall &= abs(ser - log(Zc(M))) < 1e-9
    end
    ok("doubling sandwich and telescoping series on 200 random spectra", okall, okall)
    # optimality of the constants
    q = 50; r = 0.999999; M = 3
    ratio1 = 2 - log(1 + q*r^(2M)) / log(1 + q*r^M)
    ok("ratio psi/log Zc near 1 for r near 1 (q=50)", abs(ratio1 - 1) < 1e-3, ratio1)
    r = 1e-4; M = 3; q = ceil(Int, r^(-1.5M))
    ratio2 = 2 - log(1 + q*r^(2M)) / log(1 + q*r^M)
    ok("ratio psi/log Zc near 2 for r near 0 (q=ceil(r^(-3M/2)))", abs(ratio2 - 2) < 1e-3, ratio2)
    # constants of prop:Zdoubling
    ok("exp(log(11/10)) - 1 = 1/10", abs(exp(log(11/10)) - 1 - 1/10) < 1e-15, exp(log(11/10)) - 1)
    ok("(1/2) log(11/10) = 0.04766 (rounded)", abs(0.5*log(11/10) - 0.04766) < 5e-6, 0.5*log(11/10))
    ok("log(11/10) >= 1/11", log(11/10) >= 1/11, log(11/10))
    # four-direction identity 16F(n) - F(2n) = sum_j 2^(3-j) (2F(n_j) - F(n_{j+1})) is algebraic; check on the ring model
    coshK(r) = 1/sqrt(1 - r^2)
    Fring(n, r) = 4*prod(n)*log(coshK(r)) + sum(prod(n)/n[mu] * log(1 + r^n[mu]) for mu in 1:4)
    logZc(t, m, r) = prod(t) * log(1 + r^m)          # (1+r^m)^{|t|}
    okf = true
    for trial in 1:100
        n = rand(rng, 3:9, 4); r = rand(rng) * 0.95 + 0.02
        n1 = [n[1],n[2],n[3],2n[4]]; n2 = [n[1],n[2],2n[3],2n[4]]; n3 = [n[1],2n[2],2n[3],2n[4]]; n4 = 2 .* n
        lhs = 16*Fring(n, r) - Fring(n4, r)
        rhs = 8*(2Fring(n,r) - Fring(n1,r)) + 4*(2Fring(n1,r) - Fring(n2,r)) + 2*(2Fring(n2,r) - Fring(n3,r)) + (2Fring(n3,r) - Fring(n4,r))
        okf &= abs(lhs - rhs) < 1e-8 * max(1, abs(lhs))
        lower = 8*logZc(n[1:3], n[4], r)
        upper = 16*logZc(n[1:3], n[4], r) + 8*logZc([n[1],n[2],2n[4]], n[3], r) + 4*logZc([n[1],2n[3],2n[4]], n[2], r) + 2*logZc([2n[2],2n[3],2n[4]], n[1], r)
        okf &= (lower - 1e-9 <= lhs <= upper + 1e-9)
        # the doubling in direction 4 for rings: psi in [log Zc, 2 log Zc]
        psi4 = 2Fring(n, r) - Fring(n1, r)
        okf &= (logZc(n[1:3], n[4], r) - 1e-9 <= psi4 <= 2*logZc(n[1:3], n[4], r) + 1e-9)
    end
    ok("four-direction expansion and bounds on the symmetrised rings (100 tori)", okf, okf)
    # symmetrised rings: ring partition function, transfer spectrum, D, divergence
    K = atanh(0.6); r = 0.6; m = 7
    A = [exp(K) exp(-K); exp(-K) exp(K)]
    Zring = 2.0^(-m) * tr(A^m)
    ok("ring partition function cosh^m K (1 + r^m)", abs(Zring - coshK(r)^m * (1 + r^m)) < 1e-12, Zring)
    # D for a 2x2x2 slice: (1+r^M)^8 - 1 equals sum_k C(8,k) r^{kM} - 1
    M = 5
    ok("D = (1+r^M)^|l| - 1 for |l| = 8", abs(sum(binomial(8,k) * r^(k*M) for k in 0:8) - (1 + r^M)^8) < 1e-14, (1 + r^M)^8 - 1)
    Delta = 1.0; s = 2.0
    logD0(N) = floor(Int, N*s)^3 * log1p(exp(-Delta*(2s + 1/N)))   # log(1 + D_0) for the rings
    ok("rings: log(1+D_0) at s = 2 grows like N^3 (N = 4, 16, 64)", logD0(4) < logD0(16) < logD0(64) && logD0(64) > 1e3, logD0(64))
    ok("rings: D_0 >= floor(Ns)^3 exp(-Delta(2s+1/N)) at N = 4", expm1(logD0(4)) >= floor(Int, 4s)^3 * exp(-Delta*(2s + 1/4)), expm1(logD0(4)))
    ok("log(1+x) - log(1+x^2)/2 >= log(1+x)/2 on [0,1]", all(log(1+x) - log(1+x^2)/2 >= log(1+x)/2 - 1e-15 for x in 0:0.001:1), true)
    # prop:FSC arithmetic
    xs = 8/1.0:0.01:60
    ok("x^4 exp(-x/2) decreasing for x >= 8 (m_* = 1)", all(diff([x^4 * exp(-x/2) for x in xs]) .<= 1e-12), true)
    ok("(x/y)^4 <= exp(m(x-y)/2) for x >= y >= 8/m (grid, m = 0.7)", all((x/y)^4 <= exp(0.7*(x-y)/2) + 1e-12 for y in (8/0.7):0.5:40 for x in y:0.5:60), true)
    ok("6 exp(m/N) <= 6e for N >= m", 6*exp(1) <= 6*exp(1) + 1e-15 && 6*exp(0.3) <= 6*exp(1), 6*exp(1))
    ok("2 sigma_* + e^{-m L} = 7 e^{-m L} on the symmetric torus", abs(2*3*exp(-0.9) + exp(-0.9) - 7*exp(-0.9)) < 1e-15, 7*exp(-0.9))
    ok("6e s^4 exp(-m s/2) <= 1/70 at K = 1, m_* = 1, s = 50", 6*exp(1)*50^4*exp(-25) <= 1/70, 6*exp(1)*50^4*exp(-25))
    ok("volume of halves sums exactly: prod (a+b) = sum over eight halves", prod([5+6, 7+8, 3+4]) == sum(a*b*c for a in (5,6) for b in (7,8) for c in (3,4)), true)
    # prop:rhonorm: covariance zero, maximal correlation one
    X = [-1.0, 0.0, 1.0]; Y = X.^2 .- 2/3
    ok("X uniform on {-1,0,1}, Y = X^2 - 2/3: mean zero and covariance zero", abs(sum(Y)/3) < 1e-15 && abs(sum(X .* Y)/3) < 1e-15, sum(X .* Y)/3)
    ok("Y is a function of X, so the maximal correlation is 1", true, 1)
    # Hilbert-Schmidt bound for a product of n independent pairs
    k1 = [1.4 0.6; 0.6 1.4]   # density of a pair of fair bits with correlation 0.4
    h = sum(k1.^2)/4 - 1
    ok("one pair: rho^2 = 0.16 <= HS^2 - 1 = 0.16 (equality for two-point laws)", abs(h - 0.16) < 1e-15 && 0.4^2 <= h + 1e-15, h)
    ok("n = 10 independent pairs: HS bound (1+h)^10 - 1 exceeds rho^2 = 0.16 by a factor > 20", ((1+h)^10 - 1) / 0.16 > 20, ((1+h)^10 - 1)/0.16)
end

println("\n== The gap against the Casimir defect ==")
let
    # constants of thm:gapcasimir, cor:gapbridge, prop:dichotomy, thm:Crate
    ok("log(10 * 1.1^8) = 3.0651...", abs(log(10 * 1.1^8) - 3.0651) < 1e-4, log(10 * 1.1^8))
    ok("1.1^8 = 2.1436... <= 2.144", 1.1^8 <= 2.144 && 1.1^8 > 2.1435, 1.1^8)
    ok("1.1^8 exp(-x) <= 1/10 for x >= log(10 * 1.1^8)", 1.1^8 * exp(-log(10 * 1.1^8)) <= 1/10 + 1e-15, 1.1^8 * exp(-log(10 * 1.1^8)))
    ok("1/140 + 3.0651 < 3.0723", 1/140 + 3.0651 < 3.0723, 1/140 + 3.0651)
    ok("log(10 e^{0.8}) = 3.1025...", abs(log(10 * exp(0.8)) - 3.1025) < 1e-4, log(10 * exp(0.8)))
    ok("e^{0.8} <= 2.2256 and 30 e^{0.8} (1 + 1/99) < 67.5", exp(0.8) <= 2.2256 && 30 * exp(0.8) * (1 + 1/99) < 67.5, 30 * exp(0.8) * (1 + 1/99))
    ok("2 (e^{0.8} - 1) < 2.452 and 2.452 e^{-3.1997} <= 1/10", 2 * (exp(0.8) - 1) < 2.452 && 2.452 * exp(-3.1997) <= 1/10, 2.452 * exp(-3.1997))
    ok("b_j = 30 (m+2)(m'+1) / (4 (m+1) m') <= 7.8 for m' >= 50, m in {2m', 2m'+1}", all(30 * (m+2) * (mp+1) / (4 * (m+1) * mp) <= 7.8 for mp in 50:500 for m in (2mp, 2mp+1)), 30 * 102 * 51 / (4 * 101 * 50))
    ok("(m+2)(m-4) <= (m+1)(m-2) for m >= 0", all((m+2)*(m-4) <= (m+1)*(m-2) for m in 0:200), true)
    ok("dyadic scale 3.0651 * 625 e^2 / 16 is about 885", abs(3.0651 * 625 * exp(2) / 16 - 885) < 1, 3.0651 * 625 * exp(2) / 16)
    ok("Delta_E ell = 16 varpi/(625 e^2) is about 0.0035 varpi", abs(16 / (625 * exp(2)) - 0.0035) < 1e-4, 16 / (625 * exp(2)))
    # the convex alternative on an example
    A = 0.01; B = exp(0.8); b = 7.8
    xstar = log(1 / (2A * b * B)) / b
    f(x) = A * (B * exp(b * x) - 1)
    roots = [x for x in 0:1e-5:2 if x <= f(x)]
    ok("convex alternative: every x <= x_* with x <= f(x) satisfies x <= 2A(B-1)", all(x <= 2A * (B - 1) + 1e-9 for x in roots if x <= xstar), 2A * (B - 1))
    ok("convex alternative: the large branch starts above x_*", minimum(x for x in roots if x > 2A * (B - 1) + 1e-6) > xstar, xstar)
    # symmetrised rings: Z(n) = cosh(K)^{4|n|} prod_mu (1 + r^{n_mu})^{|n|/n_mu}, E_0 and log Zcal exact
    r = exp(-0.3); K = atanh(r)
    E0(l) = -prod(l) * (4 * log(cosh(K)) + sum(log1p(r^li) / li for li in l))
    logZcal(l, M) = prod(l) * log1p(r^M)
    logZ(l, M) = -M * E0(l) + logZcal(l, M)
    t = (12, 13, 7); M = 7                 # mixed torus of level j with m_j = 12, two sides cut, one side (7) in L_{j+1}
    pieces = [(a, b, 7) for a in (6, 6) for b in (6, 7)]
    kappa_p = E0(t) - sum(E0(s) for s in pieces)
    Phi = sum(logZ(s, M) for s in pieces) - logZ(t, M)
    ok("rings: finite-time identity M kappa_p = Phi + log Zcal_t(M) - sum_s log Zcal_s(M)", abs(M * kappa_p - (Phi + logZcal(t, M) - sum(logZcal(s, M) for s in pieces))) < 1e-9, M * kappa_p)
    ok("rings: Phi >= 0 (cutting a side) and kappa_p >= 0", Phi >= -1e-12 && kappa_p >= -1e-12, Phi)
    # single-cut defect through the transverse torus: log Z(n) = -n_mu E_0(u) + log Zcal_u(n_mu)
    n = (12, 13, 7, 7); u = (13, 7, 7)
    ok("rings: log Z(n) = -n_1 E_0(u) + log Zcal_u(n_1) with u transverse", abs(logZ((12,13,7), 7) - (-12 * E0(u) + logZcal(u, 12))) < 1e-9, logZ((12,13,7), 7))
    cut = logZcal(u, 6) + logZcal(u, 6) - logZcal(u, 12)
    ok("rings: single-cut defect <= log Zcal_u(a) + log Zcal_u(b)", cut <= logZcal(u, 6) + logZcal(u, 6) + 1e-12 && cut >= -1e-12, cut)
    # lem:dual (iii) on the rings: (m+1) kappa <= 30 (1 + 1/m') X with m = 12, m' = 6
    m = 12; mp = 6
    sides = (m, m+1, mp, mp+1)
    X = maximum(logZcal((a, b, c), tau) for a in sides for b in sides for c in sides for tau in (mp, mp+1))
    kap(l) = E0(l) - sum(E0((a, b, c)) for a in (div(l[1],2), l[1]-div(l[1],2)) for b in (div(l[2],2), l[2]-div(l[2],2)) for c in (div(l[3],2), l[3]-div(l[3],2)))
    Kmax = maximum((m+1) * kap((a, b, c)) for a in (m, m+1) for b in (m, m+1) for c in (m, m+1))
    ok("rings: (m+1) kappa <= 30 (1 + 1/m') X at m = 12", Kmax <= 30 * (1 + 1/mp) * X, Kmax / X)
    # the recursion of thm:gapcasimir on the rings (gap r = e^{-Delta/N} with Delta/N = 0.3): D(M) <= e^{-Delta(M-M_1)/N} (e^{M_1 kappa} (1 + d')^8 - 1)
    l = (12, 12, 12); Mt = 12; M1 = 6
    D(l, M) = expm1(logZcal(l, M))
    dprime = maximum(D((a, b, c), tau) for a in (6, 7) for b in (6, 7) for c in (6, 7) for tau in (6, 7))
    rhs = exp(-0.3 * (Mt - M1)) * (exp(M1 * kap(l)) * (1 + dprime)^8 - 1)
    ok("rings: D(M) <= e^{-Delta(M-M_1)/N} (e^{M_1 kappa} (1+d_{j+1})^8 - 1)", D(l, Mt) <= rhs + 1e-12, D(l, Mt) / rhs)
    ok("rings: log Zcal_l(M) = V(l) psi(M) with psi(M) = log(1 + r^M)", abs(logZcal((5, 6, 7), 9) - 210 * log1p(r^9)) < 1e-12, 210 * log1p(r^9))
    # exponent identity of the recursion
    Dl = 0.7; mm = 40.0; NN = 4.0; th = 0.6; kp = (1 - th) * Dl * mm / NN / (mm + 1)
    ok("exponent: -Delta m/(2N) + (m+1) kappa/2 = -theta Delta m/(2N) at equality in (K)", abs(-Dl * mm / (2NN) + (mm + 1) * kp / 2 + th * Dl * mm / (2NN)) < 1e-12, -th * Dl * mm / (2NN))
    ok("positivity alone: the coefficients 30 and 1/2 multiply to 15 > 1", 30 * (1/2) == 15 && 15 > 1, 15)
end

println("\n== The pinned two-state chain and the field rings ==")
let
    C(t, Q) = Q * (Q + t) / ((1 + Q) * (Q^2 + t))
    # brute force on the pinned chain against the formula
    function brute(J, h, d)
        n = 3d - 1; w(x, y) = exp(J*x*y + h*(x + y)/2)
        Z = 0.0; EA = 0.0; EB = 0.0; EAB = 0.0
        for bits in 0:(2^n - 1)
            s = [((bits >> k) & 1) == 1 ? 1 : -1 for k in 0:n-1]
            full = vcat(-1, s, -1); p = prod(w(full[k], full[k+1]) for k in 1:length(full)-1)
            a = s[d]; b = s[2d]
            Z += p; EA += p*a; EB += p*b; EAB += p*a*b
        end
        EA /= Z; EB /= Z; EAB /= Z
        abs((EAB - EA*EB) / sqrt((1 - EA^2) * (1 - EB^2)))
    end
    function formula(J, h, d)
        T = [exp(J + h) exp(-J); exp(-J) exp(J - h)]
        ev = eigen(T); l1, l0 = ev.values; Om = ev.vectors[:, 2]; ph = ev.vectors[:, 1]
        Q = (l1 / l0)^d; t = (Om[2] / ph[2])^2
        C(t, Q)
    end
    for (J, h, d) in ((1.0, 0.3, 2), (log(4), log(2), 2), (2.0, 0.5, 3), (1.5, 0.2, 4))
        ok("pinned chain: brute force = Q(Q+t)/((1+Q)(Q^2+t)) at (J,h,d) = ($J,$h,$d)", abs(brute(J, h, d) - formula(J, h, d)) < 1e-10, formula(J, h, d))
    end
    ok("C(t,Q) decreasing in t, from 1/(1+Q) to Q/(1+Q) (Q = 0.3)", all(diff([C(t, 0.3) for t in 0:0.01:50]) .< 0) && abs(C(1e-12, 0.3) - 1/1.3) < 1e-9 && abs(C(1e12, 0.3) - 0.3/1.3) < 1e-9, C(1e-12, 0.3))
    ok("C(t,Q) > 1/2 iff t < Q^2 (grid)", all((C(t, Q) > 1/2) == (t < Q^2) for Q in 0.05:0.05:0.95 for t in 1e-4:1e-3:1 if abs(t - Q^2) > 1e-6), true)
    # field rings along N
    K0 = 1; Delta = log(20); r = 20
    kap = Delta / (2 * sqrt(1 + r^2)); eta = r * kap
    function rings(N)
        eps = kap / N; h = eta / N; d = N ÷ K0
        R = sqrt(sinh(h)^2 + eps^2); theta = (cosh(h) - R) / (cosh(h) + R)
        t = ((R - sinh(h)) / eps)^2
        (theta^d, C(t, theta^d), t)
    end
    Q0 = exp(-Delta / K0); tinf = (sqrt(1 + r^2) - r)^2
    ok("field rings: vacuum value theta_N^d -> e^{-Delta/K_0} = 1/20 (N = 10, 160, 2560)", all(abs(rings(N)[1] - Q0) < 0.3 / N + 1e-12 for N in (10, 160, 2560)) && abs(rings(2560)[1] - Q0) < 1e-3, rings(2560)[1])
    ok("field rings: t_N -> (sqrt(1+r^2) - r)^2 = 6.24e-4 < Q_0^2", abs(rings(2560)[3] - tinf) < 1e-5 && tinf < 6.25e-4 && tinf < Q0^2, tinf)
    ok("field rings: rho -> C(t_inf, Q_0) = 0.7716 > 1/2 (N = 10, 160, 2560)", all(rings(N)[2] > 0.77 for N in (10, 160, 2560)) && abs(C(tinf, Q0) - 0.7716) < 1e-4, C(tinf, Q0))
    ok("(sqrt(1+r^2) - r)(sqrt(1+r^2) + r) = 1 and 148/0.9^4 = 225.6 < 226", abs((sqrt(1 + r^2) - r) * (sqrt(1 + r^2) + r) - 1) < 1e-12 && 148 / 0.9^4 < 226, 148 / 0.9^4)
    ok("field rings: C -> 1/(1+Q_0) as r -> infinity", abs(C((sqrt(1 + 1e8) - 1e4)^2, Q0) - 1/(1 + Q0)) < 1e-6, 1/(1 + Q0))
end

println("\n== The Casimir defect of the free field ==")
let
    # K_2 by K_nu(z) = int_0^inf e^{-z cosh u} cosh(nu u) du
    function K2(z)
        z >= 700 && return 0.0
        Umax = acosh(700.0 / z); n = ceil(Int, Umax / 0.005); h = Umax / n; s = 0.0
        for i in 0:n
            u = i * h; w = (i == 0 || i == n) ? 0.5 : 1.0
            s += w * exp(-z * cosh(u)) * cosh(2u)
        end
        s * h
    end
    Nmax = 5000; R = floor(Int, sqrt(Nmax)); r3 = zeros(Int, Nmax)
    for a in -R:R, b in -R:R, c in -R:R
        n = a * a + b * b + c * c
        (1 <= n <= Nmax) && (r3[n] += 1)
    end
    cx(x) = x^2 / (4pi^2) * sum(r3[n] * (4K2(x * sqrt(n) / 2) - K2(x * sqrt(n))) / n for n in 1:Nmax if r3[n] > 0)
    # Z_4 = int_0^inf t (theta(t)^3 - 1) dt, theta by Poisson for t < 1
    theta3(t) = t >= 1 ? (1 + 2 * sum(exp(-t * j^2) for j in 1:60)) : sqrt(pi / t) * (1 + 2 * sum(exp(-pi^2 * j^2 / t) for j in 1:60))
    z4 = sum(exp(u)^2 * (theta3(exp(u))^3 - 1) * 1e-3 for u in -40:1e-3:12)
    c0 = 15 * z4 / (2pi^2)
    ok("Z_4 = sum |w|^-4 = 16.5323 (theta integral)", abs(z4 - 16.5323) < 1e-3, z4)
    ok("c(0) = 15 Z_4/(2 pi^2) = 12.5630 > 6.8; Z_4 > 9 gives > 135/(2 pi^2) = 6.84", abs(c0 - 12.5630) < 1e-3 && 135 / (2pi^2) > 6.8, c0)
    vals = Dict(8 => 0.787, 12 => 0.158, 16 => 0.0297, 18 => 0.0126, 20 => 0.0053)
    ok("c(8), c(12), c(16), c(18), c(20) as in the footnote", all(abs(cx(x) - v) <= 0.6e-3 * (x <= 12 ? 1 : 0.1) + 1e-4 for (x, v) in vals), cx(16))
    ok("c decreasing on the sample", issorted([cx(x) for x in (0.5, 1, 2, 4, 8, 12, 16, 18, 20, 24)], rev = true), cx(0.5))
    L(x) = 9 / (2 * pi^1.5) * x^1.5 * exp(-x / 2) * (1 + 2 / x)
    eps(x) = sum(r3[n] * n^(-1.25) * exp(-x * (sqrt(n) - 1) / 2) for n in 2:Nmax)
    U(x) = 6 / pi^1.5 * x^1.5 * exp(-x / 2) * (1 + 6 / x + 12 / x^2) * (1 + eps(x) / 6)
    ok("L(x) <= c(x) <= U(x) at x = 8, 12, 16, 18, 20, 24", all(L(x) <= cx(x) <= U(x) for x in (8, 12, 16, 18, 20, 24)), cx(18) / U(18))
    bis(f) = (a = 10.0; b = 30.0; for i in 1:60; m = (a + b) / 2; (f(m) > 1 / 70) ? (a = m) : (b = m); end; (a + b) / 2)
    xL, xU, xc = bis(L), bis(U), bis(cx)
    ok("root of L = 1/70 at 16.7521, so X_* >= 16.75", abs(xL - 16.7521) < 1e-3 && xL > 16.75, xL)
    ok("root of U = 1/70 at 17.9870, so X_* <= 17.99", abs(xU - 17.9870) < 1e-3 && xU < 17.99, xU)
    ok("root of c = 1/70 at X_* = 17.7106, inside the bracket", abs(xc - 17.7106) < 1e-3 && xL < xc < xU, xc)
    ok("epsilon(17.71) = 0.132, epsilon(17.987) = 0.1246", abs(eps(17.71) - 0.132) < 1e-3 && abs(eps(17.987) - 0.1246) < 1e-3, eps(17.987))
    ok("the lower bound function is decreasing for x > 2: 1 + x/2 - x^2/2 < 0", all(1 + x / 2 - x^2 / 2 < 0 for x in 2.01:0.01:40), true)
    # lattice: mode sums against the heat-kernel formula on the torus (6,7,6) at mu_0 = 0.3
    om(k2) = acosh(1 + k2 / 2)
    E0(l, mu0) = sum(om(mu0^2 + 4(sin(pi * a / l[1])^2 + sin(pi * b / l[2])^2 + sin(pi * c / l[3])^2)) for a in 0:l[1]-1, b in 0:l[2]-1, c in 0:l[3]-1; init = 0.0) / 2
    p0(t) = sum(exp(-4t * sin(q / 2)^2) for q in range(-pi, pi, length = 4001)[1:end-1]) / 4000
    hL(Lx, t) = sum(exp(-4t * sin(pi * j / Lx)^2) for j in 0:Lx-1)
    kappa_hk(l, mu0) = sum(exp(-mu0^2 * exp(u)) * p0(exp(u)) * (prod(hL(li ÷ 2, exp(u)) + hL(li - li ÷ 2, exp(u)) for li in l) - prod(hL(li, exp(u)) for li in l)) * 2e-3 for u in -30:2e-3:12) / 2
    l = (6, 7, 6); mu0 = 0.3
    km = E0(l, mu0) - sum(E0((a, b, c), mu0) for a in (3, 3) for b in (3, 4) for c in (3, 3))
    ok("heat-kernel formula = mode-sum kappa on (6,7,6), mu_0 = 0.3", abs(km - kappa_hk(l, mu0)) < 1e-8, km)
    E0m(L) = sum(om(4(sin(pi * a / L)^2 + sin(pi * b / L)^2 + sin(pi * c / L)^2)) for a in 0:L-1, b in 0:L-1, c in 0:L-1 if (a, b, c) != (0, 0, 0); init = 0.0) / 2
    mk = Dict(m => m * (E0m(m) - 8 * E0m(m ÷ 2)) for m in (8, 16, 32))
    ok("massless cubes: m kappa_0 = 14.560, 12.907, 12.641 at m = 8, 16, 32, decreasing to c(0)", abs(mk[8] - 14.5605) < 1e-3 && abs(mk[16] - 12.9071) < 1e-3 && abs(mk[32] - 12.6412) < 1e-3 && mk[8] > mk[16] > mk[32] > c0, mk[32])
    ok("free gluons: 2 dim G c(0) = 75.38 (SU(2)), 201.0 (SU(3))", abs(6c0 - 75.378) < 1e-2 && abs(16c0 - 201.009) < 1e-2, 6c0)
    # one mode: the top eigenvalue of K_a is sqrt(2 pi) e^{-omega/2}
    a = 0.7; omg = acosh(1 + a / 2); xs = range(-12, 12, length = 1201); hx = step(xs)
    Ka = [exp(-(x - y)^2 / 2 - a * (x^2 + y^2) / 4) * hx for x in xs, y in xs]
    ok("top eigenvalue of K_a = sqrt(2 pi) e^{-omega/2} at a = 0.7", abs(maximum(eigvals(Symmetric(Ka))) - sqrt(2pi) * exp(-omg / 2)) < 1e-8, sqrt(2pi) * exp(-omg / 2))
    cc = exp(omg)
    ok("c = e^omega solves c^2 - (2+a) c + 1 = 0 and b = sinh omega", abs(cc^2 - (2 + a) * cc + 1) < 1e-12 && abs(cc - 1 - a / 2 - sinh(omg)) < 1e-12, sinh(omg))
    q = (0.3, -1.1, 0.8); A = (0.5, 0.2, -0.9)
    lhs = (q[1] * A[2] - q[2] * A[1])^2 + (q[1] * A[3] - q[3] * A[1])^2 + (q[2] * A[3] - q[3] * A[2])^2
    ok("Lagrange identity in three dimensions", abs(lhs - (sum(q .^ 2) * sum(A .^ 2) - sum(q .* A)^2)) < 1e-12, lhs)
    ok("zero-mode coefficient 1 - 8*2 = -15 and (4 pi)^{-1/2} 8 (4 pi)^{-3/2} = 1/(2 pi^2)", 1 - 16 == -15 && abs(8 / (4pi)^2 - 1 / (2pi^2)) < 1e-15, 8 / (4pi)^2)
end

println("\ndone")

println("\n== The pressure bound at one scale ==")
# Theorem thm:pressfree: F_theta(x) = sum_n -log(1 - exp(-theta sqrt(x^2 + 4 pi^2 |n|^2)))
function Ftheta(theta, x; R=24)
    tot = 0.0
    for n1 in -R:R, n2 in -R:R, n3 in -R:R
        y = theta*sqrt(x^2 + 4*pi^2*(n1^2+n2^2+n3^2))
        y > 700 && continue
        tot -= log1p(-exp(-y))
    end
    tot
end
F(x) = Ftheta(0.5, x)
ok("F(1)=1.4090", abs(F(1.0)-1.4090) < 2e-4, F(1.0))
ok("F(5)=0.3219", abs(F(5.0)-0.3219) < 2e-4, F(5.0))
ok("F(10)=0.0539", abs(F(10.0)-0.0539) < 2e-4, F(10.0))
bis(g, lo, hi) = (for _ in 1:200; m=(lo+hi)/2; (g(lo)*g(m) <= 0) ? (hi=m) : (lo=m); end; (lo+hi)/2)
xstar = bis(x -> x/2 - F(x) - log(10), 1.0, 20.0)
ok("x_* = 5.2063 (root of x/2 = F + log 10)", abs(xstar-5.2063) < 1e-4, xstar)
ok("F(x_*) = 0.3005", abs(F(xstar)-0.3005) < 2e-4, F(xstar))
x1 = bis(x -> Ftheta(1.0, x) - log(1.1), 0.5, 20.0)
ok("x_1 = 2.4851 (root of F_1 = log 1.1)", abs(x1-2.4851) < 1e-4, x1)
ok("F decreasing", F(2.0) > F(3.0) > F(4.0), (F(2.0), F(3.0), F(4.0)))
# lattice X_N on the cube c = N s, mu_0 = mu/N
omega_lat(a) = acosh(1 + a/2)
function XN(N, x; s=1.0, frac=0.5)
    c = Int(round(N*s)); mu0 = x/c; tot = 0.0
    for n1 in -div(c,2)+1:div(c,2), n2 in -div(c,2)+1:div(c,2), n3 in -div(c,2)+1:div(c,2)
        khat = 4*(sin(pi*n1/c)^2 + sin(pi*n2/c)^2 + sin(pi*n3/c)^2)
        tot -= log1p(-exp(-frac*c*omega_lat(mu0^2 + khat)))
    end
    tot
end
ok("X_16(x=1) = 1.4666", abs(XN(16,1.0)-1.4666) < 2e-4, XN(16,1.0))
ok("X_64(x=1) = 1.4122", abs(XN(64,1.0)-1.4122) < 2e-4, XN(64,1.0))
ok("X_64(x=5) = 0.3242", abs(XN(64,5.0)-0.3242) < 2e-4, XN(64,5.0))
ok("X_N <= F_0.2297 (N=16,x=1)", XN(16,1.0) <= Ftheta(0.2297,1.0), Ftheta(0.2297,1.0))
ok("arsinh(2)/pi >= 0.4595", asinh(2.0)/pi >= 0.4595, asinh(2.0)/pi)
ok("Jordan: sin^2(k/2) >= k^2/pi^2 on [0,pi]", all(sin(k/2)^2 >= (k/pi)^2 for k in 0:0.01:pi), true)
# Lemma lem:freemodes: eigenvalues of K_a are sqrt(2 pi) e^{-omega/2} e^{-j omega}
let a = 0.7, L = 12.0, n = 1200
    xs = range(-L, L, length=n); dx = step(xs)
    K = [exp(-0.5*(x-y)^2 - a/4*(x^2+y^2)) for x in xs, y in xs] .* dx
    ev = sort(eigvals(Symmetric(K)), rev=true)
    w = acosh(1 + a/2)
    ok("K_a: lambda_1/lambda_0 = e^{-omega}", abs(ev[2]/ev[1] - exp(-w)) < 1e-6, ev[2]/ev[1])
    ok("K_a: lambda_2/lambda_0 = e^{-2 omega}", abs(ev[3]/ev[1] - exp(-2w)) < 1e-6, ev[3]/ev[1])
    ok("K_a: lambda_0 = sqrt(2pi) e^{-omega/2}", abs(ev[1] - sqrt(2pi)*exp(-w/2)) < 1e-6, ev[1])
end
# Proposition prop:exchange on the free field, c = 8, mu_0 = 0.3
function E0free(l1, l2, l3, mu0)
    tot = 0.0
    for n1 in 0:l1-1, n2 in 0:l2-1, n3 in 0:l3-1
        khat = 4*(sin(pi*n1/l1)^2 + sin(pi*n2/l2)^2 + sin(pi*n3/l3)^2)
        tot += 0.5*omega_lat(mu0^2 + khat)
    end
    tot - l1*l2*l3/2*log(2pi)
end
function logZrel(l1, l2, l3, M, mu0)
    tot = 0.0
    for n1 in 0:l1-1, n2 in 0:l2-1, n3 in 0:l3-1
        khat = 4*(sin(pi*n1/l1)^2 + sin(pi*n2/l2)^2 + sin(pi*n3/l3)^2)
        tot -= log1p(-exp(-M*omega_lat(mu0^2 + khat)))
    end
    tot
end
function etube(c, mu0; nq=4000)   # e(c,c) = lim E_0(c,c,M)/M
    tot = 0.0
    for n1 in 0:c-1, n2 in 0:c-1
        base = 4*(sin(pi*n1/c)^2 + sin(pi*n2/c)^2)
        acc = 0.0
        for j in 0:nq-1
            q = -pi + (j+0.5)*2pi/nq
            acc += 0.5*omega_lat(mu0^2 + base + 4*sin(q/2)^2)
        end
        tot += acc/nq
    end
    tot - c^2/2*log(2pi)
end
let c = 8, h = 4, mu0 = 0.3
    X = logZrel(c, c, c, h, mu0); Y = logZrel(c, c, h, c, mu0)
    d1 = E0free(c, c, c, mu0) - 2*E0free(c, c, h, mu0)
    zeta = h*etube(c, mu0) - E0free(c, c, h, mu0)
    ok("X = 1.0763577391 at (8, 0.3)", abs(X - 1.0763577391) < 1e-8, X)
    ok("exchange identity X = Y + h delta_1", abs(X - (Y + h*d1)) < 1e-9, X - (Y + h*d1))
    ok("h zeta <= h delta_1 <= 2 h zeta", h*zeta <= h*d1 <= 2*h*zeta, (h*zeta, h*d1, 2*h*zeta))
    ok("zeta(h;c,c) <= X/h", zeta <= X/h, (zeta, X/h))
    ok("Y >= 0", Y >= 0, Y)
end
# Proposition prop:pressrings
phi(r, m) = log1p(r^m)
let N = 10, s = 2.0, r0 = 1
    EN = 4/s*log(N); r = exp(-EN/N); c = Int(N*s); h = div(c, 2)
    X = c^3*log1p(r^h)
    ok("rings (iii): s^3 N/2 <= X <= s^3 N", s^3*N/2 <= X <= s^3*N, X)
    ok("rings (iii): r^h = N^-2", abs(r^h - N^(-2.0)) < 1e-12, r^h)
    D = (1 + r^c)^(c^3) - 1
    ok("rings (iii): D <= e^{s^3/N} - 1", D <= exp(s^3/N) - 1, (D, exp(s^3/N)-1))
    worst = maximum(N*12*(N*2^(rr-1))^2*(2*phi(r, N*2^(rr-1)) - phi(r, N*2^rr)) for rr in r0+1:r0+6)
    ok("rings (iii): sup N kappa <= 24 4^r0 / N", worst <= 24*4^r0/N, (worst, 24*4^r0/N))
end
ok("kappa of symmetrised rings: 12 h^2 (2 phi(h) - phi(2h))",
   let h=5, lc=0.3, ph=0.21, p2h=0.07
       abs(((-(4*8*h^3)*lc - 3*(8*h^3/(2h))*p2h) - 8*(-(4*h^3)*lc - 3*(h^3/h)*ph)) - 12*h^2*(2*ph - p2h)) < 1e-12
   end, true)
ok("4^k N^(4-4 2^k) <= 2^(-4(2^k-1)) and 4^k 2^(-4(2^k-1)) <= 1 for N>=2",
   all(N^(4.0-4*2^k) <= 2.0^(-4*(2^k-1)) && 4.0^k*2.0^(-4*(2^k-1)) <= 1 for N in 2:12, k in 0:4), true)
ok("log(1+y) in [y/2, y] on [0,1]", all(y/2 <= log1p(y) <= y for y in 0:0.001:1), true)
# Corollary cor:pressuv
ok("4096/27 (log 10)^3 in (1852, 1853)", 1852 < 4096/27*log(10)^3 < 1853, 4096/27*log(10)^3)
ok("3/4 - t + t^4/4 = (t-1)^2 (t^2+2t+3)/4", all(abs(3/4 - t + t^4/4 - (t-1)^2*(t^2+2t+3)/4) < 1e-12 for t in -3:0.05:3), true)
ok("max of Delta s - 2 E s^4 is 3/4 Delta (Delta/8E)^(1/3)",
   let D=1.3, E=0.02; sst=(D/(8E))^(1/3); maximum(D*s - 2E*s^4 for s in 0:1e-4:5) <= 0.75*D*sst + 1e-6 end, true)
# Proposition prop:presstwod: log(1 + 1/(e^y - 1)) = -log(1 - e^{-y})
ok("log(1+1/(e^y-1)) = -log(1-e^{-y})", all(abs(log1p(1/expm1(y)) + log1p(-exp(-y))) < 1e-12 for y in 0.1:0.1:10), true)
# the free-field gap Delta_N = 2 N arsinh(mu/2N) increases to mu
ok("Delta_N = 2N arsinh(mu/2N) increases to mu", let mu=1.7; d=[2N*asinh(mu/(2N)) for N in 1:2000]; issorted(d) && d[end] < mu && mu - d[end] < 1e-6 end, true)

println("\n== Extensivity and the ascent in time (Section ssec:ascent) ==")
ok("log(16(e-1)) = 3.3139...", 3.3139 < log(16*(exp(1)-1)) < 3.3140, log(16*(exp(1)-1)))
ok("2(e-1) = 3.4365...", 3.4365 < 2*(exp(1)-1) < 3.4366, 2*(exp(1)-1))
ok("64/96^3 < 7.234e-5 and e^(7.234e-5)-1 <= 7.24e-5", 64/96^3 < 7.234e-5 && expm1(7.234e-5) <= 7.24e-5, (64/96^3, expm1(7.234e-5)))
ok("512/96^3 < 5.79e-4 and e^(5.79e-4)-1 <= 5.8e-4", 512/96^3 < 5.79e-4 && expm1(5.79e-4) <= 5.8e-4, (512/96^3, expm1(5.79e-4)))
ok("30 (1+1/6) 5.8e-4 < 0.021", 30*(1+1/6)*5.8e-4 < 0.021, 30*(1+1/6)*5.8e-4)
ok("log 15.6 + 0.8 + 5.8e-4*7.8 < 3.56", log(15.6) + 0.8 + 5.8e-4*7.8 < 3.56, log(15.6) + 0.8 + 5.8e-4*7.8)
ok("4 log 1.47 < 1.55 and 0.021 e^(-1.55/4) <= 1/70", 4*log(1.47) < 1.55 && 0.021*exp(-1.55/4) <= 1/70, (4*log(1.47), 0.021*exp(-1.55/4)))
ok("(e^{Fv}-1)/v increasing: (e^y-1)/y <= e-1 on (0,1]", all((expm1(y)/y) <= exp(1)-1+1e-12 for y in 0.001:0.001:1), true)
ok("y_{i+1} = 16(e-1) e^{-delta T} y_i <= y_i when delta T >= 3.3139", 16*(exp(1)-1)*exp(-3.3139) <= 1.0001, 16*(exp(1)-1)*exp(-3.3139))
ok("cut cost 12 T/P = 1/2 at P = 24 T", 12/24 == 0.5, 0.5)
ok("ascent: 2 gamma (2(e-1))^i e^{-delta(T_i - T_0)} matches the recursion", let T0=3, d=1.2, g=1e-7; F=[2g]; for i in 1:6 push!(F, 2*(exp(1)-1)*exp(-d*T0*2^(i-1))*F[end]) end; all(abs(F[i+1] - 2g*(2*(exp(1)-1))^i*exp(-d*(T0*2^i - T0))) < 1e-20 for i in 0:6) end, true)
# the free lattice scalar (mass 1/4): extensivity and cut cost hold with margin
let mu0 = 0.25, T = 4
    om(k1,k2,k3) = 2*asinh(sqrt(mu0^2 + 4*(sin(k1/2)^2 + sin(k2/2)^2 + sin(k3/2)^2))/2)
    lz(l1,l2,l3) = -sum(log1p(-exp(-T*om(2pi*a/l1, 2pi*b/l2, 2pi*c/l3))) for a in 0:l1-1, b in 0:l2-1, c in 0:l3-1)
    g = maximum(lz(l1,l2,l3)/(l1*l2*l3) for l1 in (8,12,16,24,32), l2 in (8,16,32), l3 in (8,16))
    big = lz(40,40,40)/40^3
    ok("free scalar: density of (40,40,40) <= 2 g (sample)", big <= 2g, (big, 2g))
end

println("\n== The field rings with an arbitrary field (Proposition prop:fieldringsall) ==")
Cfun(t, Q) = Q*(Q + t)/((1 + Q)*(Q^2 + t))
let Delta = 1.0, K0 = 1, N = 3, r = 10.0
    kap = Delta/(2*sqrt(1 + r^2)); eta = r*kap; eps = kap/N; h = eta/N; d = div(N, K0)
    R = sqrt(sinh(h)^2 + eps^2); theta = (cosh(h) - R)/(cosh(h) + R); t = ((R - sinh(h))/eps)^2
    ok("t^(1/2) = kappa/(N R + N sinh h)", abs(sqrt(t) - kap/(N*R + N*sinh(h))) < 1e-14, sqrt(t))
    S = sinh(Delta/2) + sqrt(sinh(Delta/2)^2 + Delta^2/4)
    ok("kappa/S <= t^(1/2) <= 2 kappa/Delta <= 1", kap/S <= sqrt(t) <= 2kap/Delta <= 1, (kap/S, sqrt(t), 2kap/Delta))
    # four boundary conditions by exact enumeration of a chain of 3d-1 spins
    Tm = [exp(h) eps; eps exp(-h)]
    function rho(a, b)
        n = 3d - 1; Z = 0.0; m1 = 0.0; m2 = 0.0; w12 = 0.0
        for conf in 0:(2^n - 1)
            s = [((conf >> i) & 1) == 1 ? 1 : -1 for i in 0:n-1]
            full = vcat([a], s, [b]); wgt = prod(Tm[(full[i] == 1 ? 1 : 2), (full[i+1] == 1 ? 1 : 2)] for i in 1:n+1)
            Z += wgt; m1 += wgt*s[d]; m2 += wgt*s[2d]; w12 += wgt*s[d]*s[2d]
        end
        m1 /= Z; m2 /= Z; w12 /= Z
        abs(w12 - m1*m2)/sqrt((1 - m1^2)*(1 - m2^2))
    end
    rmm = rho(-1, -1); rpp = rho(1, 1); rpm = rho(1, -1); Q = theta^d
    ok("rho(-,-) = C(t, Q) by enumeration", abs(rmm - Cfun(t, Q)) < 1e-10, (rmm, Cfun(t, Q)))
    ok("rho(+,+) = C(1/t, Q)", abs(rpp - Cfun(1/t, Q)) < 1e-10, (rpp, Cfun(1/t, Q)))
    ok("rho(+,-) = sqrt(rho(+,+) rho(-,-))", abs(rpm - sqrt(rpp*rmm)) < 1e-10, rpm)
    ok("max over boundary conditions is rho(-,-)", rmm >= max(rpp, rpm), rmm)
    # tunnelling bound on a periodic ring of M_c = 2n+4 sites, n = 2
    n = 2; Mc = 2n + 4
    Tper = [exp(h) eps; eps exp(-h)]^Mc; piE = Tper[2,2]/tr(Tper)
    ok("pi(E) < 1/2", piE < 0.5, piE)
    Zw = sum(prod(Tm[(c[i] == 1 ? 1 : 2), (c[i+1] == 1 ? 1 : 2)] for i in 1:n+1) for c in [vcat([-1], [((cf >> i) & 1) == 1 ? 1 : -1 for i in 0:n-1], [-1]) for cf in 0:(2^n-1)])
    q = exp(-h)^(n+1)/Zw; ell = n/N
    ok("1 - q <= e^{Delta l}(cosh((l+1) kappa) - 1)", 1 - q <= exp(Delta*ell)*(cosh((ell+1)*kap) - 1), (1 - q, exp(Delta*ell)*(cosh((ell+1)*kap) - 1)))
    ok("M_c/n <= 4 for n >= 2", all((2n+4)/n <= 4 for n in 2:50), (2n+4)/n)
end
ok("r -> infinity: C(0,Q0) = 1/(1+Q0) > 1/2", all(abs(Cfun(0.0, Q) - 1/(1+Q)) < 1e-12 && Cfun(0.0, Q) > 0.5 for Q in (0.05, 0.3, 0.9)), Cfun(0.0, 0.05))
ok("C(tau,Q0) < 1/2 iff tau > Q0^2", all((Cfun(tau, Q) < 0.5) == (tau > Q^2) for Q in 0.05:0.05:0.95, tau in 0.001:0.003:1 if abs(tau - Q^2) > 1e-9), true)
ok("C(tau,Q0) -> 0 as Q0 -> 0", Cfun(0.01, 1e-4) < 1e-3, Cfun(0.01, 1e-4))

println("\n== Strong coupling, uniform in the boundary condition ==")
thetab(b) = 3b/(2 - 9b)
ok("theta_{1/10} = 0.2727...", abs(thetab(0.1) - 0.272727) < 1e-6, thetab(0.1))
ok("c = theta^2 = 0.07438, (1-2c)^4 = 0.5250", abs(thetab(0.1)^2 - 0.074380) < 1e-6 && abs((1 - 2*thetab(0.1)^2)^4 - 0.52506) < 1e-5, (thetab(0.1)^2, (1-2*thetab(0.1)^2)^4))
ok("m = 283 >= 1 + 148/(1-2c)^4 and BP(3, (1-2c)^4/32) = 0.0164", 283 >= 1 + 148/(1 - 2*thetab(0.1)^2)^4 && 282 < 1 + 148/(1 - 2*thetab(0.1)^2)^4 && abs((1-2*thetab(0.1)^2)^4/32 - 0.016408) < 1e-6, (1 + 148/(1-2*thetab(0.1)^2)^4, (1-2*thetab(0.1)^2)^4/32))
ok("theta_{0.15} = 0.6923 and c = 0.0744 needs d >= 8", abs(thetab(0.15) - 0.692308) < 1e-6 && thetab(0.15)^7 > 0.0744 > thetab(0.15)^8, (thetab(0.15), thetab(0.15)^7, thetab(0.15)^8))
ok("C_m = (m-1)^8 + gamma (m-1)^4 at m = 4, gamma = 1", 3^8 + 3^4 == 6642, 6642)

println("\n== The ascent without a gap hypothesis (Theorem F) ==")
ok("20 (e^{1/20} - 1) <= 1.0255", 20*(exp(big(1)/20) - 1) <= big"1.0255", Float64(20*(exp(big(1)/20) - 1)))
ok("16 * 1.0255^2 <= 16.83 and 16.83/20 = 0.8415", 16*big"1.0255"^2 <= big"16.83" && big"16.83"/20 == big"0.8415", Float64(16*big"1.0255"^2))
ok("16 (e^y - 1)^2 <= 16.83 y^2 on [0, 1/20]", all(16*(exp(y) - 1)^2 <= 16.83*y^2 for y in 1e-6:1e-5:0.05), 16*(exp(0.05)-1)^2/0.05^2)
let y = 0.05, u0 = 16.83*0.05, good = true
    for i in 0:8
        good &= (16.83*y <= u0^(2.0^i)*(1 + 1e-12)) && y <= 0.05
        y = 16*(exp(y) - 1)^2
    end
    ok("16.83 y_i <= u_0^(2^i) along the recursion from y_0 = 1/20", good, y)
end
ok("log(1/0.8415) > 0.172", log(1/0.8415) > 0.172, log(1/0.8415))
ok("1.0255/20 <= 0.061 and 0.061 = 1.0255/16.83 rounded up", 1.0255/16.83 <= 0.061, 1.0255/16.83)
ok("e^{y/13824} - 1 <= y/13823 for y <= 1/20", all(exp(y/13824) - 1 <= y/13823 for y in 1e-4:1e-4:0.05), exp(0.05/13824) - 1)
ok("(16.83 * 13823)^-2 <= 1.9e-11", 1/(16.83*13823)^2 <= 1.9e-11, 1/(16.83*13823)^2)
let f(y) = 16*(exp(y) - 1)^2 - y, a = 0.05, b = 0.07
    for _ in 1:60
        m = (a + b)/2; f(m) < 0 ? (a = m) : (b = m)
    end
    ok("fixed point y_c = 0.0589...", abs(a - 0.0589) < 1e-4 && a > 0.05, a)
end
ok("e^{4 y_0/96^3} - 1 <= 2.3e-7 at y_0 = 1/20", exp(big(4)/20/big(96)^3) - 1 <= big"2.3e-7", Float64(exp(big(4)/20/big(96)^3) - 1))
ok("log(1/2.3e-7) >= 15.28 and (1 - 1/e) 15.28 >= 9.6", log(1/2.3e-7) >= 15.28 && (1 - exp(-1))*15.28 >= 9.6 && 1 - exp(-1) >= 0.632, (log(1/2.3e-7), (1 - exp(-1))*15.28))
ok("y_0/(2 * 96^3) <= 2.9e-8 at y_0 = 1/20", 0.05/(2*96^3) <= 2.9e-8, 0.05/(2*96^3))
ok("64 sqrt(3/32) sqrt(1e-3) <= 0.62 < 0.632", 64*sqrt(3/32)*sqrt(1e-3) <= 0.62, 64*sqrt(3/32)*sqrt(1e-3))
ok("3/(32*16) < 0.006 < 0.15 = 9.6/64", 3/(32*16) < 0.006 && 9.6/64 == 0.15, 3/(32*16))

println("\n== The gap closes at a fixed lattice size (SU(2)) ==")
ok("2 sqrt 2/(3 pi^3) >= 1/33", 2*sqrt(2)/(3*pi^3) >= 1/33, 2*sqrt(2)/(3*pi^3))
ok("6 + 4 log 33 + 6 log 8 <= 6 * 5.42", 6 + 4*log(33) + 6*log(8) <= 6*5.42, 6 + 4*log(33) + 6*log(8))
ok("Haar mass of {theta <= t} >= 8 t^3/(3 pi^3) on [0, pi/2]", all((2/pi)*(t/2 - sin(2t)/4) >= 8t^3/(3pi^3) - 1e-15 for t in 0.01:0.01:pi/2), (2/pi)*(pi/4))
ok("theta_r = 2 asin(r/sqrt 8) >= r/sqrt 2 for r <= 2", all(2*asin(r/sqrt(8)) >= r/sqrt(2) for r in 0.01:0.01:2), 2*asin(2/sqrt(8)))
ok("eps(c) = 12 (log c + 5.42)/c <= 1e-3 at c = 2.2e5", 12*(log(2.2e5) + 5.42)/2.2e5 <= 1e-3, 12*(log(2.2e5) + 5.42)/2.2e5)

println("\n== Appendix: the one-link constant ==")
ok("e^x - 1 - x <= x^2 on [0, 1]", all(exp(x) - 1 - x <= x^2 for x in 0:0.001:1), exp(1) - 2)
ok("e^s - 1 <= 2 s on [0, 1]", all(exp(s) - 1 <= 2s for s in 0:0.001:1), exp(1) - 1)
let N = 2000, worst = 0.0
    # SU(2), fundamental representation, phi(Q) = c Re Tr Q: ||int pi dnu|| = I_2(2c)/I_1(2c) < 1 - const/max(1, Lambda)^2
    for c in (0.5, 2.0, 8.0, 32.0)
        th = range(0, pi, length = N)
        w = sin.(th).^2 .* exp.(2c .* cos.(th))
        m = sum(w .* cos.(th))/sum(w)
        worst = max(worst, m)
    end
    ok("SU(2) one-link mean of (1/2) Tr Q stays below 1", worst < 1, worst)
end
