from fold_intervals import *
import numpy as np
# floating-point reference evaluation (from the exploration script)
# Floating-point exploration of the branch structure of the cut-off average
# of one pair with its seven non-path links frozen (Prop. prop:branches data),
# all other links at I.  Quaternion conventions: (a,b,c,d) = a+bi+cj+dk,
# torus = {cos q + i sin q} <-> diag(e^{iq}, e^{-iq}); chordal HS distance
# delta(v,I) = 2 sqrt2 |sin(psi/2)|, psi = eigenangle; rho_G = 1.
import numpy as np
from numpy import cos, sin, sqrt, pi

def qmul(p, q):
    a1,b1,c1,d1 = p; a2,b2,c2,d2 = q
    return np.array([a1*a2-b1*b2-c1*c2-d1*d2,
                     a1*b2+b1*a2+c1*d2-d1*c2,
                     a1*c2-b1*d2+c1*a2+d1*b2,
                     a1*d2+b1*c2-c1*b2+d1*a2])
def qinv(p): return np.array([p[0],-p[1],-p[2],-p[3]])
def qlog(v):            # principal log, pure imaginary part, |.|=psi in [0,pi]
    n = np.linalg.norm(v[1:]); psi = np.arctan2(n, v[0])
    return (psi/n)*v[1:] if n>1e-300 else np.zeros(3)
def qexp(w):
    n = np.linalg.norm(w)
    return np.array([cos(n), *((sin(n)/n)*w)]) if n>1e-300 else np.array([1.,0,0,0])
def chordal(v): return 2*sqrt(2)*sqrt(max(0.,1-v[0]**2)) if False else sqrt(2)*np.linalg.norm(v-np.array([1,0,0,0]))
def phi(m): return min(1., max(0., 2-2*m))
def torus(q): return np.array([cos(q), sin(q), 0, 0])

theta = np.array([-0.5] + [0.49]*6)          # seven non-path links
L = [torus(t) for t in theta]

def Vsharp(g, L=L):
    vs = [qmul(qinv(g), l) for l in L]
    m = max(chordal(v) for v in vs)
    chi = phi(m)
    N = sum(qlog(v) for v in vs)/8.0
    return qmul(g, qexp(chi*N)), chi, m

def f0(q):               # scalar angle of V on the torus
    V, chi, m = Vsharp(torus(q)); return np.arctan2(V[1], V[0])


def gS(q, s): return np.array([cos(s)*cos(q), cos(s)*sin(q), 0., sin(s)])
# compare interval V with floating V at a few points, and derivative with finite differences
for (q, s) in [(0.0, 0.0), (0.05, 0.3), (-0.1, 0.1), (0.2, 0.6), (-0.2, 0.02)]:
    V, chi, m = V_of(I(q), D(s, 1))
    Vf, chif, mf = Vsharp(gS(q, s)); Vf2, _, _ = Vsharp(gS(q, s + 1e-7))
    for i in range(4):
        assert lo(V[i].v) - 1e-12 <= Vf[i] <= hi(V[i].v) + 1e-12, (q, s, i, V[i].v, Vf[i])
        fd = (Vf2[i] - Vf[i]) / 1e-7
        assert lo(V[i].d) - 1e-5 <= fd <= hi(V[i].d) + 1e-5, (q, s, i, V[i].d, fd)
    assert lo(chi.v) - 1e-12 <= chif <= hi(chi.v) + 1e-12 and lo(m.v) - 1e-12 <= mf <= hi(m.v) + 1e-12, (chi.v, chif, m.v, mf)
    print("ok", q, s, [float(hi(V[i].v) - lo(V[i].v)) for i in range(4)], float(hi(V[3].d) - lo(V[3].d)))
# box containment test: random points in a box
import random
q0, q1, s0, s1 = -0.01, 0.0, 0.0, 0.02
V, chi, m = V_of(I(q0, q1), D(I(s0, s1), 1))
for _ in range(200):
    q = random.uniform(q0, q1); s = random.uniform(s0, s1)
    Vf, _, _ = Vsharp(gS(q, s))
    for i in range(4): assert lo(V[i].v) - 1e-12 <= Vf[i] <= hi(V[i].v) + 1e-12
print("box ok; V_k box:", V[3].v, " dV_k/ds box:", V[3].d)
