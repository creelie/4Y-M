# Appendix B: Barta computations on S^3 and the Riccati barriers for A = I2/I1.
import sympy as sp
th,k,a=sp.symbols('theta kappa a',real=True)
b=2*sp.cot(th)-k*sp.sin(th)
psi=sp.sin(th)*sp.exp(a*k*sp.cos(th))
L1=lambda f: sp.diff(f,th,2)+b*sp.diff(f,th)-2*f/sp.sin(th)**2
r1=sp.simplify(sp.expand(-L1(psi)/psi))
print('-L1 psi/psi =',sp.factor(sp.simplify(r1)))
print('a=-1/5:',sp.simplify(r1.subs(a,-sp.Rational(1,5))-3))
r0=sp.simplify(r1+k*sp.cos(th))
print('a=-2/5 (zonal, with +kappa cos):',sp.simplify(r0.subs(a,-sp.Rational(2,5))-3))
# Riccati: A' = 1 - A^2 - 3A/k
F=lambda A: 1-A**2-3*A/k
s=sp.sqrt(k**2+4); u=sp.sqrt(k**2+9)
lo=k/(2+s); hi=k/(1+u)
print('sub residual:',sp.simplify(sp.diff(lo,k)-F(lo)), ' expect',sp.simplify(-k**2/(s*(s+2)**2)))
print('super residual:',sp.simplify(sp.diff(hi,k)-F(hi)-(u-3)**2/(u*(1+u)**2)))
print('series lo, hi:',sp.series(lo,k,0,4),sp.series(hi,k,0,4))
print('series A:',sp.series(sp.besseli(2,k)/sp.besseli(1,k),k,0,6))
