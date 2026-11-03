"""Rational arithmetic certificate for CANONICAL_CHART.md constants.

Taylor series use explicit positive tails. This checks scalar inequalities;
it does not replace the analytic proof of the chart or its Hessian.
"""
from fractions import Fraction as F


def exp_upper(t, n=16):
    term = total = F(1)
    for k in range(1, n + 1):
        term *= t / k
        total += term
    nxt = term * t / (n + 1)
    return total + nxt / (1 - t / (n + 2))


def minus_log_upper(t, n=16):
    assert 0 <= t < 1
    total = sum((t**k / k for k in range(1, n + 1)), F(0))
    return total + t**(n + 1) / ((n + 1) * (1 - t))

rho = F(1, 10000)
R = 44 * rho
s = exp_upper(46 * rho) - 1
contraction = F(7, 8) * (1 + s) / (1 - s)
image = F(7, 8) * minus_log_upper(s)
error = 14 * s / (1 - s)
logdet = 3 * minus_log_upper(error)
assert contraction < F(9, 10)
assert image < R
assert error < F(7, 100)
assert logdet < 1
assert exp_upper(rho) - 1 < F(1, 100)
first = R / (rho / 2)
second = 2 * R / (rho / 2)**2
accel = second + 2 * first**2
assert first == 88 and second == 3520000 and accel < 4000000
kappa = F(1, 10**10)
accel_total = 2 + 2 * 4000000
lower_Q = F(98, 100) * (F(1, 26) - 54000 * kappa) - 12 * kappa * accel_total
upper_Q = (48 + 84 * kappa) * (1 + 2 * first**2) + 12 * kappa * accel_total
assert lower_Q > F(1, 52)
assert upper_Q < 800000
logdet_hess = 16 / rho**2
assert logdet_hess == 1600000000
lam = F(2 * 10**11)
assert lam / 52 - logdet_hess >= lam / 104
print('Rational upper bounds (displayed as decimals; comparisons use Fractions):')
for name, val in [('complex contraction', contraction), ('fixed-point image radius', image),
                  ('normalised derivative error', error), ('logdet residual', logdet),
                  ('path acceleration bound', accel), ('Q Hessian lower bound', lower_Q),
                  ('Q Hessian absolute upper bound', upper_Q)]:
    print(f'{name}: {float(val):.12g}')
print('Sum logdet Hessian absolute bound: 1600000000')
print('lambda >= 200000000000 implies Hess V_eff >= lambda/104')
print('All rational inequalities passed.')

alpha = F(1,10**14)
neighbours = 45 * 7**4
weighted_error = neighbours * 10**6 * (exp_upper(3*alpha)-1)
assert weighted_error < F(1,208)
assert upper_Q*lam+2400000000 < 10**6*lam
assert 208*(1+2*88**2)<3300000
print(f'Weighted Hessian error / lambda <= {float(weighted_error):.12g} < 1/208')
print('Restricted-chart covariance amplitude <= 3300000/lambda, rate 1e-14')
