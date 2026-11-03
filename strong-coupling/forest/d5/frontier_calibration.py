"""Small exact calibrations for A/B/C and finite-dimensional D8/E examples.

Analytic statements and their hypotheses are in FRONTIER_LEMMAS.md.
No numerical output is a Yang--Mills continuum theorem.
"""
from fractions import Fraction as F
from math import factorial, exp
import json
from pathlib import Path


def bessel_interval(order, kappa, terms=18):
    z=kappa/2
    total=sum((z**(2*k+order)/F(factorial(k)*factorial(k+order)) for k in range(terms)),F(0))
    nxt=z**(2*terms+order)/F(factorial(terms)*factorial(terms+order))
    ratio=z*z/F((terms+1)*(terms+order+1))
    return total,total+nxt/(1-ratio)


def oscillator_diag(n, power):
    states={n:F(1)}
    # Unnormalised number basis: creation coefficient 1, annihilation coefficient n.
    for _ in range(power):
        out={}
        for j,c in states.items():
            out[j+1]=out.get(j+1,F(0))+c
            if j:out[j-1]=out.get(j-1,F(0))+j*c
        states=out
    assert power%2==0
    return states.get(n,F(0))/2**(power//2)

result={'classification':'exact finite checks and numerical toy examples; no continuum or full block bound'}
responses=[]
for beta in [F(351,1000),F(355,1000)]:
    l1,u1=bessel_interval(1,beta);l2,u2=bessel_interval(2,beta)
    lo,hi=(1+l2/u1)/4,(1+u2/l1)/4
    assert hi-lo<F(1,10**40)
    linear=(1+beta/4)/4
    assert hi<linear
    responses.append({'beta_W':str(beta),'attainable_response_interval':[float(lo),float(hi)],
                      'first_order_formula':float(linear),'response_needed_for_contraction':float(1/(10*beta))})
result['A_zero_exterior']=responses
assert oscillator_diag(0,4)==F(3,4)
assert oscillator_diag(2,4)==F(39,4)
result['B_oscillator']={'vacuum_first_order':str(oscillator_diag(0,4)),
                        'two_particle_first_order':str(oscillator_diag(2,4)),
                        'gap_first_order':str(oscillator_diag(2,4)-oscillator_diag(0,4))}
# Abstract overlap counterexample: all H labels share one Bernoulli witness.
p=F(1,10);H=6
assert p>p**H
result['C_witness_limit']={'labels':H,'probability_all':str(p),'invalid_independent_bound':str(p**H)}
# Dense rank-one perturbation of identity: form control, but no uniform exponential locality.
epsilon=F(1,2)
result['D8_dense_inverse']=[{'n':n,'endpoint_inverse':str(-epsilon/((1+epsilon)*n))} for n in [16,64,256]]
# Vanishing observable normalisation with a fixed positive physical gap.
result['E_trivial_normalization']=[{'a':a,'mass':1,'correlation_at_t1':a*a*exp(-1)} for a in [1,.1,.01,.001]]
# Jensen spectral bound on a finite positive measure, illustrating the proof.
weights=[F(1,4),F(3,4)];energies=[F(1),F(5)];v=sum(weights);M=sum(w*e for w,e in zip(weights,energies))
for t in [.1,1,3]:
    actual=sum(float(w)*exp(-t*float(e)) for w,e in zip(weights,energies))
    lower=float(v)*exp(-t*float(M/v))
    assert actual>=lower-1e-15
result['E_spectral_budget']={'variance':str(v),'first_moment':str(M),'jensen_checks':'passed'}
Path('result.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,indent=2,sort_keys=True))
