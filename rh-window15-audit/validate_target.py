import json, pathlib
root=pathlib.Path(__file__).resolve().parent
even=json.loads((root/'enrich/even_enriched_receipt.json').read_text())
odd=json.loads((root/'enrich/odd_enriched_receipt.json').read_text())
trial=json.loads((root/'spectral/even_trial_receipt.json').read_text())
target='5.26e-97'
def pick(rec):
    xs=[x for x in rec['results'] if x.get('s_label')==target]
    if len(xs)!=1: raise SystemExit(f'missing unique {target} result')
    return xs[0]
e=pick(even); o=pick(odd)
print('EVEN',e['inertia'])
print('ODD',o['inertia'])
print('TRIAL',trial['rayleigh_interval'],trial['U_even'])
if e['inertia']!={'positive':863,'negative':1,'zero':0}: raise SystemExit('TARGET FAIL: even inertia is not (863,1,0)')
if o['inertia']!={'positive':864,'negative':0,'zero':0}: raise SystemExit('TARGET FAIL: odd inertia is not (864,0,0)')
from decimal import Decimal
if not Decimal(str(trial['U_even_approx'])) < Decimal(target): raise SystemExit('TARGET FAIL: trial upper bound not below threshold')
print('TARGET PASS: simple isolated even ground state at L=1.5 certified relative to s*=5.26e-97')
