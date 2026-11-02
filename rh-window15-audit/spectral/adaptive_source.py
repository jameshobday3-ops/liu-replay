#!/usr/bin/env python3
"""Propose positive weights, then certify a continuum shift-operator bound.

NumPy is used only to propose weights. Acceptance consists of integer
inequalities and Arb enclosures of log, sqrt and digamma. No finite-grid
eigenvalue is substituted for an infinite-dimensional operator bound.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math, time
import numpy as np
from flint import arb, acb, ctx
from physical_even_trial import rat, prime_powers

OUT=Path(__file__).resolve().parent

def shifted_pair_max(v,k):
    n=len(v);r=np.zeros(n)
    if k<n:
        plus=np.zeros(n);plus[:n-k]=v[k:]
        if k+1<n:plus[:n-k-1]=np.maximum(v[k:-1],v[k+1:])
        minus=np.zeros(n);minus[k:]=v[:n-k]
        if k+1<n:minus[k+1:]=np.maximum(v[1:n-k],v[:n-k-1])
        r=plus+minus
    return r

def input_enclosures(L,n,coefficient_bits=80):
    ell=rat(L);h=2*ell/n;items=[]
    for nn,p in prime_powers(L):
        d=arb(nn).log();a=arb(p).log()/arb(nn).sqrt()
        k=int((d/h).floor().unique_fmpz())
        assert k<d/h<k+1
        ai=int((a*(1<<coefficient_bits)).upper().ceil().unique_fmpz())
        assert a<rat(F(ai,1<<coefficient_bits))
        items.append(dict(prime_power=nn,prime=p,shift_cells=k,
                          coefficient_upper_numerator=ai,
                          coefficient_bits=coefficient_bits,
                          coefficient_interval=a.str(45),
                          shift_ratio_interval=(d/h).str(45)))
    return items

def proposal(items,n,iterations,weight_bits):
    v=np.ones(n)
    coeff=[float(F(it['coefficient_upper_numerator'],1<<it['coefficient_bits'])) for it in items]
    for itn in range(iterations):
        pv=sum((a*shifted_pair_max(v,it['shift_cells']) for a,it in zip(coeff,items)),np.zeros(n))
        v=pv+v
        v=(v+v[::-1])/2
        v/=v.max()
    weights=[max(1,int(round(float(x)*(1<<weight_bits)))) for x in v]
    assert weights==weights[::-1]
    return weights

def verify(L,weights,baseline=None):
    n=len(weights);assert n>=2 and all(type(x) is int and x>0 for x in weights)
    items=input_enclosures(L,n);cb=items[0]['coefficient_bits'];scale=1<<cb
    totals=[0]*n
    def w(i):return weights[i] if 0<=i<n else 0
    for it in items:
        k=it['shift_cells'];ai=it['coefficient_upper_numerator']
        for i in range(n):
            totals[i]+=ai*(max(w(i+k),w(i+k+1))+max(w(i-k),w(i-k-1)))
    worst=max(range(n),key=lambda i:F(totals[i],scale*weights[i]))
    ratio=F(totals[worst],scale*weights[worst])
    B=F(math.ceil(100*ratio),100)
    if baseline is None:baseline=B+F(1,2)
    slacks=[baseline.numerator*scale*weights[i]-baseline.denominator*totals[i] for i in range(n)]
    assert min(slacks)>0
    floor=baseline-B
    assert floor>=F(1,2)
    assert all(totals[i]*B.denominator<=B.numerator*scale*weights[i] for i in range(n))
    ell=rat(L)
    crude=2*sum((arb(it['prime']).log()/arb(it['prime_power']).sqrt() for it in items),arb(0))
    chain=arb(0)
    for it in items:
        d=arb(it['prime_power']).log()
        r=int((2*ell/d).ceil().unique_fmpz())
        assert r-1<2*ell/d<r
        chain+=2*arb(it['prime']).log()/arb(it['prime_power']).sqrt()*(arb.pi()/(r+1)).cos()
    constant_rayleigh=2*sum((arb(it['prime']).log()/arb(it['prime_power']).sqrt()*(1-arb(it['prime_power']).log()/(2*ell)) for it in items),arb(0))
    result={
      'status':'PASS_CONTINUUM_WEIGHTED_SCHUR_SOURCE', 'L':str(L),
      'cells':n,'arithmetic_precision_bits':ctx.prec,
      'weights_sha256':hashlib.sha256(json.dumps(weights,separators=(',',':')).encode()).hexdigest(),
      'min_weight':min(weights),'max_weight':max(weights),
      'active_prime_powers':[it['prime_power'] for it in items],
      'input_enclosures':items,'worst_cell':worst,
      'exact_max_row_ratio_upper':str(ratio),
      'max_row_ratio_upper_decimal':rat(ratio).str(45),
      'B_norm_upper':str(B),'baseline_c':str(baseline),
      'source_floor_mu':str(floor),
      'source_ceiling':str(baseline+B),
      'minimum_integer_slack':str(min(slacks)),
      'constant_trial_P_expectation':constant_rayleigh.str(45),
      'crude_triangle_norm_upper':crude.str(45),
      'individual_shift_chain_norm_upper':chain.str(45),
      'interpretation':'-B I <= P_L <= B I; (c-B) I <= c I-P_L <= (c+B) I; does not assert Weil positivity',
    }
    result['kernel_budget']=kernel_budget(L,baseline)
    half_omega=result['kernel_budget']['omega']//2
    necessary_c=constant_rayleigh+rat(F(1,2))
    half_margin=necessary_c-a_of_t(half_omega)
    result['scalar_cutoff_obstruction']={
      'assumed_source_floor':'1/2',
      'necessary_c_lower_bound':necessary_c.str(45),
      'tested_cutoff':half_omega,
      'a_at_tested_cutoff':a_of_t(half_omega).str(45),
      'necessary_c_minus_a':half_margin.str(45),
      'tested_cutoff_impossible':bool(half_margin>0),
      'meaning':'If c I-P >= (1/2) I and a(t)-c >= 0 for all t >= Omega, then a(Omega) >= <1,P1>/||1||^2 + 1/2. A positive margin rules out the tested cutoff.'}
    return result

def a_of_t(t):return acb(arb(1)/4,rat(t)/2).digamma().real-arb.pi().log()

def kernel_budget(L,c,target_bits=390):
    ell=rat(L);cc=rat(c);omega=1
    while not a_of_t(omega)>cc:omega*=2
    endpoint=a_of_t(omega);a0=a_of_t(0)
    D=max(cc-a0,endpoint-cc)
    log_arch=(24*ell*omega*D/arb.pi()).log()+5*omega*ell/12
    log_pole=(48*(ell*ell.sinh()).sqrt()).log()+13*ell/24
    log_coefficient=arb(2).log()+max(log_arch,log_pole)
    need=(log_coefficient+target_bits*arb(2).log())/rat(F(3,2)).log()
    N=int(need.upper().ceil().unique_fmpz())
    log_error=log_coefficient-N*rat(F(3,2)).log()
    assert log_error < -target_bits*arb(2).log()
    return {'omega':omega,'a_omega':endpoint.str(45),
      'high_frequency_margin':(endpoint-cc).str(45),
      'D_bound':D.str(45),'target_error':'2^(-%d)'%target_bits,
      'polynomial_modes_N':N,'degree_max':N-1,
      'even_modes':(N+1)//2,'odd_modes':N//2,
      'log_error_bound':log_error.str(45),
      'bound':'[48 sqrt(L sinh L) exp(13L/24) + (24 L Omega D/pi) exp(5 Omega L/12)] (2/3)^N',
      'status':'analytic operator-norm tail budget only; finite moment matrices and inertia not computed'}

def main():
    pa=argparse.ArgumentParser();pa.add_argument('--L',type=F,default=F(3,2));pa.add_argument('--cells',type=int,default=16384)
    pa.add_argument('--iterations',type=int,default=600);pa.add_argument('--replay',type=Path)
    args=pa.parse_args();ctx.prec=256;t=time.time()
    if args.replay:
        raw=json.loads(args.replay.read_text());L=F(raw['L']);weights=raw['weights']
    else:
        L=args.L;items=input_enclosures(L,args.cells);weights=proposal(items,args.cells,args.iterations,40)
        raw={'L':str(L),'weights':weights,'weight_scale_bits':40,'proposal_iterations':args.iterations}
    tag=str(L).replace('/','_')
    (OUT/f'source_weights_L{tag}.json').write_text(json.dumps(raw,separators=(',',':'))+'\n')
    receipt=verify(L,weights);receipt['elapsed_seconds']=round(time.time()-t,3)
    (OUT/f'source_receipt_L{tag}.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:receipt[k] for k in ['status','L','B_norm_upper','baseline_c','source_floor_mu','max_row_ratio_upper_decimal','kernel_budget','elapsed_seconds']},indent=2),flush=True)

if __name__=='__main__':main()
