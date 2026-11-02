"""Whole-space source floor and compact-kernel error bounds."""
import json
from pathlib import Path
from fractions import Fraction
from flint import arb,acb,ctx
import adaptive_source
OUT=Path(__file__).resolve().parent
def run():
    ctx.prec=512
    w=json.loads((OUT/'inputs/source_weights_L3_2.json').read_text())
    source=adaptive_source.verify(Fraction(w['L']),w['weights'])
    assert source['baseline_c']=='24/5' and source['source_floor_mu']=='1/2' and source['B_norm_upper']=='43/10'
    L=arb(3)/2;Omega=arb(768);c=arb(24)/5;N=1728
    a0=acb(arb(1)/4).digamma().real-arb.pi().log()
    aO=acb(arb(1)/4,Omega/2).digamma().real-arb.pi().log()
    assert aO>c
    D=max(c-a0,aO-c);rho=arb(5)/2;Ar=2/(1-1/rho);ar=(rho+1/rho)/2;br=(rho-1/rho)/2
    tail=(8*Ar*(L*L.sinh()).sqrt()*(ar*L/2).exp()+4*Ar*L*Omega*D/arb.pi()*(br*L*Omega).exp())*rho**-N
    assert tail<arb(2)**-500
    degree=6720;x=2*L*Omega
    trunc=(2*L.cosh()*L**degree+Omega*D/arb.pi()*x**degree)/arb.fac_ui(degree)
    midpoint=arb(2)**-4097*x.cosh()
    rawerror=(2*L)**2*(trunc+midpoint)
    assert rawerror<arb(2)**-650
    result={'status':'PASS_SOURCE_AND_KERNEL_BUDGET','L':'3/2','source_floor':'1/2','source_upper':'91/10','source_receipt':source,'Omega':768,'polynomial_total_modes':1728,'modes_per_parity':864,'aOmega_minus_c':(aO-c).str(45),'projection_error':tail.str(45),'projection_debit':'2^-500','kernel_taylor_error':trunc.str(45),'moment_midpoint_error':midpoint.str(45),'raw_K_approximation_error':rawerror.str(45),'raw_K_error_upper':'2^-650'}
    (OUT/'budget_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='source_receipt'},indent=2),flush=True)
    return result
if __name__=='__main__':run()
