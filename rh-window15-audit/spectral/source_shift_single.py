#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction
import argparse,json,time
from flint import arb,arb_mat,acb,ctx
import certify_spectrum as base

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--parity',choices=['even','odd'],default='even')
    pa.add_argument('--threshold',required=True)
    pa.add_argument('--precision',type=int,default=1024)
    args=pa.parse_args()
    ctx.prec=args.precision;ctx.threads=4
    p=int(args.parity=='odd');start=time.monotonic()
    A=base.load_source('A',p);B=base.load_source('B',p);K=base.load_source('K',p)
    n=A.nrows();m=arb(1)/2
    gamma=[arb(3)/(4*i+2*p+1) for i in range(n)]
    Gamma=arb_mat(n,n);GammaInv=arb_mat(n,n)
    for i,g in enumerate(gamma):Gamma[i,i]=g;GammaInv[i,i]=1/g
    D=A-m*Gamma;H=B-m*A;Y=H.solve(D)
    U=base.symmetric(Gamma-D*Y);V=base.symmetric(-Y.transpose()*D*Y)
    R=base.symmetric(B-A*GammaInv*A)
    sp=Fraction(args.threshold);u=base.rat(sp)+arb(2)**-500
    c1=arb(24)/5;c2=arb(48059)/10000;Bnorm=arb(43)/10
    delta=c2-c1;alpha=delta/((c1-u+Bnorm)*(c2-u+Bnorm))
    aOmega=acb(arb(1)/4,arb(768)/2).digamma().real-arb.pi().log()
    if not aOmega>c2:raise ArithmeticError('high-frequency floor fail')
    F=base.symmetric(U+u*V)
    T=base.symmetric(F*K*F+(m-u)*F+alpha*F*R*F)
    tag='sourceshift_single_p%d_s%s'%(p,args.threshold.replace('-','m').replace('.','p'))
    result,_=base.full_inertia(T,tag,want_witness=False)
    rec={'status':'CERTIFIED_SOURCE_SHIFT_SINGLE_THRESHOLD','L':'3/2','parity':args.parity,
         'physical_threshold':args.threshold,'projection_debit':'2^-500','c2':'48059/10000',
         'alpha_lower':alpha.str(65),'a768_minus_c2':(aOmega-c2).str(55),
         'inertia':result['inertia'],'result':result,'elapsed_seconds':round(time.monotonic()-start,3)}
    Path('single_threshold_receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2),flush=True)
if __name__=='__main__':main()
