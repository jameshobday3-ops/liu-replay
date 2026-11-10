#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction
import argparse,json,time
from flint import arb,arb_mat,acb,ctx
import certify_spectrum as base

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--parity',choices=['even','odd'],required=True)
    pa.add_argument('--precision',type=int,default=1024)
    pa.add_argument('--thresholds',default='0,5e-97,1e-95,1e-92,1e-90')
    args=pa.parse_args()
    ctx.prec=args.precision;ctx.threads=4
    p=int(args.parity=='odd');start=time.monotonic()

    A=base.load_source('A',p);B=base.load_source('B',p);K=base.load_source('K',p)
    n=A.nrows();m=arb(1)/2
    gamma=[arb(3)/(4*i+2*p+1) for i in range(n)]
    Gamma=arb_mat(n,n);GammaInv=arb_mat(n,n)
    for i,g in enumerate(gamma):
        Gamma[i,i]=g;GammaInv[i,i]=1/g

    D=A-m*Gamma;H=B-m*A;Y=H.solve(D)
    U=base.symmetric(Gamma-D*Y)
    V=base.symmetric(-Y.transpose()*D*Y)
    R=base.symmetric(B-A*GammaInv*A)

    c1=arb(24)/5;c2=arb(48059)/10000;Bnorm=arb(43)/10
    delta=c2-c1
    aOmega=acb(arb(1)/4,arb(768)/2).digamma().real-arb.pi().log()
    if not aOmega>c2: raise ArithmeticError('c2 exceeds high-frequency floor')
    if not c2-Bnorm>m: raise ArithmeticError('new source floor not > 1/2')

    results=[]
    for st in args.thresholds.split(','):
        sp=Fraction(st)
        u=base.rat(sp)+arb(2)**-500
        F=base.symmetric(U+u*V)
        alpha=delta/((c1-u+Bnorm)*(c2-u+Bnorm))
        T=base.symmetric(F*K*F+(m-u)*F+alpha*F*R*F)
        tag='sourceshift_scan_p%d_s%s'%(p,st.replace('-','m').replace('.','p'))
        result,_=base.full_inertia(T,tag,want_witness=False)
        row={'physical_threshold':st,'alpha_lower':alpha.str(60),
             'inertia':result['inertia'],'result':result}
        results.append(row)
        print('SCAN',args.parity,st,result['inertia'],flush=True)

    receipt={
      'status':'CERTIFIED_SOURCE_SHIFT_THRESHOLD_SCAN',
      'L':'3/2','parity':args.parity,'dimension':n,
      'projection_debit':'2^-500','c1':'24/5','c2':'48059/10000','Bnorm':'43/10',
      'new_source_floor':'5059/10000','a768_minus_c2':(aOmega-c2).str(55),
      'results':results,'elapsed_seconds':round(time.monotonic()-start,3)
    }
    Path(f'{args.parity}_source_shift_scan_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)

if __name__=='__main__':
    main()
