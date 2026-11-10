#!/usr/bin/env python3
"""Certified Window15 source-shift lower comparison at c2=4.8059.

Let H_c=(c-u)I-P and decompose relative to the finite polynomial space E:
H_c = [[A_c, B*],[B,D_c]].  The exact shorted finite source is
S_c=A_c-B*D_c^{-1}B.  Since K_c+A_c is independent of c,

T_exact(c2)-T_exact(c1)=B*(D_c1^{-1}-D_c2^{-1})B
=delta B*D_c1^{-1}D_c2^{-1}B.

With ||P||<=Bnorm, D_cj <= (cj-u+Bnorm)I, hence
T_exact(c2)-T_exact(c1) >= alpha R,
alpha=delta/[(c1-u+Bnorm)(c2-u+Bnorm)],
R=B_source-A_source Gamma^{-1} A_source = B*B.

The existing c1 inverse-majorant lower target is congruent via F to
F K F +(m-u)F.  Therefore a certified lower target for c2 is
F K F +(m-u)F + alpha F R F.
"""
from pathlib import Path
from fractions import Fraction
import argparse,json,time
from flint import arb,arb_mat,acb,ctx
import certify_spectrum as base

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--parity',choices=['even','odd'],required=True)
    pa.add_argument('--precision',type=int,default=1024)
    args=pa.parse_args()
    ctx.prec=args.precision;ctx.threads=4
    p=int(args.parity=='odd');start=time.monotonic()
    A=base.load_source('A',p);B=base.load_source('B',p);K=base.load_source('K',p)
    n=A.nrows();m=arb(1)/2
    gamma=[arb(3)/(4*i+2*p+1) for i in range(n)]
    Gamma=arb_mat(n,n);GammaInv=arb_mat(n,n)
    for i,g in enumerate(gamma):Gamma[i,i]=g;GammaInv[i,i]=1/g

    s=Fraction('5.26e-97');u=base.rat(s)+arb(2)**-500
    D=A-m*Gamma;H=B-m*A;Y=H.solve(D)
    U=base.symmetric(Gamma-D*Y);V=base.symmetric(-Y.transpose()*D*Y)
    F=base.symmetric(U+u*V)

    R=base.symmetric(B-A*GammaInv*A)
    c1=arb(24)/5;c2=arb(48059)/10000;Bnorm=arb(43)/10
    delta=c2-c1
    alpha=delta/((c1-u+Bnorm)*(c2-u+Bnorm))

    aOmega=acb(arb(1)/4,arb(768)/2).digamma().real-arb.pi().log()
    if not aOmega>c2: raise ArithmeticError('c2 exceeds high-frequency archimedean floor')
    if not c2-Bnorm>m: raise ArithmeticError('source floor did not improve')

    T0=base.symmetric(F*K*F+(m-u)*F)
    T=base.symmetric(T0+alpha*F*R*F)
    tag=f'sourceshift_p{p}_s5p26em97'
    result,_=base.full_inertia(T,tag,want_witness=False)
    receipt={
      'status':'CERTIFIED_SOURCE_SHIFT_LOWER_COMPARISON',
      'L':'3/2','parity':args.parity,'dimension':n,
      'physical_threshold':'5.26e-97','projection_debit':'2^-500',
      'c1':'24/5','c2':'48059/10000','Bnorm':'43/10',
      'new_source_floor':'5059/10000',
      'a768_minus_c2':(aOmega-c2).str(55),
      'alpha_lower':alpha.str(70),
      'shorting_identity':'T_exact(c2)-T_exact(c1) = delta B*D1^-1 D2^-1 B >= alpha R',
      'comparison_inertia':result['inertia'],
      'result':result,
      'elapsed_seconds':round(time.monotonic()-start,3)
    }
    Path(f'{args.parity}_source_shift_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)

if __name__=='__main__':main()
