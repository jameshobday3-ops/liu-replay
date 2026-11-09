#!/usr/bin/env python3
"""Certified diagnostic of the perfect-inverse ceiling at L=3/2.

For H_u=M-uI and polynomial embedding E, operator convexity gives
(E* H_u^{-1} E)^{-1} <= E* H_u E.
Thus the exact shorted target C + (E*H_u^{-1}E)^{-1} is bounded above
by the direct polynomial compression C + E*H_u E.

If this ceiling has too many negative directions, no inverse-majorant /
Krylov refinement at the same polynomial compression can close the target.
"""
from pathlib import Path
from fractions import Fraction
import argparse, json, time
from flint import arb, arb_mat, ctx
import certify_spectrum as cs

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("--parity",choices=["even","odd"],required=True)
    pa.add_argument("--thresholds",default="0,5.26e-97")
    pa.add_argument("--precision",type=int,default=1024)
    args=pa.parse_args()
    ctx.prec=args.precision
    ctx.threads=4
    p=0 if args.parity=="even" else 1
    start=time.monotonic()
    A=cs.load_source("A",p)
    K=cs.load_source("K",p)
    n=A.nrows()
    gamma=[arb(3)/(4*i+2*p+1) for i in range(n)]
    Gamma=arb_mat(n,n)
    for i in range(n):
        Gamma[i,i]=gamma[i]
    results=[]
    for s_text in args.thresholds.split(","):
        s=Fraction(s_text)
        # Match the conservative finite-to-full comparison convention.
        u=cs.rat(s) + (arb(0) if s==0 else arb(2)**-500)
        T=cs.symmetric(K + A - u*Gamma)
        tag=f"ceiling_p{p}_s{s_text.replace('-','m').replace('.','p')}"
        try:
            result,_=cs.full_inertia(T,tag,want_witness=False)
            status="CERTIFIED_FULL_INERTIA"
        except Exception as exc:
            # A failure of the positive-tail factorization is itself informative,
            # but is not converted into a sign claim.
            result={"error":repr(exc)}
            status="INERTIA_ROUTINE_DID_NOT_CLOSE"
        result.update({
            "s_label":s_text,
            "physical_threshold":str(s),
            "comparison_shift": u.str(50),
            "operator":"K + A - u*Gamma",
            "interpretation":"upper ceiling for exact shorted target at fixed 864-mode polynomial compression"
        })
        results.append(result)
        print(json.dumps(result,indent=2),flush=True)
    receipt={
        "status":"WINDOW15_PERFECT_INVERSE_CEILING_DIAGNOSTIC",
        "L":"3/2",
        "parity":args.parity,
        "precision_bits":ctx.prec,
        "results":results,
        "elapsed_seconds":round(time.monotonic()-start,3)
    }
    Path(f"{args.parity}_perfect_inverse_ceiling.json").write_text(json.dumps(receipt,indent=2)+"\n")
if __name__=="__main__":
    main()
