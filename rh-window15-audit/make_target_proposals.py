#!/usr/bin/env python3
"""Freeze one target-adapted residual proposal per parity at u=1e-95+2^-500.

This is a proposal generator only.  Acceptance remains entirely in the
subsequent Arb higher-moment and full-inertia certificates.
"""
from pathlib import Path
from fractions import Fraction
import argparse, json, sys

pa=argparse.ArgumentParser()
pa.add_argument("--baseline",type=Path,required=True)
pa.add_argument("--out",type=Path,required=True)
pa.add_argument("--precision",type=int,default=1024)
args=pa.parse_args()
basepath=args.baseline.resolve(); out=args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(basepath))
import certify_spectrum as cs
from flint import arb,arb_mat,ctx
ctx.prec=args.precision; ctx.threads=4
cs.OUT=basepath

s=Fraction("1e-95")
for p,name in [(0,"even"),(1,"odd")]:
    A=cs.load_source("A",p); B=cs.load_source("B",p); K=cs.load_source("K",p)
    n=A.nrows(); m=arb(1)/2
    gamma=[arb(3)/(4*i+2*p+1) for i in range(n)]
    Gamma=arb_mat(n,n)
    for i in range(n): Gamma[i,i]=gamma[i]

    D=A-m*Gamma
    H=B-m*A
    Y=H.solve(D)
    J0=cs.symmetric(Gamma-D*Y)
    J1=cs.symmetric(-Y.transpose()*D*Y)
    for i in range(n):
        for j in range(n):
            J0[i,j]/=gamma[i]*gamma[j]
            J1[i,j]/=gamma[i]*gamma[j]

    u=cs.rat(s)+arb(2)**-500
    J=J0+u*J1
    T=cs.symmetric(J*K*J+(m-u)*J)
    tag=f"targetproposal_p{p}_s1em95"
    result,z=cs.full_inertia(T,tag,want_witness=True)
    expected=2 if p==0 else 1
    if result["inertia"]["negative"] != expected:
        raise SystemExit(f"unexpected baseline target inertia {name}: {result['inertia']}")

    # Shifted analogue of the zero-shift proposal a=Y*J_0*z/m.
    # Here W_u(X_0)=Gamma^(1/2) J_u Gamma^(1/2)/(m-u).
    a=Y*J*z/(m-u)
    bits=700
    coeff=[int((a[i,0].mid()*(1<<bits)+arb(1)/2).floor().unique_fmpz()) for i in range(n)]
    if not any(coeff):
        raise SystemExit(f"zero target proposal in {name}")

    rec={
        "bits":bits,
        "raw_Legendre_coefficients":[str(x) for x in coeff],
        "status":"TARGET_ADAPTED_PROPOSAL_ONLY; acceptance requires fresh higher moments and full Arb inertia",
        "parity":name,
        "proposal_threshold":"1e-95",
        "proposal_shift":"1e-95 + 2^-500",
        "baseline_target_inertia":result["inertia"],
        "target_negative_witness_energy":result["negative_witness_energy"],
        "formula":"a = Y J_u z / (m-u)"
    }
    path=out/f"p{p}_s0_residual_polynomial_proposal.json"
    path.write_text(json.dumps(rec,indent=2)+"\n")
    print("TARGET_PROPOSAL",name,result["inertia"],result["negative_witness_energy"],path)
