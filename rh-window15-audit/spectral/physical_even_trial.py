#!/usr/bin/env python3
"""Physical even sine-trial certificate. Candidate generation is not acceptance."""
from pathlib import Path
from fractions import Fraction as F
import argparse,ast,json,hashlib,time
from flint import arb,arb_mat,ctx

OUT=Path(__file__).resolve().parent
def rat(q):
    q=F(q);return arb(q.numerator)/q.denominator
def prime_powers(L):
    """Exhaustive sieve, followed by strict Arb inclusion/exclusion tests."""
    ell=rat(L)
    assert F(0)<L<=F(3), 'this executable is scoped to 0<L<=3'
    stop=int((2*ell).exp().upper().ceil().unique_fmpz())
    prime=[True]*(stop+1);prime[0]=prime[1]=False
    for p in range(2,stop+1):
        if prime[p]:
            for m in range(p*p,stop+1,p):prime[m]=False
    result=[]
    for p in range(2,stop+1):
        if not prime[p]:continue
        nn=p
        while nn<=stop:
            if arb(nn).log()<2*ell:result.append((nn,p))
            else:assert arb(nn).log()>2*ell, 'ambiguous cutoff; increase precision'
            nn*=p
    return sorted(result)
def gauss(n):
    return [arb.legendre_p_root(n,k,weight=True) for k in range(n)]
def setup(L,n):
    ell=rat(L); a=[(2*i+1)*arb.pi()/2 for i in range(n)]
    powers=prime_powers(L)
    return ell,a,powers
def overlap(t,a,i,j):
    if i==j:return ((2-t)*(a[i]*t).cos()+(a[i]*t).sin()/a[i])/2
    return (a[i]*(a[j]*t).sin()-a[j]*(a[i]*t).sin())/(a[i]**2-a[j]**2)
def matrix_candidate(L,n,nq):
    ell,a,powers=setup(L,n);J=[arb(0) for _ in a];D=[arb(0) for _ in a]
    for x,w in gauss(nq):
        t=x+1;u=ell*t;den=u.sinh();em=(-u).exp();ep=(u/2).exp()
        for i,aa in enumerate(a):
            sn,cs=(aa*t).sin_cos()
            J[i]+=w*ell*ep*sn/den
            hii=((2-t)*cs+sn/aa)/2
            D[i]+=w*ell*(em-ep*hii)/den
    v=[2*aa*ell.sqrt()*(ell/2).cosh()/(aa*aa+ell*ell/4) for aa in a]
    c=-arb.const_euler()-arb.pi().log()-(1-(-4*ell).exp()).log()
    Q=arb_mat(n,n)
    for i in range(n):
        for j in range(i,n):
            val=2*v[i]*v[j]
            if i==j:val+=c+D[i]
            else:val-=(a[i]*J[j]-a[j]*J[i])/(a[i]**2-a[j]**2)
            for nn,p in powers:
                val-=2*arb(p).log()/arb(nn).sqrt()*overlap(arb(nn).log()/ell,a,i,j)
            Q[i,j]=Q[j,i]=val
    return Q
def propose(L,n,nq,bits):
    Q=matrix_candidate(L,n,nq)
    print("candidate matrix assembled",n,flush=True)
    v=arb_mat([[arb(1) if i==0 else arb(0)] for i in range(n)])
    for k in range(4):
        v=Q.solve(v);v=v/v[0,0]
        print("inverse iteration",k+1,flush=True)
    ints=[int(((v[i,0].mid()*(1<<bits))+arb(1)/2).floor().unique_fmpz()) for i in range(n)]
    return ints
def certify(L,ints,bits,nq):
    ell,a,powers=setup(L,len(ints));C=[F(x,1<<bits) for x in ints]
    norm=sum((c*c for c in C),F(0));N=rat(norm);c=[rat(x) for x in C]
    dd=[c[i]**2/(2*a[i]) for i in range(len(c))]
    for i in range(len(c)):
        for j in range(i+1,len(c)):
            fac=2*c[i]*c[j]/(a[i]**2-a[j]**2)
            dd[i]-=fac*a[j];dd[j]+=fac*a[i]
    def h(t):
        s=arb(0)
        for i,aa in enumerate(a):
            sn,cs=(aa*t).sin_cos()
            s+=(2-t)*c[i]**2*cs/2+dd[i]*sn
        return s
    pole_linear=sum((c[i]*2*aa*ell.sqrt()*(ell/2).cosh()/(aa*aa+ell*ell/4) for i,aa in enumerate(a)),arb(0))
    pole=2*pole_linear**2
    const=(-arb.const_euler()-arb.pi().log()-(1-(-4*ell).exp()).log())*N
    prime=sum((-2*arb(p).log()/arb(nn).sqrt()*h(arb(nn).log()/ell) for nn,p in powers),arb(0))
    quad=arb(0)
    for x,w in gauss(nq):
        t=x+1;u=ell*t
        quad+=w*ell*(N*(-u).exp()-(u/2).exp()*h(t))/u.sinh()
    bb=rat(F(7,24));T=rat(F(7,3));H=arb(0);Ht=arb(0)
    for i,aa in enumerate(a):
        ei=(aa*bb).exp()
        H+=abs(c[i])**2*(T+1/aa)*ei/2
        Ht+=abs(c[i])**2*T*aa*ei/2
        for j in range(i+1,len(c)):
            aj=a[j];ej=(aj*bb).exp();den=aj**2-aa**2
            H+=2*abs(c[i]*c[j])*(aa*ej+aj*ei)/den
            Ht+=2*abs(c[i]*c[j])*aa*aj*(ei+ej)/den
    xmin=-ell/24;xmax=ell*rat(F(49,24));Y=ell*bb
    assert Y<arb.pi()
    deriv=N*(-xmin).exp()+(xmax/2).exp()*(Ht/ell+H/2)
    mellipse=ell*deriv/(Y.sin()/Y)
    k=0
    while not mellipse<arb(2)**k:k+=1
    rem=F(32*(1<<k))*F(3,4)**(2*nq)
    q=(pole+const+prime+quad)/N
    total=q+arb(0,rat(rem/norm))
    import math
    ef=math.ceil(math.log2(abs(float(total.upper()))))
    U=F(2)**ef
    assert total<rat(U)
    receipt={
      "status":"PASS_PHYSICAL_EVEN_TRIAL_UPPER_BOUND",
      "L":str(L),"basis":"sin((2j+1)*pi*(u+L)/(2L))/sqrt(L)",
      "modes":len(ints),"coefficient_bits":bits,
      "coefficients_sha256":hashlib.sha256(json.dumps(ints,separators=(',',':')).encode()).hexdigest(),
      "norm_squared":str(norm),"active_prime_powers":[nn for nn,p in powers],
      "arithmetic_precision_bits":ctx.prec,"quadrature_nodes":nq,
      "ellipse_rho":"4/3","ellipse_bound_power_of_two":k,
      "quadrature_error_before_norm":str(rem),
      "rayleigh_interval":total.str(70),
      "quadrature_midpoint_value":q.str(70),
      "U_even":str(U),"U_even_power_of_two":ef,
      "U_even_approx":float(U),
      "spectral_interpretation":"lambda_1^+(L) <= Q(f)/||f||^2 < U_even; no lower bound on lambda_1^+",
      "physical_operator":"geometric Weil explicit formula in inherited [-L,L] normalization",
      "scope":"physical Rayleigh upper bound only; no ground-state lower bound or spectral gap at this window",
      "contribution_diagnostics":{"pole":pole.str(25),"constant":const.str(25),"prime":prime.str(25),"archimedean":quad.str(25)}
    }
    return receipt
def main():
    pa=argparse.ArgumentParser();pa.add_argument("--L",type=F,default=F(3,2));pa.add_argument("--modes",type=int,default=128);pa.add_argument("--nq",type=int,default=512)
    pa.add_argument("--precision",type=int,default=1024);pa.add_argument("--replay",type=Path);pa.add_argument("--regression",type=Path)
    args=pa.parse_args();ctx.prec=args.precision;t=time.time()
    if args.regression:
        tree=ast.parse(args.regression.read_text());d={}
        for item in tree.body:
            if isinstance(item,ast.Assign) and len(item.targets)==1 and isinstance(item.targets[0],ast.Name) and item.targets[0].id in ["COEFF","BITS"]:
                d[item.targets[0].id]=ast.literal_eval(item.value)
        rec=certify(F(59,64),d["COEFF"],d["BITS"],args.nq)
        assert arb(rec["rayleigh_interval"])<rat(F(7,2**80))
        rec["regression_historical_upper_bound"]="7/2^80"
        rec["regression_pass"]=True
        print("REGRESSION PASS",rec["rayleigh_interval"],flush=True)
        (OUT/"regression_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
        return
    if args.replay:
        trial=json.loads(args.replay.read_text());L=F(trial["L"]);ints=trial["coefficients"];bits=trial["bits"]
    else:
        L=args.L;bits=512;ints=propose(L,args.modes,args.nq,bits)
        trial={"L":str(L),"bits":bits,"coefficients":ints}
        (OUT/"even_trial.json").write_text(json.dumps(trial,indent=2)+"\n")
    rec=certify(L,ints,bits,args.nq);rec["elapsed_seconds"]=round(time.time()-t,3)
    (OUT/"even_trial_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
    print(json.dumps({k:rec[k] for k in ["status","modes","rayleigh_interval","U_even","U_even_approx","elapsed_seconds"]},indent=2),flush=True)
if __name__=="__main__":main()
