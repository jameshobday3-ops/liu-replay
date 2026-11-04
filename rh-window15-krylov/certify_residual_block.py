"""Adaptive low-rank block certificate inside span(E, Q M E A_seed).

All seed vectors are exact dyadic polynomial coefficient vectors.  For each
seed a_j, v_j=(I-Pi) M E a_j.  The script consumes certified C3 columns
E^*M^3E a_j and a certified C4 Gram matrix a_i^*E^*M^4E a_j, forms the
optimized inverse majorant for the whole residual block, and certifies the
full 864-dimensional comparison inertia.

If requested, a new proposal is generated from the most negative certified
comparison witness and Gamma-orthogonalized against all existing seeds.
Proposal generation is not a proof step; all acceptance happens only after
fresh moment and inertia certificates.
"""
from pathlib import Path
from fractions import Fraction
import argparse,json,time,sys
import numpy as np
import scipy.linalg as sla
import mpmath as mp
from flint import arb,arb_mat,ctx


def sym(A): return (A+A.transpose())/2

def load_interval_vector(path):
    d=json.loads(Path(path).read_text()); den=1 << (d['bits']+1)
    return [arb(int(lo)+int(hi))/den + arb(0,arb(int(hi)-int(lo))/den) for lo,hi in d['rows']]

def load_seed(path):
    d=json.loads(Path(path).read_text()); bits=int(d['bits'])
    vals=[arb(int(x))/(1<<bits) for x in d['raw_Legendre_coefficients']]
    if len(vals)!=864: raise ValueError((path,len(vals)))
    return vals,d

def save_seed(vec,path,p,parents,diagnostics,bits=700):
    ints=[int((x.mid()*(1<<bits)+arb(1)/2).floor().unique_fmpz()) for x in vec]
    obj={'parity':p,'bits':bits,'raw_Legendre_coefficients':[str(x) for x in ints],
         'status':'FROZEN_PROPOSAL_ONLY; requires fresh certified C3/C4 moments and full Arb inertia',
         'source':'adaptive Gamma-orthogonalized defect proposal in first block-Krylov residual space Q M E',
         'parents':[str(x) for x in parents], 'diagnostics':diagnostics}
    Path(path).write_text(json.dumps(obj,indent=2)+'\n')
    return obj

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--baseline',type=Path,required=True)
    pa.add_argument('--parity',choices=['even','odd'],required=True)
    pa.add_argument('--seed',action='append',required=True,help='exact dyadic seed JSON; repeat')
    pa.add_argument('--c3',action='append',required=True,help='certified C3 column JSON matching seeds')
    pa.add_argument('--c4-gram',type=Path,required=True,help='certified interval r x r C4 Gram JSON')
    pa.add_argument('--threshold',default='5.26e-97')
    pa.add_argument('--precision',type=int,default=1024)
    pa.add_argument('--tag',required=True)
    pa.add_argument('--next-seed',type=Path)
    args=pa.parse_args()
    if len(args.seed)!=len(args.c3): raise SystemExit('seed/C3 count mismatch')
    start=time.monotonic();ctx.prec=args.precision;ctx.threads=4
    p=int(args.parity=='odd');rblk=len(args.seed);n=864;m=arb(1)/2
    sys.path.insert(0,str(args.baseline.resolve()))
    import certify_spectrum as base
    base.OUT=args.baseline
    A=base.load_source('A',p);B=base.load_source('B',p);K=base.load_source('K',p)
    # Freeze the source directory as read-only input; all new certificates go here.
    base.OUT=Path.cwd()
    gamma=[arb(3)/(4*i+2*p+1) for i in range(n)];Gamma=arb_mat(n,n)
    for i in range(n): Gamma[i,i]=gamma[i]
    seeds=[];seedmeta=[]
    for path in args.seed:
        v,d=load_seed(path);seeds.append(arb_mat([[x] for x in v]));seedmeta.append(d)
    C3=[arb_mat([[x] for x in load_interval_vector(path)]) for path in args.c3]
    gd=json.loads(args.c4_gram.read_text());
    if gd['dimension']!=rblk: raise SystemExit(('C4 dimension',gd['dimension'],rblk))
    den=1 << (gd['bits']+1);C4=arb_mat(rblk,rblk)
    for i,row in enumerate(gd['rows']):
        for j,(lo,hi) in enumerate(row): C4[i,j]=arb(int(lo)+int(hi))/den+arb(0,arb(int(hi)-int(lo))/den)
    C4=sym(C4)
    Aa=[A*a for a in seeds];Ba=[B*a for a in seeds]
    b=[arb_mat([[Aa[k][i,0]/gamma[i]] for i in range(n)]) for k in range(rblk)]
    R=[Ba[k]-A*b[k] for k in range(rblk)]
    Tcol=[C3[k]-B*b[k] for k in range(rblk)]
    Nu=arb_mat(rblk,rblk);Eta=arb_mat(rblk,rblk);Theta=arb_mat(rblk,rblk)
    for i in range(rblk):
        for j in range(rblk):
            Nu[i,j]=(seeds[i].transpose()*B*seeds[j])[0,0]-(Aa[i].transpose()*b[j])[0,0]
            Eta[i,j]=(seeds[i].transpose()*C3[j])[0,0]-(seeds[i].transpose()*B*b[j])[0,0]-(b[i].transpose()*B*seeds[j])[0,0]+(b[i].transpose()*A*b[j])[0,0]
            Theta[i,j]=C4[i,j]-(b[j].transpose()*C3[i])[0,0]-(b[i].transpose()*C3[j])[0,0]+(b[i].transpose()*B*b[j])[0,0]
    Nu=sym(Nu);Eta=sym(Eta);Theta=sym(Theta)
    s=Fraction(args.threshold);u=base.rat(s)+arb(2)**-500
    if not (0<u<m): raise SystemExit('bad shift')
    D=A-m*Gamma;H0=B-(m+u)*A+m*u*Gamma
    Hcols=arb_mat(n,rblk)
    for k in range(rblk):
        hk=Tcol[k]-(m+u)*R[k]
        for i in range(n):Hcols[i,k]=hk[i,0]
    G=sym(Theta-(m+u)*Eta+m*u*Nu)
    rhs=arb_mat(n,n+rblk)
    for i in range(n):
        for j in range(n):rhs[i,j]=D[i,j]
        for k in range(rblk):rhs[i,n+k]=Hcols[i,k]
    solved=H0.solve(rhs);Y=base.sub(solved,0,n,0,n);YH=base.sub(solved,0,n,n,n+rblk)
    S=sym(G-Hcols.transpose()*YH)
    srec,SV,negcols=base.finite_head(S,args.tag+'_residual_block')
    if srec['negative'] or srec['zero'] or srec['positive']!=rblk:
        raise ArithmeticError(f'residual block not positive: {srec}')
    RR=arb_mat(n,rblk)
    for k in range(rblk):
        for i in range(n):RR[i,k]=R[k][i,0]
    W=RR-D*YH
    numerator=base.symmetric(Gamma-D*Y-W*S.solve(W.transpose()))
    J=arb_mat(n,n)
    for i in range(n):
        for j in range(n):J[i,j]=numerator[i,j]/(gamma[i]*gamma[j])
    T=base.symmetric(J*K*J+(m-u)*J)
    result,_=base.full_inertia(T,args.tag,want_witness=False)
    result.update({'s_label':args.threshold,'physical_threshold':str(s),'projection_debit':'2^-500',
                   'inverse_trial_dimension':n+rblk,'constraint_dimension':n,'residual_block_dimension':rblk,
                   'residual_block_certificate':srec,'tag':args.tag})
    receipt={'status':'PASS_ADAPTIVE_RESIDUAL_BLOCK_COMPARISON','L':'3/2','parity':args.parity,
             'arithmetic_precision_bits':ctx.prec,'seed_count':rblk,'result':result,
             'elapsed_seconds':round(time.monotonic()-start,3)}
    Path(args.tag+'_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2),flush=True)
    target_neg = 1 if p==0 else 0
    if p==0 and result['inertia']['negative'] < 1:
        raise ArithmeticError('even comparison has no negative direction despite certified physical trial below threshold')
    if args.next_seed and result['inertia']['negative']>target_neg:
        # Discovery only: use the high-precision midpoint Schur eigenvector and
        # a float tail backsolve, then freeze the resulting direction dyadically.
        # No sign claim depends on these floating calculations; the next rank is
        # accepted only after fresh directed C3/C4 moments and Arb inertia.
        sd=json.loads(Path(str(args.tag)+'_schur.json').read_text());mp.mp.dps=180
        denmp=mp.mpf(2)**(sd['bits']+1);Sm=mp.matrix(sd['nrows'])
        for i,row in enumerate(sd['rows']):
            for j,(lo,hi) in enumerate(row):Sm[i,j]=(mp.mpf(lo)+mp.mpf(hi))/denmp
        vals,Q=mp.eigsy(Sm);head=np.array([float(Q[i,0]) for i in range(64)],dtype=float)
        def as_float(M): return np.array([[float(M[i,j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])
        Tf=as_float(T);Rtf=(Tf[64:,64:]+Tf[64:,64:].T)/2;crossf=Tf[:64,64:]
        lower=-sla.cho_solve(sla.cho_factor(Rtf,lower=True,check_finite=True),crossf.T@head,check_finite=False)
        z=np.concatenate([head,lower]);z/=np.linalg.norm(z)
        Af=as_float(A);Bf=as_float(B);Jf=as_float(J);gamf=np.array([float(x) for x in gamma])
        Gf=np.diag(gamf);Df=Af-0.5*Gf
        Y0=sla.solve(Bf-0.5*Af,Df,assume_a='sym',check_finite=False)
        anew=Y0@(Jf@z)/0.5
        projections=[]
        for aa in seeds:
            av=np.array([float(aa[i,0]) for i in range(n)])
            coef=float(np.sum(gamf*av*anew)/np.sum(gamf*av*av));projections.append(coef);anew-=coef*av
        norm=float(np.sqrt(np.sum(gamf*anew*anew)))
        if not (norm>0 and np.isfinite(norm)):raise ArithmeticError('proposal collapsed after Gamma orthogonalization')
        anew/=norm;bits=700;ints=[int(round(float(x)*(2**bits))) for x in anew]
        obj={'parity':p,'bits':bits,'raw_Legendre_coefficients':[str(x) for x in ints],
             'status':'FROZEN_PROPOSAL_ONLY; acceptance requires fresh directed C3/C4 moments and full Arb inertia',
             'source':'adaptive Gamma-orthogonalized lowest-Schur defect proposal in first block-Krylov residual space Q M E',
             'parents':[str(x) for x in args.seed],
             'diagnostics':{'current_inertia':result['inertia'],'lowest_schur_eigenvalue':mp.nstr(vals[0],50),
                            'next_schur_eigenvalue':mp.nstr(vals[1],50) if len(vals)>1 else None,
                            'gamma_projection_coefficients':projections,'gamma_norm_after':1.0,
                            'max_abs_coefficient':float(np.max(np.abs(anew)))}
             }
        args.next_seed.write_text(json.dumps(obj,indent=2)+'\n')
        print('NEXT_SEED',args.next_seed,obj['diagnostics'],flush=True)

if __name__=='__main__': main()
