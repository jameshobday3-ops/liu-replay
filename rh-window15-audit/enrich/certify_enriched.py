"""Full interval comparison, optimized over one added residual per parity."""
from pathlib import Path
from fractions import Fraction
import json,argparse,time
from flint import arb,arb_mat,ctx
from inputs import baseline_spectrum as base
from build_targeted_moments import proposals,load_intervals,save_intervals
OUT=Path(__file__).resolve().parent
def negative_witness(T,tag):
    S=base.load_ball(OUT/f'{tag}_schur.json');data=json.loads((OUT/f'{tag}_head_congruence.json').read_text())
    V=arb_mat([[arb(int(x))/(1<<data['bits']) for x in row] for row in data['coefficients']])
    diagonal=V.transpose()*S*V;j=next(i for i in range(64) if diagonal[i,i]<0)
    head=base.sub(V,0,64,j,j+1);R=base.sub(T,64,864,64,864);cross=base.sub(T,0,64,64,864)
    lower=-R.solve(cross.transpose()*head);bits=700;ints=[]
    for i in range(864):
        x=head[i,0] if i<64 else lower[i-64,0]
        ints.append(int((x.mid()*(1<<bits)+arb(1)/2).floor().unique_fmpz()))
    z=arb_mat([[arb(x)/(1<<bits)] for x in ints]);value=(z.transpose()*T*z)[0,0];assert value<0
    receipt={'bits':bits,'coefficients':[str(x) for x in ints],'raw_comparison_energy':value.str(65),'interpretation':'enlarged inverse-comparison direction; not a physical negative eigenvalue'}
    (OUT/f'{tag}_negative_witness.json').write_text(json.dumps(receipt,indent=2)+'\n')
    return value.str(65)
def main():
    pa=argparse.ArgumentParser();pa.add_argument('--baseline',type=Path,required=True);pa.add_argument('--parity',choices=['even','odd'],required=True);pa.add_argument('--thresholds',default='0,1e-95');pa.add_argument('--precision',type=int,default=1024)
    args=pa.parse_args();start=time.monotonic();ctx.prec=args.precision;ctx.threads=4;p=int(args.parity=='odd')
    assert json.loads((OUT/'moment_validation.json').read_text())['status']=='PASS_TARGETED_MOMENT_CROSSCHECKS'
    base.OUT=args.baseline;A=base.load_source('A',p);B=base.load_source('B',p);K=base.load_source('K',p);base.OUT=OUT
    n=864;gamma=[arb(3)/(4*i+2*p+1) for i in range(n)];Gamma=arb_mat(n,n)
    for i in range(n):Gamma[i,i]=gamma[i]
    m=arb(1)/2;a=arb_mat([[x] for x in proposals()[p]])
    Aa=A*a;Ba=B*a;b=arb_mat([[Aa[i,0]/gamma[i]] for i in range(n)])
    C3a=arb_mat([[x] for x in load_intervals(OUT/f'moment3_column_p{p}.json')]);scalars=load_intervals(OUT/f'scalar_moments_p{p}.json')
    r=Ba-A*b;M2v=C3a-B*b
    vv=(a.transpose()*Ba-Aa.transpose()*b)[0,0]
    vMv=scalars[1]-2*(b.transpose()*Ba)[0,0]+(b.transpose()*A*b)[0,0]
    vM2v=scalars[2]-2*(b.transpose()*C3a)[0,0]+(b.transpose()*B*b)[0,0]
    assert vv>0 and vMv>=m*vv and vM2v>=m*m*vv
    results=[];D=A-m*Gamma
    print(args.parity,'new residual scalars',[x.str(30) for x in [vv,vMv,vM2v]],flush=True)
    save_intervals([vv,vMv,vM2v],OUT/f'{args.parity}_residual_scalars.json',{'quantities':['norm_squared','M_energy','M2_energy'],'parity':p},require_narrow=False)
    for s_text in args.thresholds.split(','):
        s=Fraction(s_text);ideal=s==0;u=base.rat(s)+(arb(0) if ideal else arb(2)**-500);assert 0<=u<m
        H=B-(m+u)*A+m*u*Gamma;h=M2v-(m+u)*r;g=vM2v-(m+u)*vMv+m*u*vv
        rhs=arb_mat(n,n+1)
        for i in range(n):
            for j in range(n):rhs[i,j]=D[i,j]
            rhs[i,n]=h[i,0]
        solved=H.solve(rhs);Y=base.sub(solved,0,n,0,n);yh=base.sub(solved,0,n,n,n+1)
        sigma=g-(h.transpose()*yh)[0,0];assert sigma>0
        w=r-D*yh;numerator=base.symmetric(Gamma-D*Y-(w*w.transpose())/sigma)
        J=arb_mat(n,n)
        for i in range(n):
            for j in range(n):J[i,j]=numerator[i,j]/(gamma[i]*gamma[j])
        T=base.symmetric(J*K*J+(m-u)*J)
        tag=f'enriched_p{p}_s{s_text.replace("-","m").replace(".","p")}'
        print(tag,'augmented inverse positive; sigma',sigma.str(35),'seconds',round(time.monotonic()-start,2),flush=True)
        result,_=base.full_inertia(T,tag,want_witness=False)
        if ideal:assert result['inertia']['negative']<=1, 'enlargement cannot worsen ideal comparison inertia'
        if s==Fraction('1e-95'):
            assert result['inertia']['negative']<=(2 if p==0 else 1), 'enlargement cannot worsen the saved target comparison'
            if p==0:assert result['inertia']['negative']>=1, 'a positive target matrix would contradict the certified even trial'
        result.update({'s_label':s_text,'physical_threshold':None if ideal else str(s),'ideal_zero_shift':ideal,'projection_debit':'0' if ideal else '2^-500','inverse_trial_dimension':865,'constraint_dimension':864,'residual_schur_sigma':sigma.str(65),'tag':tag})
        if ideal and result['inertia']['negative']:result['negative_witness_energy']=negative_witness(T,tag)
        results.append(result)
        receipt={'status':'PASS_ENRICHED_COMPARISON_INERTIAS','L':'3/2','parity':args.parity,'arithmetic_precision_bits':ctx.prec,'results':results,'elapsed_seconds':round(time.monotonic()-start,3)}
        (OUT/f'{args.parity}_enriched_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    print('COMPLETE',args.parity,round(time.monotonic()-start,3),flush=True)
if __name__=='__main__':main()
