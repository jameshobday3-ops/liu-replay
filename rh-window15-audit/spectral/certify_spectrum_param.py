"""Full interval inertia of the optimized residual inverse comparison.

Floating-point factors/eigenvectors only propose exact dyadic congruences.
Acceptance uses Arb, strict diagonal dominance, and proven invertibility.
"""
from pathlib import Path
from fractions import Fraction
import argparse,json,gzip,time,hashlib,math
import numpy as np
import scipy.linalg as sla
import mpmath as mp
from flint import arb,arb_mat,ctx
OUT=Path(__file__).resolve().parent
def rat(x):
    x=Fraction(x);return arb(x.numerator)/x.denominator
def symmetric(A):return (A+A.transpose())/2
def identity(n):return arb_mat(n,n,[int(i==j) for i in range(n) for j in range(n)])
def sub(A,r0,r1,c0,c1):return arb_mat([[A[i,j] for j in range(c0,c1)] for i in range(r0,r1)])
def exact_float(x):
    num,den=float(x).as_integer_ratio();return arb(num)/den
def from_float(Q):return arb_mat([[exact_float(v) for v in row] for row in Q])
def as_float(A):return np.array([[float(A[i,j].mid()) for j in range(A.ncols())] for i in range(A.nrows())])
def save_ball(A,path,bits=700):
    scale=arb(2)**bits;rows=[]
    for i in range(A.nrows()):
        row=[]
        for j in range(A.ncols()):
            v=A[i,j]*scale
            row.append([str(int(v.lower().floor().unique_fmpz())),str(int(v.upper().ceil().unique_fmpz()))])
        rows.append(row)
    obj={'bits':bits,'nrows':A.nrows(),'ncols':A.ncols(),'rows':rows}
    path.write_text(json.dumps(obj,separators=(',',':'))+'\n')
def load_ball(path):
    obj=json.loads(path.read_text());scale=1<<(obj['bits']+1)
    return arb_mat([[arb(int(lo)+int(hi))/scale+arb(0,arb(int(hi)-int(lo))/scale) for lo,hi in row] for row in obj['rows']])
def load_source(name,p):
    path=OUT/f'{name}{p}.jsonl.gz';h=hashlib.sha256()
    with gzip.open(path,'rt',encoding='utf-8') as f:
        line=f.readline();h.update(line.encode());header=json.loads(line)
        if header.get('name')!=name or header.get('parity')!=p or header.get('bits')!=512 or header.get('radius_exponent')!=-510:
            raise ValueError(f'bad source header: {header}')
        n=int(header['n']); A=arb_mat(n,n)
        for i in range(n):
            line=f.readline();h.update(line.encode());row=json.loads(line);assert len(row)==864-i
            if len(row)!=n-i: raise ValueError((name,p,i,len(row),n-i))
            for j,value in enumerate(row,start=i):
                v=arb(int(value))/(1<<512)+arb(0,arb(2)**-510);A[i,j]=A[j,i]=v
        if f.read(): raise ValueError('trailing source data')
    receipt=json.loads((OUT/'source_build_receipt.json').read_text())
    assert h.hexdigest()==receipt['source_uncompressed_hashes'][path.name]
    return A
def dominance(A):
    positives=negatives=0;margin=None
    for i in range(A.nrows()):
        radius=sum((abs(A[i,j]) for j in range(A.ncols()) if i!=j),arb(0))
        if A[i,i]>radius:positives+=1;v=A[i,i]-radius
        elif -A[i,i]>radius:negatives+=1;v=-A[i,i]-radius
        else:raise ArithmeticError(f'row {i} not signed dominant: diagonal {A[i,i].str(10)}, radius {radius.str(10)}')
        margin=v if margin is None else min(margin,v)
    return {'positive':positives,'negative':negatives,'zero':0,'minimum_dominance_margin':margin.str(35)}
def finite_head(S,tag):
    n=S.nrows();mp.mp.dps=220
    raw=mp.matrix([[mp.mpf(S[i,j].mid().str(235,radius=False)) for j in range(n)] for i in range(n)])
    eig,Q=mp.eigsy(raw)
    bits=700;coeff=[[int(mp.nint(Q[i,j]*mp.mpf(2)**bits)) for j in range(n)] for i in range(n)]
    V=arb_mat([[arb(x)/(1<<bits) for x in row] for row in coeff])
    gram=dominance(V.transpose()*V)
    assert gram['positive']==n and gram['negative']==0
    C=symmetric(V.transpose()*S*V)
    powers=[]
    for i in range(n):
        assert eig[i]!=0
        powers.append(int(mp.floor(-mp.log(abs(eig[i]),2)/2)))
    D=arb_mat(n,n)
    for i,power in enumerate(powers):D[i,i]=arb(2)**power
    result=dominance(D*C*D)
    witness={'bits':bits,'coefficients':[[str(v) for v in row] for row in coeff],'scale_powers':powers,'gram_certificate':gram}
    (OUT/f'{tag}_head_congruence.json').write_text(json.dumps(witness,separators=(',',':'))+'\n')
    save_ball(S,OUT/f'{tag}_schur.json')
    negative_columns=[i for i in range(n) if C[i,i]<0]
    return result,V,negative_columns
def full_inertia(T,tag,want_witness=False):
    n=T.nrows();head=min(64,n//4);tail=n-head
    H=sub(T,0,head,0,head);cross=sub(T,0,head,head,n);R=sub(T,head,n,head,n)
    approx=as_float(R);approx=(approx+approx.T)/2
    chol=sla.cholesky(approx,lower=True,check_finite=True)
    Qf=sla.solve_triangular(chol.T,np.eye(tail),lower=False)
    assert np.count_nonzero(np.tril(Qf,-1))==0 and np.all(np.diag(Qf)!=0)
    np.savez_compressed(OUT/f'{tag}_tail_congruence.npz',Q=Qf)
    Q=from_float(Qf)
    transformed=symmetric(Q.transpose()*R*Q)
    tail_receipt=dominance(transformed)
    assert tail_receipt['positive']==tail and tail_receipt['negative']==0
    print(tag,'tail PASS',tail,flush=True)
    crossQ=cross*Q;solution=transformed.solve(crossQ.transpose())
    S=symmetric(H-crossQ*solution)
    head_receipt,V,negative_columns=finite_head(S,tag)
    print(tag,'head PASS',head_receipt,flush=True)
    result={'dimension':n,'head_dimension':head,'tail_dimension':tail,'tail':tail_receipt,'head':head_receipt,'inertia':{'positive':tail+head_receipt['positive'],'negative':head_receipt['negative'],'zero':0}}
    if want_witness:
        assert negative_columns
        j=negative_columns[0];v=sub(V,0,head,j,j+1);lower=-Q*solution*v
        bits=700;coeff=[]
        for i in range(n):
            x=v[i,0] if i<head else lower[i-head,0]
            coeff.append(int((x.mid()*(1<<bits)+arb(1)/2).floor().unique_fmpz()))
        z=arb_mat([[arb(x)/(1<<bits)] for x in coeff]);energy=(z.transpose()*T*z)[0,0]
        assert energy<0
        data={'bits':bits,'coefficients':[str(x) for x in coeff],'raw_comparison_energy':energy.str(65),'interpretation':'negative direction of comparison matrix ONLY; not a physical Rayleigh quotient'}
        (OUT/f'{tag}_negative_witness.json').write_text(json.dumps(data,indent=2)+'\n')
        result['negative_witness_energy']=energy.str(65)
        return result,z
    return result,None
def main():
    pa=argparse.ArgumentParser();pa.add_argument('--parity',choices=['even','odd'],required=True);pa.add_argument('--thresholds',default='0,1e-95');pa.add_argument('--precision',type=int,default=1024)
    args=pa.parse_args();ctx.prec=args.precision;ctx.threads=4;p=int(args.parity=='odd');start=time.monotonic()
    A=load_source('A',p);B=load_source('B',p);K=load_source('K',p);n=A.nrows();m=arb(1)/2
    gamma=[arb(3)/(4*i+2*p+1) for i in range(n)];Gamma=arb_mat(n,n)
    for i in range(n):Gamma[i,i]=gamma[i]
    D=A-m*Gamma;H=B-m*A
    Y=H.solve(D)
    print(args.parity,'inverse solve complete',round(time.monotonic()-start,2),flush=True)
    U=symmetric(Gamma-D*Y);V=symmetric(-Y.transpose()*D*Y)
    for i in range(n):
        for j in range(n):U[i,j]/=gamma[i]*gamma[j];V[i,j]/=gamma[i]*gamma[j]
    results=[]
    for text_s in args.thresholds.split(','):
        s=Fraction(text_s);ideal=s==0;u=rat(s)+(arb(0) if ideal else arb(2)**-500)
        assert 0<=u<m
        F=U+u*V;T=symmetric(F*K*F+(m-u)*F)
        tag=f'p{p}_s{text_s.replace("-","m").replace(".","p")}'
        result,z=full_inertia(T,tag,want_witness=ideal)
        result.update({'physical_threshold':str(s) if not ideal else None,'ideal_zero_shift':ideal,'s_label':text_s,'tag':tag,'projection_debit':'0' if ideal else '2^-500'})
        if ideal:
            a=Y*U*z/m;bits=700
            coeff=[str(int((a[i,0].mid()*(1<<bits)+arb(1)/2).floor().unique_fmpz())) for i in range(n)]
            (OUT/f'{tag}_residual_polynomial_proposal.json').write_text(json.dumps({'bits':bits,'raw_Legendre_coefficients':coeff,'status':'PROPOSAL_ONLY; residual moments not yet computed','parity':args.parity},indent=2)+'\n')
        results.append(result)
        receipt={'status':'PASS_COMPARISON_INERTIAS_NOT_GAP','L':'3/2','parity':args.parity,'precision_bits':ctx.prec,'source_radius':'2^-510','results':results,'elapsed_seconds':round(time.monotonic()-start,3)}
        (OUT/f'{args.parity}_spectral_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps(result,indent=2),flush=True)
    print('COMPLETE',args.parity,round(time.monotonic()-start,2),flush=True)
if __name__=='__main__':main()
