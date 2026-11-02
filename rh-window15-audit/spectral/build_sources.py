"""Build all raw Legendre moment entries with Arb and analytic kernel errors.
The source recurrences are inherited from the saved 28 October verifier.
M^2 includes all intermediate support restrictions and its boundary term.
"""
from pathlib import Path
import json,time,hashlib,gzip
from flint import arb,ctx
from physical_even_trial import prime_powers
from fractions import Fraction
OUT=Path(__file__).resolve().parent
def main():
    ctx.prec=8192;start=time.monotonic();L=arb(3)/2;N=1728;top=2*N-2;last=top+1
    mom=json.loads((OUT/'moments.json').read_text())
    assert mom['count']==3360 and mom['output_bits']==4096 and mom['omega']==768 and mom['comparison_constant']=='24/5'
    assert all(0<=int(hi)-int(lo)<=1 for lo,hi in mom['rows'])
    budget=json.loads((OUT/'budget_receipt.json').read_text());assert budget['status']=='PASS_SOURCE_AND_KERNEL_BUDGET'
    shifts=[(arb(0),arb(24)/5,'0')]
    powers=prime_powers(Fraction(3,2))
    assert [n for n,p in powers]==[2,3,4,5,7,8,9,11,13,16,17,19]
    for n,p in powers:
        d=arb(n).log();w=-arb(p).log()/arb(n).sqrt();shifts.extend([(d,w,str(n)),(-d,w,'-'+str(n))])
    shifts.sort(key=lambda x:float(x[0]));zero=12
    assert shifts[zero][2]=='0' and len(shifts)==25
    assert all(shifts[i][0]<shifts[i+1][0] for i in range(24))
    cache={}
    def leg_at(key):
        if key not in cache:
            sign,ib,edge=key;x=(sign*L+shifts[ib][0]-shifts[edge][0])/L;v=[arb(1),x]
            for degree in range(2,last+1):v.append(((2*degree-1)*x*v[-1]-(degree-1)*v[-2])/degree)
            cache[key]=v
        return cache[key]
    a0=[arb(0) for _ in range(last+1)];a1=list(a0);b0=list(a0);b1=list(a0);pairs=0;active=0
    for ia,(da,wa,_) in enumerate(shifts):
        for ib,(db,wb,_) in enumerate(shifts):
            pairs+=1;low=min(zero,ia,ib);high=max(zero,ia,ib);left=-L-shifts[low][0];right=L-shifts[high][0]
            if right<left:continue
            assert right>left;active+=1
            vl=leg_at((-1,ib,low));vr=leg_at((1,ib,high))
            row=[right-left]+[L*(vr[j+1]-vr[j-1]-vl[j+1]+vl[j-1])/(2*j+1) for j in range(1,top+1)]
            for j in range(0,top+1,2):
                b0[j]+=wa*wb*row[j]
                if ia==zero:a0[j]+=wb*row[j]
            for j in range(1,top,2):
                value=((j+1)*row[j+1]+j*row[j-1])/(2*j+1)+(da-db)*row[j]/L
                b1[j]+=wa*wb*value
                if ia==zero:a1[j]+=wb*value
    assert pairs==625 and active==379
    endpoint=[arb(0) for _ in range(top+1)]
    for k,(d,w,_) in enumerate(shifts[:zero]):
        vals=leg_at((1,k,zero))
        for j in range(1,top):endpoint[j]+=w*(vals[j+1]-vals[j-1])
    print('prime seeds complete',round(time.monotonic()-start,3),flush=True)
    checked={};max_errors={};hashes={}
    def write_matrix(first,second,name):
        previous=current=None;count=0;maxerr=arb(0);handles=[];hashers=[]
        for p in [0,1]:
            handles.append(gzip.open(OUT/f'{name}{p}.jsonl.gz','wt',encoding='utf-8'))
            hashers.append(hashlib.sha256())
            line=json.dumps({'name':name,'parity':p,'n':864,'bits':512,'radius_exponent':-510},separators=(',',':'))+'\n'
            handles[-1].write(line);hashers[-1].update(line.encode())
        for i in range(N):
            if i==0:row={j:first[j] for j in range(0,top+1,2)}
            elif i==1:row={j:second[j] for j in range(1,top,2)}
            else:
                row={}
                for j in range(i,top-i+1,2):
                    value=previous[j]+arb(2*i-1)*(current[j-1]-current[j+1])/(2*j+1)
                    if name=='B':value+=2*L*endpoint[i-1]*endpoint[j]/(2*j+1)
                    row[j]=value
            p=i%2;centers=[]
            for j in range(i,N,2):
                value=row[j]*(L*L if name=='K' else 1)
                center=int((value.mid()*(1<<512)+arb(1)/2).floor().unique_fmpz())
                err=(value-arb(center)/(1<<512)).abs_upper()
                assert err<arb(2)**-511,(name,i,j,err)
                centers.append(str(center));count+=1;maxerr=max(maxerr,err)
            line=json.dumps(centers,separators=(',',':'))+'\n';handles[p].write(line);hashers[p].update(line.encode())
            previous,current=current,row
        assert count==N*(N+2)//4
        for p in [0,1]:handles[p].close();hashes[f'{name}{p}.jsonl.gz']=hashers[p].hexdigest()
        checked[name]=count;max_errors[name]=maxerr.str(30)
        print('PASS source entries',name,count,round(time.monotonic()-start,3),flush=True)
    write_matrix(a0,a1,'A');write_matrix(b0,b1,'B')
    moments=[arb(int(lo)+int(hi))/(1<<4097) for lo,hi in mom['rows']]
    k0=[arb(0) for _ in range(top+1)];k1=list(k0);fac=arb(1);lp=arb(1);bp=arb(1)
    for q in range(3360):
        degree=2*q
        if q:fac*=degree*(degree-1);lp*=L*L;bp*=(1536*L)**2
        coeff=(2*lp+(-1 if q%2 else 1)*bp*moments[q])/fac
        f=arb(8)/((degree+1)*(degree+2))
        for j in range(0,min(degree,top)+1,2):
            k0[j]+=coeff*f;f*=arb((degree+1-j)*(degree-j))/((degree+j+3)*(degree+j+4))
        if degree:
            v=-arb(8)*degree/((degree+1)*(degree+2)*(degree+4))
            for j in range(1,min(degree,top+1),2):
                k1[j]+=coeff*v;v*=arb((degree-j-1)*(degree-j+2))/((degree+j+2)*(degree+j+5))
    write_matrix(k0,k1,'K')
    rec={'status':'PASS_ALL_2242080_SOURCE_ENTRIES','L':'3/2','dimension':N,'precision_bits':ctx.prec,'ordered_shift_pairs':pairs,'nonempty_shift_pairs':active,'entries_checked':checked,'max_raw_arithmetic_errors':max_errors,'raw_ball_radius':'2^-510','raw_arithmetic_error_upper':'2^-511','kernel_additional_raw_error_upper':'2^-650','moments_sha256':hashlib.sha256((OUT/'moments.json').read_bytes()).hexdigest(),'source_uncompressed_hashes':hashes,'elapsed_seconds':round(time.monotonic()-start,3)}
    (OUT/'source_build_receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2),flush=True)
if __name__=='__main__':main()
