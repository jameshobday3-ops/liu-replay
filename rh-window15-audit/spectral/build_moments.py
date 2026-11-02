"""Integrated moments, rebuilt from the saved 28 October implementation.
Finite Binet remainder is included. This file does not certify positivity.
"""
import json,time,hashlib
from pathlib import Path
from flint import arb,acb,fmpq,ctx
OUT=Path(__file__).resolve().parent
def main():
    ctx.prec=12288;t=time.monotonic();count=3360;b=arb(384);s=arb(3073)/4
    shift=[arb(0) for _ in range(count)]
    for j in range(769):
        c=arb(4*j+1)/4;v=(b/c).atan()/b;row=[]
        for q in range(count):
            if q:v=c/(b*b*(2*q-1))-(c/b)**2*v
            if j<768:shift[q]+=v
            else:row.append(v)
        if j==768:js=row
    endpoint=[(acb(s,b)**(-r)).real for r in range(1536)]
    bern=[arb(fmpq.bernoulli(2*k))/(2*k) for k in range(1,769)]
    rem=arb(4)*arb.fac_ui(1537)/(arb(6)**1538*s**1538)
    assert rem<arb(2)**-4600
    rows=[];top=arb(2)**4096;maxwidth=0
    for q in range(count):
        d=2*q;h=js[q];correction=arb(0)
        for r in range(1,1536):
            h=(endpoint[r]-(d+1-r)*h)/(r*s)
            if (r+1)%2==0:correction+=bern[(r+1)//2-1]*h
        value=((s*s+b*b).log()/2-arb.pi().log()-arb(24)/5)/(d+1)-arb(1)/(d+1)**2+(s/(d+1)-arb(1)/2)*js[q]-correction-shift[q]
        value*=768/arb.pi()
        value+=arb(0,(768*rem/(arb.pi()*(d+1))).upper())
        assert value.is_finite()
        scaled=value*top
        lo=int(scaled.lower().floor().unique_fmpz());hi=int(scaled.upper().ceil().unique_fmpz())
        assert 0<=hi-lo<=3,(q,hi-lo)
        rows.append([str(lo),str(hi)]);maxwidth=max(maxwidth,hi-lo)
        if q%512==0:print('moment',q,'width',hi-lo,'seconds',round(time.monotonic()-t,2),flush=True)
    obj={'status':'INTERVAL_MOMENTS_ONLY_NOT_POSITIVITY','precision_bits':ctx.prec,'output_bits':4096,'count':count,'omega':768,'comparison_constant':'24/5','binet_remainder_included':True,'max_width_units':maxwidth,'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT/'moments.json').write_text(json.dumps(obj,separators=(',',':'))+'\n')
    print('MOMENTS PASS',count,'seconds',round(time.monotonic()-t,3),flush=True)
if __name__=='__main__':main()
