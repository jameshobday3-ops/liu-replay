"""Rigorous targeted M^3 columns and M^4 scalars for exact dyadic polynomials.
M^k uses support-restricted shift paths. Columns use polynomial Hankel
products; M^4 is checked by a separate piecewise M^2 f norm calculation.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from collections import defaultdict
from functools import lru_cache
import json,time,argparse,gc
from flint import arb,arb_poly,ctx
OUT=Path(__file__).resolve().parent
POWERS=[(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11),(13,13),(16,2),(17,17),(19,19)]
DEG=1727;N=864
def rational(q):q=Q(q);return arb(q.numerator)/q.denominator
def lograt(q):return rational(q).log()
def proposals():
    result=[]
    for p in [0,1]:
        raw=json.loads((OUT/f'inputs/p{p}_s0_residual_polynomial_proposal.json').read_text())
        coeff=[arb(int(v))/(1<<raw['bits']) for v in raw['raw_Legendre_coefficients']]
        assert len(coeff)==N and raw['bits']==700
        result.append(coeff)
    return result
def polynomial_pair(coeff):
    x=arb_poly([0,1]);previous=arb_poly([1]);current=x
    f=[previous*coeff[0][0],current*coeff[1][0]]
    for degree in range(2,DEG+1):
        previous,current=current,((2*degree-1)*x*current-(degree-1)*previous)*(arb(1)/degree)
        f[degree%2]+=current*coeff[degree%2][degree//2]
    assert f[0].degree()==1726 and f[1].degree()==1727
    return f
def paths(k):
    shifts=[(Q(1),arb(24)/5)]
    for n,p in POWERS:
        w=-arb(p).log()/arb(n).sqrt();shifts.extend([(Q(n),w),(Q(1,n),w)])
    groups=defaultdict(dict);active=0
    for path in product(shifts,repeat=k):
        q=Q(1);partial=[q];weight=arb(1)
        for r,w in path:q*=r;partial.append(q);weight*=w
        low=min(partial);high=max(partial);length=arb(3)-lograt(high/low)
        if length<0:continue
        assert length>0;active+=1
        for key,w in [((1,high),weight),((-1,low),-weight)]:
            groups[q][key]=groups[q].get(key,arb(0))+w
    return groups,active
def endpoint(key):side,q=key;return arb(side)-lograt(q)*arb(2)/3
def column_pair(k):
    start=time.monotonic();ctx.prec=8192;coeff=proposals();f=polynomial_pair(coeff)
    groups,active=paths(k);selected=sorted(q for q in groups if q>=1)
    @lru_cache(maxsize=96)
    def endpoint_primitive(key):
        x=endpoint(key);v=x;values=[]
        for j in range(2*DEG+1):values.append(v/(j+1));v*=x
        return arb_poly(values)
    moments=[[arb(0) for _ in range(DEG+1)] for _ in [0,1]]
    for ii,q in enumerate(selected):
        boundary=arb_poly([])
        for key,w in groups[q].items():boundary+=endpoint_primitive(key)*w
        displacement=lograt(q)*arb(2)/3
        for p in [0,1]:
            shifted=f[p](arb_poly([displacement,1]));degree=shifted.degree()
            convolution=arb_poly(list(reversed(shifted.coeffs())))*boundary
            fac=arb(3)/2*(1 if q==1 else 2)
            for j in range(p,DEG+1,2):moments[p][j]+=fac*convolution[degree+j]
        if ii%16==0:print('column',k,ii+1,'/',len(selected),'seconds',round(time.monotonic()-start,2),flush=True)
    columns=[[],[]];x=arb_poly([0,1]);previous=arb_poly([1]);current=x
    for degree in range(DEG+1):
        if degree==0:P=previous
        elif degree==1:P=current
        else:
            previous,current=current,((2*degree-1)*x*current-(degree-1)*previous)*(arb(1)/degree);P=current
        p=degree%2
        value=sum((P[j]*moments[p][j] for j in range(p,degree+1,2)),arb(0))
        assert value.is_finite();columns[p].append(value)
    metadata={'power':k,'precision_bits':ctx.prec,'ordered_paths':25**k,'nonempty_paths':active,'total_distinct_displacements':len(groups),'displacements_computed_using_parity':len(selected),'elapsed_seconds':round(time.monotonic()-start,3)}
    for p in [0,1]:save_intervals(columns[p],OUT/f'moment{k}_column_p{p}.json',metadata|{'parity':p})
    endpoint_primitive.cache_clear();print('COLUMN PASS',metadata,flush=True);gc.collect()
    return columns
def save_intervals(values,path,metadata,require_narrow=True):
    bits=700;rows=[];maxwidth=0
    for i,v in enumerate(values):
        scaled=v*(1<<bits);lo=int(scaled.lower().floor().unique_fmpz());hi=int(scaled.upper().ceil().unique_fmpz())
        assert hi>=lo
        if require_narrow:assert hi-lo<=3,(path.name,i,hi-lo,v.rad())
        rows.append([str(lo),str(hi)]);maxwidth=max(maxwidth,hi-lo)
    path.write_text(json.dumps(metadata|{'status':'PASS_TARGETED_MOMENT_INTERVALS','bits':bits,'max_width_units':maxwidth,'rows':rows},separators=(',',':'))+'\n')
def load_intervals(path):
    r=json.loads(path.read_text());den=1<<(r['bits']+1)
    return [arb(int(lo)+int(hi))/den+arb(0,arb(int(hi)-int(lo))/den) for lo,hi in r['rows']]
def scalar_pair():
    start=time.monotonic();ctx.prec=12288;coeff=proposals();f=polynomial_pair(coeff)
    events={};groups_by_k={};counts={}
    for k in [1,2]:
        groups,active=paths(k);groups_by_k[k]=groups;counts[k]=active
        for q in groups:
            for key in groups[q]:
                if key not in events:events[key]=[[arb_poly([]),arb_poly([])] for _ in [0,1]]
    keys=sorted(events,key=lambda key:float(endpoint(key)))
    assert all(endpoint(keys[i])<endpoint(keys[i+1]) for i in range(len(keys)-1))
    assert endpoint(keys[0])==-1 and endpoint(keys[-1])==1
    allq=sorted({q for g in groups_by_k.values() for q in g if q>=1})
    for ii,q in enumerate(allq):
        disp=lograt(q)*arb(2)/3;positive=[f[p](arb_poly([disp,1])) for p in [0,1]]
        for qq in ([q] if q==1 else [q,1/q]):
            for p in [0,1]:
                shifted=positive[p] if qq==q else arb_poly([v*((-1)**(p+j)) for j,v in enumerate(positive[p].coeffs())])
                for k in [1,2]:
                    if qq not in groups_by_k[k]:continue
                    for key,w in groups_by_k[k][qq].items():events[key][k-1][p]-=shifted*w
        if ii%16==0:print('scalar shifts',ii+1,'/',len(allq),'seconds',round(time.monotonic()-start,2),flush=True)
    running=[[arb_poly([]),arb_poly([])] for _ in [0,1]];sums=[[arb(0),arb(0),arb(0)] for _ in [0,1]]
    for i,key in enumerate(keys[:-1]):
        for k in [0,1]:
            for p in [0,1]:running[k][p]+=events[key][k][p]
        left=endpoint(key);right=endpoint(keys[i+1])
        if right<=0:continue
        if left<0:left=arb(0)
        for p in [0,1]:
            h=running[0][p];g=running[1][p]
            for j,integrand in enumerate([h*h,h*g,g*g]):
                primitive=integrand.integral();sums[p][j]+=3*(primitive(right)-primitive(left))
        if i%16==0:print('scalar pieces',i+1,'/',len(keys)-1,'seconds',round(time.monotonic()-start,2),flush=True)
    for p in [0,1]:
        assert sums[p][0]>0 and sums[p][2]>0
        save_intervals(sums[p],OUT/f'scalar_moments_p{p}.json',{'parity':p,'quantities':['a_star_B_a','a_star_C3_a','a_star_C4_a'],'method':'piecewise polynomial Mf and M^2f; half-interval by parity; antiderivatives with Arb','precision_bits':ctx.prec,'segments':len(keys)-1,'nonempty_M2_paths':counts[2],'elapsed_seconds':round(time.monotonic()-start,3)})
        print('SCALAR PASS',p,[v.str(25) for v in sums[p]],flush=True)
def main():
    pa=argparse.ArgumentParser();pa.add_argument('--stage',choices=['column2','column3','scalars','all'],default='all');args=pa.parse_args()
    if args.stage in ['all','column2']:column_pair(2)
    if args.stage in ['all','column3']:column_pair(3)
    if args.stage in ['all','scalars']:scalar_pair()
if __name__=='__main__':main()
